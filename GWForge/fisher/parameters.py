import numpy

# Parameters whose derivatives can be obtained analytically
ANALYTIC_PARAMETERS = (
    "luminosity_distance",
    "geocent_time",
    "psi",
    "ra",
    "dec",
)

# Fisher parameters used when the configuration does not name any.
DEFAULT_PARAMETERS = (
    "chirp_mass",
    "symmetric_mass_ratio",
    "luminosity_distance",
    "ra",
    "dec",
    "theta_jn",
    "psi",
    "geocent_time",
    "phase",
    "chi_1",
    "chi_2",
)

# Parameters for a precessing waveform
# Note the Fisher is singular at exactly zero tilt,
# Same is true for say inclination for π/2
PRECESSING_PARAMETERS = (
    "chirp_mass",
    "symmetric_mass_ratio",
    "luminosity_distance",
    "ra",
    "dec",
    "theta_jn",
    "psi",
    "geocent_time",
    "phase",
    "a_1",
    "a_2",
    "tilt_1",
    "tilt_2",
    "phi_12",
    "phi_jl",
)

# Parameters for a precessing *and* eccentric waveform
ECCENTRIC_PRECESSING_PARAMETERS = PRECESSING_PARAMETERS + (
    "eccentricity",
    "mean_anomaly",
)

# Hard physical limits. Steps are shrunk so that ``value +/- 2 * step`` stays
# strictly inside these, which matters most for ``symmetric_mass_ratio``: a
# source near equal mass sits right against 0.25, and stepping past it makes
# bilby's eta-to-q inversion return NaN.
_PARAMETER_BOUNDS = {
    "symmetric_mass_ratio": (0.0, 0.25),
    "mass_ratio": (0.0, 1.0),
    "chi_1": (-1.0, 1.0),
    "chi_2": (-1.0, 1.0),
    "a_1": (0.0, 1.0),
    "a_2": (0.0, 1.0),
    "tilt_1": (0.0, numpy.pi),
    "tilt_2": (0.0, numpy.pi),
    "eccentricity": (0.0, 1.0),
    "chirp_mass": (0.0, None),
    "total_mass": (0.0, None),
    "mass_1": (0.0, None),
    "mass_2": (0.0, None),
    "luminosity_distance": (0.0, None),
    "lambda_1": (0.0, None),
    "lambda_2": (0.0, None),
}

# Starting guesses for the step search. Angles and spins get an absolute step
# (they are O(1) and may pass through zero); masses and distances get a
# relative one. :func:`calibrate_step` rescales these per source, so they only
# need to be within an order of magnitude or two.
DEFAULT_STEP_SIZES = {
    "chirp_mass": 1e-6,
    "symmetric_mass_ratio": 1e-6,
    "mass_1": 1e-6,
    "mass_2": 1e-6,
    "mass_ratio": 1e-6,
    "total_mass": 1e-6,
    "chi_1": 1e-5,
    "chi_2": 1e-5,
    "a_1": 1e-5,
    "a_2": 1e-5,
    "tilt_1": 1e-5,
    "tilt_2": 1e-5,
    "phi_12": 1e-5,
    "phi_jl": 1e-5,
    "theta_jn": 1e-5,
    "phase": 1e-5,
    "lambda_1": 1e-3,
    "lambda_2": 1e-3,
    # An eccentric waveform responds far more sharply to eccentricity than to
    # anything else
    "eccentricity": 1e-5,
    "mean_anomaly": 1e-4,
}

# Parameters whose step size is absolute rather than a fraction of the value.
# An angle near zero has no meaningful relative step, and a spin may be
# identically zero.
_ABSOLUTE_STEP_PARAMETERS = frozenset(
    {
        "chi_1",
        "chi_2",
        "a_1",
        "a_2",
        "tilt_1",
        "tilt_2",
        "phi_12",
        "phi_jl",
        "theta_jn",
        "phase",
        "psi",
        "ra",
        "dec",
        "eccentricity",
        "mean_anomaly",
    }
)


def step_size(name, value, step_sizes=None):
    """Starting finite-difference step for one parameter.

    Parameters
    ----------
    name : str
        Parameter name.
    value : float
        Fiducial value.
    step_sizes : dict or None
        Overrides for :data:`DEFAULT_STEP_SIZES`.

    Returns
    -------
    float
    """
    sizes = dict(DEFAULT_STEP_SIZES)
    if step_sizes:
        sizes.update(step_sizes)
    scale = sizes.get(name, 1e-6)
    if name in _ABSOLUTE_STEP_PARAMETERS or value == 0.0:
        return bounded_step(name, value, scale)
    return bounded_step(name, value, scale * abs(value))


def bounded_step(name, value, step):
    """Shrink a step so that ``value +/- 2 * step`` stays physical.

    The five-point stencil reaches twice the step either side, so a source close
    to a boundary -- an almost equal-mass binary against ``eta = 0.25``, an
    extremal spin against 1 -- would otherwise be evaluated at parameters the
    waveform cannot represent.

    A value sitting exactly *on* a bound is physical, not an error -- a
    non-spinning binary has ``a_1 = 0`` and a circular one has
    ``eccentricity = 0``. There is simply no two-sided derivative there, and the
    Fisher matrix is singular in that parameter because the waveform cannot
    respond to a displacement in the forbidden direction. That gets its own
    message, because "``a_1 = 0.0`` is outside its physical range" is a
    thoroughly confusing way to describe a perfectly ordinary injection.

    Parameters
    ----------
    name : str
        Parameter name; unbounded names are returned unchanged.
    value : float
        Fiducial value.
    step : float
        Requested step.

    Returns
    -------
    float
        The step, possibly reduced. Never returns zero for a positive input.

    Raises
    ------
    ValueError
        If the value is on or outside a bound, with different wording for the
        two cases.
    """
    bounds = _PARAMETER_BOUNDS.get(name)
    if bounds is None:
        return step
    lower, upper = bounds
    room = []
    if lower is not None:
        room.append(value - lower)
    if upper is not None:
        room.append(upper - value)
    if not room:
        return step

    if min(room) < 0.0:
        raise ValueError(
            "{} = {} is outside its physical range {}".format(name, value, bounds)
        )
    if min(room) == 0.0:
        raise ValueError(
            "{} = {} sits exactly on the edge of its physical range {}, so it "
            "cannot be differentiated and the Fisher matrix would be singular "
            "in it. Either move the injection off the boundary or drop {} from "
            "the Fisher parameters.".format(name, value, bounds, name)
        )
    # Keep a 20% margin so the outermost stencil point is comfortably inside.
    return min(step, 0.4 * min(room))


REDUNDANT_PARAMETER_GROUPS = (
    (
        "mass_1",
        "mass_2",
        "chirp_mass",
        "total_mass",
        "symmetric_mass_ratio",
        "mass_ratio",
        "mass_1_source",
        "mass_2_source",
        "chirp_mass_source",
        "total_mass_source",
    ),
    ("chi_1", "a_1", "tilt_1", "cos_tilt_1", "spin_1z", "chi_1_in_plane"),
    ("chi_2", "a_2", "tilt_2", "cos_tilt_2", "spin_2z", "chi_2_in_plane"),
)


def strip_shadowed_parameters(parameters, fisher_parameters):
    """Remove source parameters that would silently pin a Fisher parameter.

    Parameters
    ----------
    parameters : dict
        Source parameters as read from the configuration or an injection file.
    fisher_parameters : sequence of str
        The names the Fisher matrix is taken over.

    Returns
    -------
    tuple
        ``(parameters, removed)``. ``parameters`` is a new dict; ``removed`` is
        the sorted list of keys dropped, which the caller should log so the
        change is visible rather than silent.

    """
    varied = set(fisher_parameters)
    parameters = dict(parameters)
    removed = []
    for group in REDUNDANT_PARAMETER_GROUPS:
        if not varied.intersection(group):
            continue
        for name in group:
            if name in parameters and name not in varied:
                del parameters[name]
                removed.append(name)
    return parameters, sorted(removed)
