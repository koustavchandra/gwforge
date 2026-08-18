"""Which parameters the Fisher matrix is taken over, and how far to step them.

The named sets exist because choosing the wrong one is a *silent* mistake rather
than a loud one -- ``convert_to_lal_binary_black_hole_parameters`` will happily
turn ``chi_1``/``chi_2`` into zero-tilt spins for a precessing waveform, so the
forecast comes out looking fine and describes the wrong binary.
"""

import numpy

#: Parameters whose derivatives are exact because the waveform does not depend
#: on them. ``luminosity_distance`` is an overall scaling; the rest enter only
#: through the detector projection.
ANALYTIC_PARAMETERS = (
    "luminosity_distance",
    "geocent_time",
    "psi",
    "ra",
    "dec",
)

#: Fisher parameters used when the configuration does not name any. These are
#: the eleven that map one-to-one onto GWFast's ``ParNums``, which keeps the
#: cross-validation a pure Jacobian rotation.
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

#: Parameters for a precessing waveform such as ``SEOBNRv5PHM`` or
#: ``IMRPhenomXPHM``.
#:
#: Using :data:`DEFAULT_PARAMETERS` for a precessing model is a silent mistake
#: rather than a loud one: ``convert_to_lal_binary_black_hole_parameters`` turns
#: ``chi_1``/``chi_2`` into ``a_1``/``a_2`` with both tilts set to zero, so the
#: waveform generates happily and forecasts an aligned-spin binary.
#:
#: Note the Fisher is singular at exactly zero tilt, where ``phi_12`` and
#: ``phi_jl`` do not affect the waveform at all.
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

#: Parameters for a precessing *and* eccentric waveform such as ``pyEFPEHM``.
#: As above, the Fisher is singular at exactly zero eccentricity, where
#: ``mean_anomaly`` has no effect.
ECCENTRIC_PRECESSING_PARAMETERS = PRECESSING_PARAMETERS + (
    "eccentricity",
    "mean_anomaly",
)

#: Hard physical limits. Steps are shrunk so that ``value +/- 2 * step`` stays
#: strictly inside these, which matters most for ``symmetric_mass_ratio``: a
#: source near equal mass sits right against 0.25, and stepping past it makes
#: bilby's eta-to-q inversion return NaN.
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

#: Starting guesses for the step search. Angles and spins get an absolute step
#: (they are O(1) and may pass through zero); masses and distances get a
#: relative one. :func:`calibrate_step` rescales these per source, so they only
#: need to be within an order of magnitude or two.
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
    # anything else: for pyEFPEHM a step of 1e-3 changes the waveform by 19%,
    # 500 times more than the same step in mean_anomaly. Starting small keeps
    # the calibration from having to climb down several decades.
    "eccentricity": 1e-5,
    "mean_anomaly": 1e-4,
}

#: Parameters whose step size is absolute rather than a fraction of the value.
#: An angle near zero has no meaningful relative step, and a spin may be
#: identically zero.
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
