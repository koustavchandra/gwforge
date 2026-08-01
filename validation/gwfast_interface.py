"""Shared plumbing for cross-validating GWForge's Fisher matrix against GWFast.

Two independent codes only agree if they are asked the same question, so most of
this module is about removing every avoidable difference between them: the same
detector geometry, the same PSD files, the same frequency grid, the same band,
and the static long-wavelength response on both sides.

Conventions
-----------
The mapping was established by ``tests/test_antenna_response.py`` and is reused
rather than rediscovered:

* ``det_xax = xarm_azimuth + 45`` -- GWFast orients a detector by the bisector
  of its arms, bilby by the x arm.
* ``theta = pi/2 - dec`` and ``phi = ra``.
* ``dL`` is in Gpc, GWForge's ``luminosity_distance`` in Mpc.
* GWFast writes :math:`h \\sim e^{+i\\Psi}` where bilby writes :math:`e^{-i\\Psi}`.
  That is a global conjugation, which cancels in
  :math:`\\langle \\partial h \\mid \\partial h \\rangle`, so the Fisher matrix is
  unaffected.

Detector choice
---------------
The network is three L-shaped sites at zero elevation. That is not cosmetic:
GWFast models a spherical Earth and ignores elevation, and it collapses a
triangle's three arms onto a single vertex where bilby walks them apart -- so ET
cannot be compared cleanly and is deliberately excluded. Three sites are needed
because the sky-position block of a two-detector Fisher is singular.
"""

import os

import numpy

#: GWForge's bundled noise curves.
NOISE_CURVES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "GWForge",
    "ifo",
    "noise_curves",
)

#: ``(name, latitude, longitude, x-arm azimuth, arm length in km, ASD file)``.
#: The y arm is always the x arm plus 90 degrees, because GWFast assumes
#: perpendicular arms.
DETECTORS = (
    ("CE40", 46.0, -125.0, 260.0, 40.0, "CE40-asd.txt"),
    ("CE20", 29.0, -94.0, 200.0, 20.0, "CE20-asd.txt"),
    ("CE20AU", -31.34, 115.91, 2.0, 20.0, "CE20-asd.txt"),
)

#: GWForge parameter name, GWFast parameter name, and
#: ``d(GWFast) / d(GWForge)``. Every GWFast parameter depends on exactly one
#: GWForge parameter, so the Fisher change of basis is a diagonal congruence.
#:
#: The four intrinsic parameters carry a minus sign, and it is a relative sign,
#: not a global one. GWFast takes ``hp, hc`` from LAL -- the same convention
#: bilby uses -- and then multiplies by
#: ``exp(+1j * (phiD + phiL + 2 pi f tcoal - Phicoal))``, where bilby applies
#: ``exp(-2j pi f tau)``. So GWFast's extrinsic phase is conjugated relative to
#: bilby's while the intrinsic LAL phase is not, and derivatives with respect to
#: intrinsic and extrinsic parameters acquire opposite signs. Flipping the four
#: intrinsic parameters fixes that (flipping the seven extrinsic ones instead
#: would do equally well -- an overall sign cancels in a congruence).
#:
#: This was measured, not assumed, and it only affects correlations: a congruence
#: by a diagonal of plus and minus ones leaves the diagonal untouched, which is
#: why the sigmas agreed to 0.3% before the signs were corrected. Candidate flip
#: sets were scored by the largest disagreement in the normalised correlation
#: matrix, over three sources and both approximants:
#:
#: =========================  =========================
#: flip set                   worst correlation residual
#: =========================  =========================
#: none                       1.46 - 1.97
#: extrinsic angles only      0.32 - 0.92
#: intrinsic (this one)       0.035 - 0.094
#: =========================  =========================
#:
#: ``Phicoal`` additionally carries a magnitude of two, because bilby's
#: quadrupole picks up :math:`e^{-2i\phi}` against GWFast's single global
#: factor; with a non-higher-mode approximant the ratio of the two ``phase``
#: Fisher diagonals is 4.0014, i.e. exactly :math:`2^2`.
#:
#: ``tcoal`` has unit magnitude because GWFast already divides that row of the
#: Fisher by 86400, leaving it as a derivative with respect to a coalescence
#: time in seconds.
PARAMETER_MAP = (
    ("chirp_mass", "Mc", -1.0),
    ("symmetric_mass_ratio", "eta", -1.0),
    ("luminosity_distance", "dL", 1e-3),
    ("ra", "phi", 1.0),
    ("dec", "theta", -1.0),
    ("theta_jn", "iota", 1.0),
    ("psi", "psi", 1.0),
    ("geocent_time", "tcoal", 1.0),
    ("phase", "Phicoal", 2.0),
    ("chi_1", "chi1z", -1.0),
    ("chi_2", "chi2z", -1.0),
)

#: Base step for GWFast's ``numdifftools`` derivatives.
#:
#: GWFast defaults to ``1e-5``, which sits below the noise floor of a LAL
#: waveform and gives a Fisher matrix that has not converged: with the default,
#: GWFast's ``iota`` Fisher entry for IMRPhenomXHM disagrees with GWForge by a
#: factor of 2.05. Raising the step to 1e-4 or beyond puts GWFast on the same
#: plateau GWForge calibrates onto, and every numerically differentiated
#: parameter then agrees to 1-2%:
#:
#: ===========  ==========  ==========  ==========  ==========
#: base_step    Mc          eta         iota        chi2z
#: ===========  ==========  ==========  ==========  ==========
#: 1e-5         1.056       1.014       2.048       1.018
#: 1e-4         1.019       1.012       0.981       1.014
#: 1e-3         1.019       1.012       0.980       1.014
#: 1e-2         1.019       1.012       0.980       1.013
#: ===========  ==========  ==========  ==========  ==========
#:
#: This is the same plateau GWForge finds; see
#: :data:`GWForge.fisher.derivatives.TARGET_FRACTIONAL_CHANGE`.
GWFAST_BASE_STEP = 1e-3

#: The GWForge parameters compared, in matrix order.
PARAMETERS = tuple(name for name, _, _ in PARAMETER_MAP)


def patch_gwfast():
    """Make GWFast 1.1.2 importable and usable on a modern scientific stack.

    Two things are broken out of the box and neither is GWForge's fault:

    * ``scipy.integrate.cumtrapz`` was removed in SciPy 1.14; GWFast's
      ``GWSignal`` constructor still calls it.
    * ``numdifftools`` cannot handle the complex derivatives a waveform
      produces. GWFast ships a patch for this
      (``gwfast/.patch/patch_ndt_complex_0-9-41.patch``) that the user is
      expected to apply to their site-packages by hand; every LAL waveform
      crashes without it, and LAL waveforms are exactly what we need for
      IMRPhenomXHM.

    Both are applied in memory so nothing installed is modified.
    """
    import scipy.integrate

    if not hasattr(scipy.integrate, "cumtrapz"):
        scipy.integrate.cumtrapz = scipy.integrate.cumulative_trapezoid

    from numdifftools.limits import _Limit

    if getattr(_Limit._add_error_to_outliers, "_gwforge_patched", False):
        return

    original = _Limit._add_error_to_outliers

    def add_error_to_outliers(der, trim_fact=10):
        if numpy.iscomplexobj(der):
            return numpy.sqrt(
                original(numpy.real(der), trim_fact) ** 2
                + original(numpy.imag(der), trim_fact) ** 2
            )
        return original(der, trim_fact)

    add_error_to_outliers._gwforge_patched = True
    _Limit._add_error_to_outliers = staticmethod(add_error_to_outliers)


def bilby_network(minimum_frequency, maximum_frequency):
    """Build the validation network as bilby interferometers.

    Constructed explicitly rather than through
    :class:`GWForge.ifo.detectors.Network` so the geometry handed to bilby and
    the geometry handed to GWFast come from the same literal numbers.

    Parameters
    ----------
    minimum_frequency, maximum_frequency : float
        Band limits in Hz, applied identically on both sides.

    Returns
    -------
    bilby.gw.detector.InterferometerList
    """
    import bilby
    from bilby.gw.detector.psd import PowerSpectralDensity

    interferometers = []
    for name, latitude, longitude, xarm, length, asd in DETECTORS:
        interferometers.append(
            bilby.gw.detector.Interferometer(
                name=name,
                power_spectral_density=PowerSpectralDensity.from_amplitude_spectral_density_file(
                    os.path.join(NOISE_CURVES, asd)
                ),
                minimum_frequency=minimum_frequency,
                maximum_frequency=maximum_frequency,
                length=length,
                latitude=latitude,
                longitude=longitude,
                elevation=0.0,
                xarm_azimuth=xarm,
                yarm_azimuth=xarm + 90.0,
            )
        )
    return bilby.gw.detector.InterferometerList(interferometers)


def gwfast_network(minimum_frequency, maximum_frequency, approximant="IMRPhenomXHM"):
    """Build the same network as a GWFast :class:`DetNet`.

    Parameters
    ----------
    minimum_frequency, maximum_frequency : float
        Band limits in Hz.
    approximant : str
        Any frequency-domain LAL approximant.

    Returns
    -------
    gwfast.network.DetNet

    Notes
    -----
    ``compute_sequence=False`` on purpose. GWFast's own docstring warns that the
    faster ``SimInspiralChooseFDWaveformSequence`` path "shows numerical issues
    with some waveform models (e.g. IMRPhenomXHM)", which is precisely the model
    under test here.
    """
    patch_gwfast()
    from gwfast.network import DetNet
    from gwfast.signal import GWSignal
    from gwfast.waveforms import LAL_WF

    # is_HigherModes only matters for models that have them; setting it for a
    # quadrupole-only approximant is harmless and keeps the call uniform.
    waveform = LAL_WF(approximant, is_HigherModes=True, compute_sequence=False)
    signals = {}
    for name, latitude, longitude, xarm, _, asd in DETECTORS:
        signals[name] = GWSignal(
            waveform,
            psd_path=os.path.join(NOISE_CURVES, asd),
            detector_shape="L",
            det_lat=latitude,
            det_long=longitude,
            # GWFast orients by the bisector of the arms, bilby by the x arm.
            det_xax=xarm + 45.0,
            is_ASD=True,
            useEarthMotion=False,
            fmin=minimum_frequency,
            fmax=maximum_frequency,
            verbose=False,
        )
    return DetNet(signals, verbose=False), waveform


def to_gwfast_parameters(parameters):
    """Translate one GWForge injection into GWFast's convention.

    Parameters
    ----------
    parameters : dict
        GWForge parameters, including ``geocent_time`` as a GPS time.

    Returns
    -------
    dict
        Arrays of length one, as GWFast expects.
    """
    patch_gwfast()
    from gwfast.gwfastUtils import GPSt_to_LMST

    values = {
        "Mc": parameters["chirp_mass"],
        "eta": parameters["symmetric_mass_ratio"],
        "dL": parameters["luminosity_distance"] * 1e-3,
        "theta": numpy.pi / 2.0 - parameters["dec"],
        "phi": parameters["ra"],
        "iota": parameters["theta_jn"],
        "psi": parameters["psi"],
        # GWFast measures the coalescence time as a sidereal fraction of a day.
        "tcoal": GPSt_to_LMST(parameters["geocent_time"], 0.0, 0.0),
        "Phicoal": parameters["phase"],
        "chi1z": parameters["chi_1"],
        "chi2z": parameters["chi_2"],
    }
    return {key: numpy.array([value]) for key, value in values.items()}


def to_gwforge_basis(fisher, gwfast_parameter_numbers):
    """Rotate a GWFast Fisher matrix into GWForge's parameters.

    A Fisher matrix transforms as :math:`F' = J^T F J` with
    :math:`J_{ij} = \\partial p^{\\rm gwfast}_i / \\partial p^{\\rm gwforge}_j`.
    Here every GWFast parameter is a function of exactly one GWForge parameter,
    so ``J`` is a permutation times a diagonal.

    Parameters
    ----------
    fisher : numpy.ndarray
        ``(n, n)`` GWFast Fisher matrix.
    gwfast_parameter_numbers : dict
        GWFast's ``ParNums``, mapping its parameter names to matrix indices.

    Returns
    -------
    numpy.ndarray
        ``(n, n)`` matrix in the order of :data:`PARAMETERS`.
    """
    count = len(PARAMETER_MAP)
    jacobian = numpy.zeros((count, count))
    for column, (_, gwfast_name, factor) in enumerate(PARAMETER_MAP):
        jacobian[gwfast_parameter_numbers[gwfast_name], column] = factor
    return jacobian.T @ numpy.asarray(fisher) @ jacobian
