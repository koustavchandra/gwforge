"""The Fisher module against waveforms it was not developed on.

:mod:`GWForge.fisher` was developed against ``IMRPhenomXHM`` — an aligned-spin,
quasi-circular LAL waveform reached through ``lal_binary_black_hole``. Its whole
point, though, is that it works with any
waveform :class:`bilby.gw.WaveformGenerator` can produce, and a single
well-behaved model is a weak test of that. These four each break something
different:

===================  ==============================  ===============================
waveform             source model                    what it breaks
===================  ==============================  ===============================
``IMRPhenomXHM``     ``lal_binary_black_hole``       nothing — the control
``SEOBNRv5HM``       ``gwsignal_binary_black_hole``  a second, non-LAL backend
``SEOBNRv5PHM``      ``gwsignal_binary_black_hole``  precession
``pyEFPEHM``         ``EFPE_binary_black_hole``      precession + eccentricity,
                                                     inspiral only
===================  ==============================  ===============================

The control matters: when a check fails for one of the other three, its passing
here is what says the harness is fine and the waveform is the interesting part.

Everything is order 1. ``SEOBNRv5HM``/``SEOBNRv5PHM`` need ``pyseobnr``
(``pip install -e .[waveforms]``) and ``pyEFPEHM`` is installed from git; both
are skipped rather than failed when absent.
"""

import warnings

import numpy
import pytest

import bilby
from GWForge.fisher import FisherMatrix, analytic_derivatives, covariance
from GWForge.fisher.derivatives import WaveformDerivatives
from GWForge.fisher.parameters import (
    DEFAULT_PARAMETERS,
    ECCENTRIC_PRECESSING_PARAMETERS,
    PRECESSING_PARAMETERS,
)
from GWForge.ifo import antenna
from GWForge.ifo.detectors import Network

bilby.core.utils.setup_logger(log_level="error")

# bilby 2.8.1 deprecates gwsignal_binary_black_hole and logs it on every call.
# A Fisher makes ~50 calls per source, which would bury everything else.
warnings.filterwarnings("ignore", message=".*gwsignal_binary_black_hole is deprecated.*")
warnings.filterwarnings("ignore", message=".*UNREVIEWED.*")

GPS_TIME = 1187008882.4

# Extrinsic parameters shared by every configuration.
EXTRINSIC = dict(
    luminosity_distance=1000.0,
    theta_jn=0.9,
    psi=0.9,
    phase=1.3,
    geocent_time=GPS_TIME,
    ra=1.375,
    dec=-1.2108,
)

# Spins chosen well away from zero. A tilt of exactly zero makes ``phi_12`` and
# ``phi_jl`` unmeasurable and the Fisher singular, so a precessing test at zero
# tilt would be testing nothing.
PRECESSING_SPINS = dict(
    a_1=0.5, tilt_1=0.6, phi_12=1.7, a_2=0.4, tilt_2=1.0, phi_jl=0.3
)
ALIGNED_SPINS = dict(chi_1=0.3, chi_2=-0.2)


def _configuration(
    label,
    approximant,
    source_model,
    parameters,
    fisher_parameters,
    detectors=("CE40", "CE20"),
    duration=8.0,
    sampling_frequency=2048.0,
    minimum_frequency=20.0,
    waveform_arguments=None,
    requires=None,
    sigma_tolerance=1.0,
):
    return dict(
        label=label,
        approximant=approximant,
        source_model=source_model,
        parameters={**EXTRINSIC, **parameters},
        fisher_parameters=list(fisher_parameters),
        detectors=list(detectors),
        duration=duration,
        sampling_frequency=sampling_frequency,
        minimum_frequency=minimum_frequency,
        waveform_arguments=waveform_arguments,
        requires=requires,
        sigma_tolerance=sigma_tolerance,
    )


CONFIGURATIONS = [
    _configuration(
        "IMRPhenomXHM",
        "IMRPhenomXHM",
        "lal_binary_black_hole",
        dict(chirp_mass=28.0, symmetric_mass_ratio=0.24, **ALIGNED_SPINS),
        DEFAULT_PARAMETERS,
        # Measured: sigma(chi_1) is flat to 0.03% across steps from 3e-4 to
        # 1e-2, so a closed-form waveform is held to a tight bound.
        sigma_tolerance=0.10,
    ),
    _configuration(
        "SEOBNRv5HM",
        "SEOBNRv5HM",
        "gwsignal_binary_black_hole",
        dict(chirp_mass=28.0, symmetric_mass_ratio=0.24, **ALIGNED_SPINS),
        DEFAULT_PARAMETERS,
        requires="pyseobnr",
    ),
    _configuration(
        "SEOBNRv5PHM",
        "SEOBNRv5PHM",
        "gwsignal_binary_black_hole",
        dict(chirp_mass=28.0, symmetric_mass_ratio=0.24, **PRECESSING_SPINS),
        PRECESSING_PARAMETERS,
        requires="pyseobnr",
    ),
    _configuration(
        "pyEFPEHM",
        "pyEFPEHM",
        "EFPE_binary_black_hole",
        dict(
            chirp_mass=28.0,
            symmetric_mass_ratio=0.24,
            eccentricity=0.15,
            mean_anomaly=0.7,
            luminosity_distance=400.0,
            **PRECESSING_SPINS,
        ),
        ECCENTRIC_PRECESSING_PARAMETERS,
        # An inspiral-only post-Newtonian model belongs on current detectors.
        detectors=("H1", "L1", "V1"),
        # Generate from 10 Hz but integrate from 20 Hz. pyEFPEHM's eccentric
        # harmonics radiate *below* the quadrupole start frequency (measured
        # support 4.3-136.5 Hz for f22_start = 10), so starting the waveform at
        # 20 Hz would switch harmonics on inside the analysis band.
        minimum_frequency=20.0,
        waveform_arguments=dict(f22_start=10.0, reference_frequency=10.0),
        # The waveform is built directly in the frequency domain, so duration
        # buys resolution rather than avoiding wrap-around; an eccentric
        # spectrum has structure that a coarse grid smears.
        duration=64.0,
        requires="pyEFPEHM",
    ),
]

IDS = [configuration["label"] for configuration in CONFIGURATIONS]


@pytest.fixture(scope="module", params=CONFIGURATIONS, ids=IDS)
def setup(request):
    """Network, waveform generator and parameters for one configuration."""
    configuration = request.param
    if configuration["requires"]:
        pytest.importorskip(
            configuration["requires"],
            reason="{} needs {}".format(
                configuration["label"], configuration["requires"]
            ),
        )

    from GWForge.fisher.config import SOURCE_MODELS

    source_model = SOURCE_MODELS[configuration["source_model"]]
    if source_model is None:
        pytest.skip("{} is unavailable".format(configuration["source_model"]))

    interferometers = Network(ifos=configuration["detectors"]).initialise_ifos()
    for interferometer in interferometers:
        interferometer.minimum_frequency = configuration["minimum_frequency"]
    interferometers.set_strain_data_from_zero_noise(
        sampling_frequency=configuration["sampling_frequency"],
        duration=configuration["duration"],
        start_time=GPS_TIME - configuration["duration"] + 2,
    )

    arguments = dict(
        waveform_approximant=configuration["approximant"],
        reference_frequency=configuration["minimum_frequency"],
        minimum_frequency=configuration["minimum_frequency"],
    )
    if configuration["waveform_arguments"]:
        arguments.update(configuration["waveform_arguments"])

    waveform_generator = bilby.gw.WaveformGenerator(
        duration=configuration["duration"],
        sampling_frequency=configuration["sampling_frequency"],
        frequency_domain_source_model=source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=arguments,
    )
    return dict(
        interferometers=interferometers,
        waveform_generator=waveform_generator,
        parameters=configuration["parameters"],
        fisher_parameters=configuration["fisher_parameters"],
        label=configuration["label"],
        sigma_tolerance=configuration["sigma_tolerance"],
    )


@pytest.fixture(scope="module")
def fisher(setup):
    """One order-1 Fisher matrix per configuration, built once."""
    return FisherMatrix(
        setup["interferometers"],
        setup["waveform_generator"],
        setup["parameters"],
        fisher_parameters=setup["fisher_parameters"],
        order=1,
    )


def relative_error(actual, expected, mask):
    scale = numpy.max(numpy.abs(expected[mask]))
    return numpy.max(numpy.abs(actual[mask] - expected[mask])) / scale


# ---------------------------------------------------------------------------
# The waveform itself
# ---------------------------------------------------------------------------


def test_waveform_generates_with_support_in_band(setup):
    """A waveform that is silently zero would make every later check vacuous."""
    interferometer = setup["interferometers"][0]
    polarizations = setup["waveform_generator"].frequency_domain_strain(
        dict(setup["parameters"])
    )
    mask = interferometer.frequency_mask
    for polarization in ("plus", "cross"):
        strain = polarizations[polarization]
        assert numpy.all(numpy.isfinite(strain))
        assert numpy.any(numpy.abs(strain[mask]) > 0), polarization


def test_waveform_errors_are_not_swallowed(setup):
    """``catch_waveform_errors`` must stay off.

    Both ``gwsignal_binary_black_hole`` and ``EFPE_binary_black_hole`` return
    ``None`` for a caught failure. In a likelihood that becomes ``-inf``; in a
    finite-difference stencil it would quietly corrupt a derivative, so the
    Fisher must see the exception instead.
    """
    assert not setup["waveform_generator"].waveform_arguments.get(
        "catch_waveform_errors", False
    )


# ---------------------------------------------------------------------------
# The claim that the analytic derivatives rest on
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("parameter", ["ra", "dec", "psi", "geocent_time"])
def test_polarizations_do_not_depend_on_the_geometric_parameters(setup, parameter):
    """Differencing only the response reproduces the full derivative exactly.

    This is the load-bearing assumption behind the closed-form derivatives, and
    it is a statement about the detector projection rather than about any
    waveform — so it has to hold for *every* source model, whatever the physics
    inside. Bit-identical, not merely close.
    """
    interferometer = setup["interferometers"][0]
    parameters = setup["parameters"]
    fiducial = setup["waveform_generator"].frequency_domain_strain(dict(parameters))
    step = 1e-4 if parameter != "geocent_time" else 1e-3

    def full(shift):
        shifted = dict(parameters)
        shifted[parameter] += shift
        return antenna.detector_response(
            interferometer,
            setup["waveform_generator"].frequency_domain_strain(dict(shifted)),
            shifted,
        )

    def geometry_only(shift):
        shifted = dict(parameters)
        shifted[parameter] += shift
        return antenna.detector_response(interferometer, fiducial, shifted)

    numpy.testing.assert_array_equal(
        (geometry_only(step) - geometry_only(-step)) / (2 * step),
        (full(step) - full(-step)) / (2 * step),
    )


@pytest.mark.parametrize("parameter", ["ra", "dec", "psi"])
def test_analytic_first_derivatives(setup, parameter):
    """Closed-form sky and polarisation derivatives against finite differences."""
    interferometer = setup["interferometers"][0]
    parameters = setup["parameters"]
    polarizations = setup["waveform_generator"].frequency_domain_strain(
        dict(parameters)
    )
    _, first, _ = analytic_derivatives(
        interferometer, polarizations, parameters, order=1
    )

    def central(width):
        def response(shift):
            shifted = dict(parameters)
            shifted[parameter] += shift
            return antenna.detector_response(interferometer, polarizations, shifted)

        return (response(width) - response(-width)) / (2 * width)

    reference = (4 * central(5e-5) - central(1e-4)) / 3
    assert (
        relative_error(first[parameter], reference, interferometer.frequency_mask)
        < 1e-7
    )


def test_exact_distance_and_polarisation_identities(setup):
    """``dh/dd_L = -h/d_L`` and ``d2h/dpsi2 = -4h``, whatever the waveform."""
    interferometer = setup["interferometers"][0]
    parameters = setup["parameters"]
    polarizations = setup["waveform_generator"].frequency_domain_strain(
        dict(parameters)
    )
    response, first, second = analytic_derivatives(
        interferometer, polarizations, parameters, order=2
    )
    distance = parameters["luminosity_distance"]
    numpy.testing.assert_allclose(
        first["luminosity_distance"], -response / distance, rtol=1e-14
    )
    numpy.testing.assert_allclose(
        second[("psi", "psi")], -4.0 * response, rtol=1e-10, atol=1e-40
    )


# ---------------------------------------------------------------------------
# The Fisher matrix
# ---------------------------------------------------------------------------


def test_fisher_is_symmetric_and_positive_definite(fisher, setup):
    numpy.testing.assert_allclose(fisher.matrix, fisher.matrix.T, rtol=1e-12)
    smallest = numpy.linalg.eigvalsh(fisher.matrix).min()
    assert smallest > 0, "{}: smallest eigenvalue {:.3e}".format(
        setup["label"], smallest
    )


def test_covariance_round_trip(fisher, setup):
    matrix, diagnostics = covariance(fisher.matrix, names=fisher.names)
    numpy.testing.assert_allclose(
        matrix @ fisher.matrix, numpy.eye(len(fisher.names)), atol=1e-5
    )
    assert numpy.all(numpy.isfinite(list(fisher.sigmas.values())))
    assert all(sigma > 0 for sigma in fisher.sigmas.values())
    assert diagnostics["condition_number"] < 1e15, setup["label"]


def test_snr_agrees_with_bilby(fisher, setup):
    """The Fisher SNR must equal what bilby reports for the same signal."""
    polarizations = setup["waveform_generator"].frequency_domain_strain(
        dict(setup["parameters"])
    )
    for interferometer in setup["interferometers"]:
        signal = antenna.detector_response(
            interferometer, polarizations, setup["parameters"]
        )
        expected = numpy.sqrt(interferometer.optimal_snr_squared(signal=signal).real)
        assert fisher.optimal_snrs[interferometer.name] == pytest.approx(
            expected, rel=1e-10
        )


def test_every_parameter_is_measured(fisher, setup):
    """No parameter may leave the waveform untouched.

    A zero on the Fisher diagonal means the derivative vanished — a parameter
    the source model silently ignored, which is exactly how a name that a model
    does not understand fails.
    """
    diagonal = numpy.diag(fisher.matrix)
    dead = [
        name for name, value in zip(fisher.names, diagonal) if value <= 0 or not numpy.isfinite(value)
    ]
    assert not dead, "{}: no response to {}".format(setup["label"], dead)


# Noise-to-signal at which a derivative stops meaning anything. Measured
# values are far below it -- 0.3% at worst, for SEOBNRv5HM's masses and spins
# -- so in practice nothing is excluded and the reproducibility check below
# applies to every parameter. It is kept as a guard for a future model that is
# genuinely undifferentiable.
NOISE_DOMINATED = 0.05

# Below this a derivative is limited by something other than the model's noise.
NEGLIGIBLE_NOISE = 1e-4


def test_noise_is_measured_and_reported(fisher, setup):
    """Whether a model can be differentiated is measured, not assumed.

    A finite difference cannot see beneath the model's own numerical noise.
    ``WaveformDerivatives.noise_ratios`` records that noise divided by the
    change the calibrated step produces -- roughly the relative error of the
    derivative. Measured here:

    ================  ==========================================================
    model             worst relative error
    ================  ==========================================================
    IMRPhenomXHM      4e-6  (essentially machine precision)
    SEOBNRv5HM        3e-3  on masses and spins; 5e-14 on theta_jn and phase
    SEOBNRv5PHM       2e-4  on the precession angles
    pyEFPEHM          2e-4  on chirp_mass
    ================  ==========================================================

    The SEOBNRv5HM pattern is the one that says the estimator is measuring
    something real rather than producing numbers: masses and spins go through
    the adaptive EOB ODE solve and pick up scatter, while ``theta_jn`` and
    ``phase`` are applied analytically outside it and come out clean by ten
    orders of magnitude.
    """
    ratios = fisher.derivatives.noise_ratios
    assert set(ratios) == set(fisher.derivatives.numerical)
    assert all(value >= 0 for value in ratios.values())


@pytest.fixture(scope="module")
def repeated_sigmas(setup):
    """Sigmas from three different starting steps, computed once per model."""
    results = {}
    for factor in (1.0, 0.1, 10.0):
        guesses = (
            None
            if factor == 1.0
            else {
                name: factor * 1e-5
                for name in setup["fisher_parameters"]
                if name
                not in ("luminosity_distance", "ra", "dec", "psi", "geocent_time")
            }
        )
        matrix = FisherMatrix(
            setup["interferometers"],
            setup["waveform_generator"],
            setup["parameters"],
            fisher_parameters=setup["fisher_parameters"],
            order=1,
            step_sizes=guesses,
        )
        results[factor] = (matrix.sigmas, dict(matrix.derivatives.noise_ratios))
    return results


def test_calibration_is_reproducible_where_the_waveform_allows_it(
    repeated_sigmas, setup
):
    """Different starting guesses must land on the same Fisher.

    This is the check that caught a 150% error on ``sigma(eta)`` for
    IMRPhenomXHM before the step calibration existed. How tight it can be is a
    property of the waveform, so the bound is declared per model and measured
    rather than assumed.

    A closed-form waveform supports a tight bound: scanning IMRPhenomXHM's
    ``chi_1`` step from 3e-4 to 1e-2, with every other step calibrated, moves
    ``sigma(chi_1)`` by 0.03%. The same scan for SEOBNRv5HM moves it by 70%,
    with only a narrow flat region around 3e-4 to 1e-3. That is **not** the
    model's numerical noise -- :func:`estimate_noise` puts that at 0.3%, ten
    times smaller. It is the derivative genuinely varying with the step, i.e.
    truncation error from structure the fourth-order stencil does not capture,
    so no amount of noise-aware step selection removes it.

    The bound therefore stays loose for those models, and its job is to catch an
    outright breakdown rather than to pretend to a precision the waveform cannot
    support.
    """
    reference, ratios = repeated_sigmas[1.0]
    tolerance = setup["sigma_tolerance"]
    worst = {
        name: max(run_ratios.get(name, 0.0) for _, run_ratios in repeated_sigmas.values())
        for name in ratios
    }
    for factor in (0.1, 10.0):
        alternative, _ = repeated_sigmas[factor]
        for name in reference:
            if worst.get(name, 0.0) >= NOISE_DOMINATED:
                continue
            assert alternative[name] == pytest.approx(
                reference[name], rel=tolerance
            ), "{}: sigma({}) moved by more than {:.0%} when the starting step was scaled by {}".format(
                setup["label"], name, tolerance, factor
            )
    assert any(
        worst.get(name, 0.0) < NOISE_DOMINATED for name in reference
    ), "{}: every derivative is noise-dominated, so this test asserted nothing".format(
        setup["label"]
    )


def test_the_control_waveform_is_smooth(repeated_sigmas, setup):
    """IMRPhenomXHM must come out at machine precision.

    A closed-form waveform has no adaptive internals to scatter, so if this ever
    reported real noise the estimator would have started seeing things, and the
    exclusions above would quietly stop testing anything.
    """
    if setup["label"] != "IMRPhenomXHM":
        pytest.skip("only the control is a closed-form waveform")
    for _, ratios in repeated_sigmas.values():
        assert all(value < NEGLIGIBLE_NOISE for value in ratios.values()), ratios


def test_analytic_parameters_cost_no_waveform_evaluations(setup):
    """The five closed-form parameters must stay free on every source model."""
    numerical = [
        name
        for name in setup["fisher_parameters"]
        if name in ("chirp_mass", "symmetric_mass_ratio")
    ]
    without = WaveformDerivatives(
        setup["interferometers"],
        setup["waveform_generator"],
        setup["parameters"],
        fisher_parameters=numerical,
        order=1,
    )
    with_analytic = WaveformDerivatives(
        setup["interferometers"],
        setup["waveform_generator"],
        setup["parameters"],
        fisher_parameters=numerical
        + ["luminosity_distance", "ra", "dec", "psi", "geocent_time"],
        order=1,
    )
    assert with_analytic.waveform_evaluations == without.waveform_evaluations
