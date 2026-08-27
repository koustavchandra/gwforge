"""Tests for the Fisher and DALI forecasting in :mod:`GWForge.fisher`.

The theme running through most of these is that five parameters --
``luminosity_distance``, ``geocent_time``, ``psi``, ``ra``, ``dec`` -- have
*exact* derivatives, because a bilby source model returns polarizations that do
not depend on them at all. That claim is what makes the module fast and, for
``geocent_time``, what makes it correct at all, so it is checked directly rather
than assumed.
"""

import numpy
import pytest

import bilby
from GWForge.fisher import (
    ANALYTIC_PARAMETERS,
    FisherMatrix,
    analytic_derivatives,
    covariance,
    inner_product,
)
from GWForge.fisher.dali import DALILikelihood, build_priors
from GWForge.fisher.antenna_derivatives import AntennaDerivatives
from GWForge.fisher.derivatives import WaveformDerivatives
from GWForge.fisher.parameters import bounded_step, strip_shadowed_parameters
from GWForge.fisher.matrix import check_order
from GWForge.ifo import antenna
from GWForge.ifo.detectors import Network

bilby.core.utils.setup_logger(log_level="error")

GPS_TIME = 1126259642.413
DURATION = 8.0
SAMPLING_FREQUENCY = 1024.0
MINIMUM_FREQUENCY = 20.0

PARAMETERS = dict(
    chirp_mass=28.0,
    symmetric_mass_ratio=0.24,
    chi_1=0.3,
    chi_2=-0.2,
    luminosity_distance=1500.0,
    theta_jn=0.7,
    psi=0.9,
    phase=1.3,
    geocent_time=GPS_TIME,
    ra=1.375,
    dec=-1.2108,
)

RESPONSE_OPTIONS = [
    pytest.param(False, False, id="static"),
    pytest.param(True, False, id="earth-rotation"),
    pytest.param(False, True, id="finite-size"),
    pytest.param(True, True, id="both"),
]


@pytest.fixture(scope="module")
def network():
    ifos = Network(ifos=["CE40", "CE20"]).initialise_ifos()
    ifos.set_strain_data_from_zero_noise(
        sampling_frequency=SAMPLING_FREQUENCY,
        duration=DURATION,
        start_time=GPS_TIME - DURATION + 2,
    )
    return ifos


@pytest.fixture(scope="module")
def waveform_generator():
    return bilby.gw.WaveformGenerator(
        duration=DURATION,
        sampling_frequency=SAMPLING_FREQUENCY,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=dict(
            waveform_approximant="IMRPhenomXHM",
            reference_frequency=MINIMUM_FREQUENCY,
            minimum_frequency=MINIMUM_FREQUENCY,
        ),
    )


def relative_error(actual, expected, mask):
    scale = numpy.max(numpy.abs(expected[mask]))
    return numpy.max(numpy.abs(actual[mask] - expected[mask])) / scale


# ---------------------------------------------------------------------------
# The claim that makes the analytic derivatives possible
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("earth_rotation, finite_size", RESPONSE_OPTIONS)
@pytest.mark.parametrize("parameter", ["ra", "dec", "psi", "geocent_time"])
def test_polarizations_do_not_depend_on_the_geometric_parameters(
    network, waveform_generator, parameter, earth_rotation, finite_size
):
    """Differencing only the response reproduces the full derivative exactly.

    This is the load-bearing assumption of the module: ``ra``, ``dec``, ``psi``
    and ``geocent_time`` reach the strain only through the detector projection.
    Holding the polarizations fixed and varying only the response must therefore
    give *bit-identical* results to varying everything, and it does -- including
    under the frequency-dependent, Earth-rotating response, where
    ``F_plus``/``F_cross`` are complex and vary across the band.
    """
    interferometer = network[0]
    options = dict(earth_rotation=earth_rotation, finite_size=finite_size)
    fiducial = waveform_generator.frequency_domain_strain(dict(PARAMETERS))
    step = 1e-4 if parameter != "geocent_time" else 1e-3

    def full(shift):
        shifted = dict(PARAMETERS)
        shifted[parameter] += shift
        return antenna.detector_response(
            interferometer,
            waveform_generator.frequency_domain_strain(dict(shifted)),
            shifted,
            **options,
        )

    def geometry_only(shift):
        shifted = dict(PARAMETERS)
        shifted[parameter] += shift
        return antenna.detector_response(interferometer, fiducial, shifted, **options)

    from_everything = (full(step) - full(-step)) / (2 * step)
    from_geometry = (geometry_only(step) - geometry_only(-step)) / (2 * step)
    numpy.testing.assert_array_equal(from_geometry, from_everything)


@pytest.mark.parametrize("earth_rotation, finite_size", RESPONSE_OPTIONS)
def test_antenna_derivative_values_match_the_response(
    network, earth_rotation, finite_size
):
    """The derivative class must rebuild the response it differentiates."""
    interferometer = network[0]
    frequencies = numpy.array([21.0, 55.0, 213.0, 400.0])
    times_to_coalescence = antenna.time_to_coalescence(frequencies, 28.0, 0.24)
    options = dict(earth_rotation=earth_rotation, finite_size=finite_size)

    expected = antenna.antenna_response(
        interferometer,
        PARAMETERS["ra"],
        PARAMETERS["dec"],
        PARAMETERS["geocent_time"],
        PARAMETERS["psi"],
        frequencies,
        interferometer.strain_data.start_time,
        times_to_coalescence=times_to_coalescence,
        **options,
    )
    derivatives = AntennaDerivatives(
        interferometer,
        ra=PARAMETERS["ra"],
        dec=PARAMETERS["dec"],
        geocent_time=PARAMETERS["geocent_time"],
        psi=PARAMETERS["psi"],
        frequencies=frequencies,
        start_time=interferometer.strain_data.start_time,
        times_to_coalescence=times_to_coalescence,
        **options,
    )
    for actual, reference in zip(
        (derivatives.plus["value"], derivatives.cross["value"], derivatives.delay["value"]),
        expected,
    ):
        numpy.testing.assert_allclose(actual, reference, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("earth_rotation, finite_size", RESPONSE_OPTIONS)
@pytest.mark.parametrize("parameter", ["ra", "dec", "psi"])
def test_analytic_first_derivatives(
    network, waveform_generator, parameter, earth_rotation, finite_size
):
    """Closed-form sky and polarisation derivatives against finite differences.

    ``geocent_time`` is excluded on purpose: no finite difference can check it,
    which is precisely why it is done analytically. See
    :func:`test_geocent_time_derivative_beats_finite_differencing`.
    """
    interferometer = network[0]
    options = dict(earth_rotation=earth_rotation, finite_size=finite_size)
    polarizations = waveform_generator.frequency_domain_strain(dict(PARAMETERS))
    _, first, _ = analytic_derivatives(
        interferometer, polarizations, PARAMETERS, order=1, **options
    )

    def response(shift):
        shifted = dict(PARAMETERS)
        shifted[parameter] += shift
        return antenna.detector_response(
            interferometer, polarizations, shifted, **options
        )

    step = 1e-4

    def central(width):
        return (response(width) - response(-width)) / (2 * width)

    # Richardson extrapolation, so the reference is fourth-order accurate.
    reference = (4 * central(step / 2) - central(step)) / 3
    assert (
        relative_error(first[parameter], reference, interferometer.frequency_mask)
        < 1e-7
    )


def test_luminosity_distance_derivatives_are_exact(network, waveform_generator):
    """The waveform scales as 1/d_L, so both derivatives are rescalings."""
    interferometer = network[0]
    polarizations = waveform_generator.frequency_domain_strain(dict(PARAMETERS))
    response, first, second = analytic_derivatives(
        interferometer, polarizations, PARAMETERS, order=2
    )
    distance = PARAMETERS["luminosity_distance"]
    numpy.testing.assert_allclose(
        first["luminosity_distance"], -response / distance, rtol=1e-14
    )
    numpy.testing.assert_allclose(
        second[("luminosity_distance", "luminosity_distance")],
        2 * response / distance**2,
        rtol=1e-14,
    )


def test_second_derivative_in_psi_is_minus_four_times_the_response(
    network, waveform_generator
):
    """(F_plus, F_cross) is a rigid rotation by 2*psi, so d2h/dpsi2 = -4h.

    This holds bin by bin even when the transfer functions make the pattern
    functions complex and frequency-dependent, because they multiply the
    contracted tensors linearly.
    """
    interferometer = network[0]
    polarizations = waveform_generator.frequency_domain_strain(dict(PARAMETERS))
    response, _, second = analytic_derivatives(
        interferometer, polarizations, PARAMETERS, order=2
    )
    numpy.testing.assert_allclose(
        second[("psi", "psi")], -4.0 * response, rtol=1e-10, atol=1e-40
    )


def test_geocent_time_derivative_beats_finite_differencing(network, waveform_generator):
    """No step size lets a finite difference reach the analytic t_c derivative.

    ``geocent_time`` is a GPS number of order 1e9, so the segment-relative delay
    is a catastrophic subtraction, while the derivative carries an oscillatory
    ``2 pi f``. Roundoff falls as the step grows and truncation rises, and the
    two never leave a usable window -- the best any step manages is around a
    part in a thousand. This test asserts that the floor exists, which is the
    justification for the closed form.
    """
    interferometer = network[0]
    mask = interferometer.frequency_mask
    polarizations = waveform_generator.frequency_domain_strain(dict(PARAMETERS))
    _, first, _ = analytic_derivatives(
        interferometer, polarizations, PARAMETERS, order=1
    )

    def response(shift):
        shifted = dict(PARAMETERS)
        shifted["geocent_time"] += shift
        return antenna.detector_response(interferometer, polarizations, shifted)

    best = min(
        relative_error(
            (response(step) - response(-step)) / (2 * step), first["geocent_time"], mask
        )
        for step in (1e-6, 1e-5, 1e-4, 1e-3)
    )
    assert best > 1e-5, (
        "a finite difference matched the analytic geocent_time derivative better "
        "than expected; if the response no longer subtracts the segment start "
        "time this test's premise has changed"
    )
    assert best < 1e-1


# ---------------------------------------------------------------------------
# Fisher matrix
# ---------------------------------------------------------------------------


def test_fisher_of_a_linear_signal_matches_its_analytic_covariance(network):
    """For a signal linear in its parameters the Fisher matrix is exact.

    Building ``h = a * u + b * v`` from two fixed frequency series makes
    ``dh/da = u`` and ``dh/db = v``, so the Fisher matrix must be the Gram
    matrix of ``u`` and ``v`` under the noise-weighted inner product. That
    checks the inner product, the accumulation and the inversion together,
    without any waveform involved.
    """
    interferometer = network[0]
    mask = interferometer.frequency_mask
    frequencies = interferometer.frequency_array

    basis = []
    for power in (-7.0 / 6.0, -1.0 / 6.0):
        series = numpy.zeros(len(frequencies), dtype=complex)
        series[mask] = 1e-23 * frequencies[mask] ** power * numpy.exp(
            1j * frequencies[mask] / 40.0
        )
        basis.append(series)

    expected = numpy.array(
        [[inner_product(left, right, interferometer) for right in basis] for left in basis]
    )
    inverse, diagnostics = covariance(expected)
    numpy.testing.assert_allclose(inverse @ expected, numpy.eye(2), atol=1e-10)
    assert diagnostics["inversion_error"] < 1e-8


def test_covariance_round_trip(network, waveform_generator):
    """Sigma times F is the identity, and sigmas are the diagonal square roots."""
    fisher = FisherMatrix(network, waveform_generator, PARAMETERS, order=1)
    matrix, diagnostics = covariance(fisher.matrix, names=fisher.names)
    numpy.testing.assert_allclose(
        matrix @ fisher.matrix, numpy.eye(len(fisher.names)), atol=1e-6
    )
    assert diagnostics["condition_number"] > 1
    for index, name in enumerate(fisher.names):
        assert fisher.sigmas[name] == pytest.approx(numpy.sqrt(matrix[index, index]))


def test_fisher_matrix_is_symmetric_and_positive_definite(network, waveform_generator):
    fisher = FisherMatrix(network, waveform_generator, PARAMETERS, order=1)
    numpy.testing.assert_allclose(fisher.matrix, fisher.matrix.T, rtol=1e-12)
    assert numpy.linalg.eigvalsh(fisher.matrix).min() > 0


def test_network_fisher_is_the_sum_over_detectors(network, waveform_generator):
    fisher = FisherMatrix(network, waveform_generator, PARAMETERS, order=1)
    total = sum(fisher.per_detector.values())
    numpy.testing.assert_allclose(fisher.matrix, total, rtol=1e-12)


def test_snr_agrees_with_bilby(network, waveform_generator):
    """The Fisher SNR must equal the one bilby reports for the same injection."""
    fisher = FisherMatrix(network, waveform_generator, PARAMETERS, order=1)
    for interferometer in network:
        signal = antenna.detector_response(
            interferometer,
            waveform_generator.frequency_domain_strain(dict(PARAMETERS)),
            PARAMETERS,
        )
        expected = numpy.sqrt(interferometer.optimal_snr_squared(signal=signal).real)
        assert fisher.optimal_snrs[interferometer.name] == pytest.approx(
            expected, rel=1e-10
        )


def test_gaussian_prior_only_touches_the_diagonal(network, waveform_generator):
    fisher = FisherMatrix(network, waveform_generator, PARAMETERS, order=1)
    combined = fisher.with_gaussian_prior({"luminosity_distance": 100.0})
    index = fisher.names.index("luminosity_distance")
    difference = combined - fisher.matrix
    assert difference[index, index] == pytest.approx(1.0 / 100.0**2)
    difference[index, index] = 0.0
    numpy.testing.assert_allclose(difference, 0.0, atol=0)
    # A prior can only shrink the uncertainty.
    assert covariance(combined)[0][index, index] < fisher.covariance[index, index]


def test_singular_fisher_matrix_is_reported_clearly():
    with pytest.raises(numpy.linalg.LinAlgError, match="non-positive diagonal"):
        covariance(numpy.array([[1.0, 0.0], [0.0, 0.0]]), names=["good", "unconstrained"])


# ---------------------------------------------------------------------------
# Step sizes
# ---------------------------------------------------------------------------


def test_steps_stay_inside_physical_bounds():
    """A near equal-mass binary must not be stepped past eta = 0.25."""
    step = bounded_step("symmetric_mass_ratio", 0.2499, 1e-2)
    assert 0.0 < step <= 0.4 * (0.25 - 0.2499)
    assert 0.2499 + 2 * step < 0.25
    # Unbounded parameters are untouched.
    assert bounded_step("phase", 1.3, 1e-2) == 1e-2


def test_calibration_makes_the_fisher_insensitive_to_the_starting_step(
    network, waveform_generator
):
    """Different starting guesses must converge on the same Fisher matrix.

    Without calibration this is badly false: the sigma on
    ``symmetric_mass_ratio`` moves by more than a factor of two across a decade
    of step size, because LAL places each higher-mode's onset at a mass-dependent
    frequency and differencing across that edge contributes a spurious term
    growing like ``1 / step``.
    """
    reference = FisherMatrix(network, waveform_generator, PARAMETERS, order=1).sigmas
    for factor in (0.1, 10.0):
        guesses = {
            name: factor * value
            for name, value in dict(
                chirp_mass=1e-6, symmetric_mass_ratio=1e-6, chi_1=1e-5, chi_2=1e-5
            ).items()
        }
        alternative = FisherMatrix(
            network, waveform_generator, PARAMETERS, order=1, step_sizes=guesses
        ).sigmas
        for name in reference:
            assert alternative[name] == pytest.approx(reference[name], rel=0.05)


# ---------------------------------------------------------------------------
# DALI
# ---------------------------------------------------------------------------


def test_order_three_is_refused_with_an_explanation():
    with pytest.raises(NotImplementedError, match="third derivatives"):
        check_order(3)
    for order in (1, 2):
        assert check_order(order) == order
    with pytest.raises(ValueError):
        check_order(4)


@pytest.fixture(scope="module")
def doublet(network, waveform_generator):
    return FisherMatrix(network, waveform_generator, PARAMETERS, order=2)


def test_dali_tensor_symmetries(doublet):
    """G is symmetric in its first pair; H is symmetric within and across pairs."""
    numpy.testing.assert_allclose(
        doublet.doublet, doublet.doublet.transpose(1, 0, 2), rtol=1e-8
    )
    numpy.testing.assert_allclose(
        doublet.quartic, doublet.quartic.transpose(1, 0, 2, 3), rtol=1e-8
    )
    numpy.testing.assert_allclose(
        doublet.quartic, doublet.quartic.transpose(2, 3, 0, 1), rtol=1e-8
    )


def test_quartic_term_is_negative_definite(doublet):
    """H_abcd d^a d^b d^c d^d is a squared norm, so the doublet is normalisable.

    The quartic contribution to the log-likelihood is ``-1/8`` of this, which
    must therefore be negative for every displacement -- otherwise the posterior
    would grow without bound.
    """
    count = len(doublet.names)
    flat = doublet.quartic.reshape(count * count, count * count)
    eigenvalues = numpy.linalg.eigvalsh(0.5 * (flat + flat.T))
    assert eigenvalues.min() > -1e-8 * numpy.abs(eigenvalues).max()

    generator = numpy.random.default_rng(7)
    offsets = generator.normal(size=(50, count)) * numpy.sqrt(
        numpy.diag(doublet.covariance)
    )
    outer = (offsets[:, :, None] * offsets[:, None, :]).reshape(len(offsets), -1)
    assert numpy.all(((outer @ flat) * outer).sum(axis=1) >= 0)


def test_dali_reduces_to_the_fisher_gaussian_at_order_one(doublet):
    likelihood = DALILikelihood(doublet, order=1)
    generator = numpy.random.default_rng(3)
    offsets = generator.normal(size=(5, len(doublet.names))) * numpy.sqrt(
        numpy.diag(doublet.covariance)
    )
    expected = -0.5 * numpy.einsum(
        "wa,ab,wb->w", offsets, doublet.matrix, offsets, optimize=True
    )
    numpy.testing.assert_allclose(
        likelihood.log_likelihood_from_offsets(offsets), expected, rtol=1e-10
    )


def test_dali_order_two_differs_and_peaks_at_the_fiducial_point(doublet):
    likelihood = DALILikelihood(doublet, order=2)
    assert likelihood.log_likelihood_from_offsets(
        numpy.zeros((1, len(doublet.names)))
    ) == pytest.approx(0.0)
    generator = numpy.random.default_rng(11)
    offsets = generator.normal(size=(20, len(doublet.names))) * numpy.sqrt(
        numpy.diag(doublet.covariance)
    )
    values = likelihood.log_likelihood_from_offsets(offsets)
    assert numpy.all(values < 0)
    gaussian = DALILikelihood(doublet, order=1).log_likelihood_from_offsets(offsets)
    assert not numpy.allclose(values, gaussian, rtol=1e-6)


def test_scalar_and_vectorised_likelihoods_agree(doublet):
    """The bilby entry point and the vectorised one must be the same function.

    The offsets are recovered from the absolute values rather than used
    directly, so both paths see bit-identical numbers. That matters because
    ``geocent_time`` is a GPS number of order 1e9 and a round trip through it
    loses precision -- quantified in
    :func:`test_absolute_geocent_time_limits_likelihood_precision`.
    """
    likelihood = DALILikelihood(doublet, order=2)
    absolute = {
        name: doublet.parameters[name]
        + 0.3 * numpy.sqrt(doublet.covariance[index, index])
        for index, name in enumerate(doublet.names)
    }
    offsets = numpy.array(
        [absolute[name] - doublet.parameters[name] for name in doublet.names]
    )
    assert likelihood.log_likelihood(absolute) == pytest.approx(
        float(likelihood.log_likelihood_from_offsets(offsets[numpy.newaxis, :])[0]),
        rel=1e-12,
    )


def test_absolute_geocent_time_limits_likelihood_precision(doublet):
    """Sampling t_c in absolute GPS quantises it, and this bounds how finely.

    bilby hands the sampler absolute parameter values, so ``geocent_time``
    carries the resolution of a float near 1e9 -- about 240 ns. Compared with a
    typical sigma of order 1e-4 s that leaves hundreds of resolvable steps per
    sigma, which is ample. This test pins the ratio so that a future network
    tight enough to make it marginal shows up as a failure rather than as
    quietly grainy posteriors.
    """
    index = doublet.names.index("geocent_time")
    sigma = numpy.sqrt(doublet.covariance[index, index])
    resolution = numpy.spacing(doublet.parameters["geocent_time"])
    steps_per_sigma = sigma / resolution
    assert steps_per_sigma > 100, (
        "geocent_time is resolved to only {:.0f} steps per sigma; the DALI "
        "posterior would be visibly quantised".format(steps_per_sigma)
    )


def test_prior_file_must_cover_every_fisher_parameter(tmp_path, doublet):
    """bilby's own prior files use different names, so this must fail loudly."""
    path = tmp_path / "incomplete.prior"
    path.write_text("chirp_mass = Uniform(minimum=27, maximum=29, name='chirp_mass')\n")
    with pytest.raises(KeyError, match="does not define a sampling prior"):
        build_priors(str(path), doublet.names, doublet.parameters)


def test_prior_must_contain_the_fiducial_point(tmp_path):
    path = tmp_path / "offset.prior"
    path.write_text("chirp_mass = Uniform(minimum=1, maximum=2, name='chirp_mass')\n")
    with pytest.raises(ValueError, match="lie outside the prior"):
        build_priors(str(path), ["chirp_mass"], {"chirp_mass": 28.0})


def test_enforce_physicality_clips_to_the_physical_range(tmp_path):
    path = tmp_path / "wide.prior"
    path.write_text(
        "symmetric_mass_ratio = Uniform(minimum=0.1, maximum=0.9, "
        "name='symmetric_mass_ratio')\n"
    )
    names, fiducial = ["symmetric_mass_ratio"], {"symmetric_mass_ratio": 0.24}
    loose = build_priors(str(path), names, fiducial, enforce_physicality=False)
    assert loose["symmetric_mass_ratio"].maximum == 0.9
    strict = build_priors(str(path), names, fiducial, enforce_physicality=True)
    assert strict["symmetric_mass_ratio"].maximum == 0.25


def test_narrow_prior_is_warned_about(tmp_path, caplog):
    path = tmp_path / "narrow.prior"
    path.write_text(
        "chirp_mass = Uniform(minimum=27.99, maximum=28.01, name='chirp_mass')\n"
    )
    with caplog.at_level("WARNING"):
        build_priors(
            str(path), ["chirp_mass"], {"chirp_mass": 28.0}, sigmas={"chirp_mass": 1.0}
        )
    assert "truncate the posterior" in caplog.text


# ---------------------------------------------------------------------------
# Bookkeeping
# ---------------------------------------------------------------------------


def test_missing_parameters_are_reported(network, waveform_generator):
    with pytest.raises(KeyError, match="missing from the injection"):
        WaveformDerivatives(
            network,
            waveform_generator,
            {"chirp_mass": 28.0},
            fisher_parameters=["chirp_mass", "symmetric_mass_ratio"],
        )


def test_analytic_parameters_cost_no_extra_waveform_evaluations(
    network, waveform_generator
):
    """Adding all five exact parameters must not call the waveform again."""
    numerical = ["chirp_mass", "symmetric_mass_ratio"]
    without = WaveformDerivatives(
        network, waveform_generator, PARAMETERS, fisher_parameters=numerical, order=1
    )
    with_analytic = WaveformDerivatives(
        network,
        waveform_generator,
        PARAMETERS,
        fisher_parameters=numerical + list(ANALYTIC_PARAMETERS),
        order=1,
    )
    assert (
        with_analytic.waveform_evaluations == without.waveform_evaluations
    ), "the closed-form parameters should be free"


# ---------------------------------------------------------------------------
# Redundant parametrisations
# ---------------------------------------------------------------------------
#
# A gwforge_population catalogue carries every mass parametrisation at once, and
# bilby's conversion prefers the component masses over chirp_mass. Feeding one
# straight to a Fisher over (chirp_mass, symmetric_mass_ratio) therefore
# displaced keys the waveform ignored: the mass derivatives came back *exactly*
# zero and sigma(chirp_mass) was wrong by six orders of magnitude, with nothing
# raised. These pin the fix.


@pytest.mark.parametrize(
    "varied, shadowed",
    [
        (("chirp_mass", "symmetric_mass_ratio"), ("mass_1", "mass_2", "mass_ratio")),
        (("chirp_mass", "symmetric_mass_ratio"), ("mass_1_source", "mass_2_source")),
        (("mass_1", "mass_2"), ("chirp_mass", "symmetric_mass_ratio")),
        (("chi_1",), ("a_1", "tilt_1", "spin_1z")),
    ],
)
def test_shadowed_parameters_are_removed(varied, shadowed):
    """Whatever duplicates a varied parameter is dropped; the varied ones stay."""
    parameters = dict(
        mass_1=30.0,
        mass_2=25.0,
        mass_1_source=25.0,
        mass_2_source=20.8,
        chirp_mass=23.8,
        symmetric_mass_ratio=0.2469,
        mass_ratio=0.8333,
        total_mass=55.0,
        redshift=0.2,
        luminosity_distance=1000.0,
        chi_1=0.3,
        a_1=0.3,
        tilt_1=0.0,
        spin_1z=0.3,
        ra=1.375,
    )
    stripped, removed = strip_shadowed_parameters(parameters, varied)
    for name in shadowed:
        assert name in removed
        assert name not in stripped
    for name in varied:
        assert stripped[name] == parameters[name]
    # Untouched groups are left alone.
    assert stripped["ra"] == parameters["ra"]
    assert stripped["luminosity_distance"] == parameters["luminosity_distance"]


def test_a_catalogue_row_still_constrains_the_masses(network, waveform_generator):
    """The regression itself: a Fisher built from a full catalogue row.

    Every redundant name a population file carries is present. Without the
    strip, ``convert_to_lal_binary_black_hole_parameters`` rebuilds mass_1 and
    mass_2 from the component or source masses, the chirp-mass derivative is
    identically zero, and sigma(chirp_mass) runs away to ~1e5 instead of ~1e-2.
    """
    row = dict(PARAMETERS)
    row.update(
        mass_1=32.06,
        mass_2=32.06,
        mass_ratio=1.0,
        total_mass=64.12,
        redshift=0.3,
        mass_1_source=24.66,
        mass_2_source=24.66,
        spin_1z=0.3,
        spin_2z=-0.2,
    )
    fisher = FisherMatrix(
        network,
        waveform_generator,
        row,
        fisher_parameters=["chirp_mass", "symmetric_mass_ratio", "luminosity_distance"],
    )
    assert "mass_1" not in fisher.parameters
    assert "mass_1_source" not in fisher.parameters
    assert fisher.parameters["chirp_mass"] == PARAMETERS["chirp_mass"]

    sigma = fisher.sigmas["chirp_mass"] / PARAMETERS["chirp_mass"]
    # A CE40+CE20 source at 1.5 Gpc is loud; the point is only that the mass is
    # measured at all, which the shadowed version failed by ~1e6.
    assert 0.0 < sigma < 1e-2
