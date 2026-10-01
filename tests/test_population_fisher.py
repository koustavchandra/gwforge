"""Tests for :mod:`GWForge.population_fisher`.

Three claims are being defended.

**The densities are GWForge's own.** The mass model here calls
:mod:`GWForge.population._smoothed_mass` rather than reimplementing it, and the
tests check that what comes out is the same distribution
:class:`GWForge.population.mass.Mass` samples from. A population Fisher run on a
catalogue drawn from a *different* density than the one being differentiated is
not wrong by a little; it is meaningless.

**The scores are analytic and exact.** Every model's ``analytic_score`` is
checked against a centred finite difference. The comparison is scaled by the
column's root-mean-square rather than by each entry, because a score column
legitimately passes through zero and a per-entry relative error there says
nothing.

**The estimator is right, not merely self-consistent.** ``test_fisher_matches_
maximum_likelihood_scatter`` draws many independent catalogues, fits each by
maximum likelihood, and checks the Term-I sigma against the actual spread of the
estimates. Every other test in this file would pass if the whole formalism were
subtly wrong.
"""

import configparser
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy
import pytest
from scipy import stats
from scipy.optimize import minimize

from GWForge.population._smoothed_mass import (
    BrokenPowerLawTwoPeakSmoothedMassDistribution,
)
from GWForge.population.mass import Mass
from GWForge.population_fisher import (
    BrokenPowerLawTwoPeakMass,
    DefaultSpin,
    JointPopulationModel,
    MadauDickinsonRedshift,
    PopulationFisherResult,
    SpectralSirenModel,
    fit_model,
    load_catalogue,
    population_fisher,
)
from GWForge.population_fisher import base
from GWForge.population_fisher.catalogue import expand_detectors, network_snr_from
from GWForge.population_fisher.config import build_model
from GWForge.population.mass import BGP_PARAMETERS
from GWForge.population.spin import DEFAULT_BBH_SPIN_PARAMETERS, truncated_normal
from GWForge.cosmology import FlatwCDM
from GWForge.plotting import POPULATION_LATEX_LABELS

MAXIMUM_REDSHIFT = 10.0

# The BGP fiducial the shipped configs use: O4b Default BBH posterior medians.
#
# Peak *1* is the narrow low-mass one (9.99 +/- 0.66 Msun) carrying 52% of the
# weight; peak 2 is the broad 33.3 Msun one carrying 5.7%. The primary and
# secondary have independent tapers -- (4.49, 3.12) against (3.49, 5.60) -- so
# these values exercise the whole model rather than a degenerate corner of it.
MASS_PARAMETERS = dict(
    alpha_1=1.456442737,
    alpha_2=5.100400428,
    m_break=37.9806382,
    mmin=4.489078237,
    delta_m=3.123416382,
    mmin_2=3.487749145,
    delta_m_2=5.60056759,
    m_high=300.0,
    lam_0=0.4206178692,
    lam_1=0.5228817702,
    mpp_1=9.989384122,
    sigpp_1=0.6601839103,
    mpp_2=33.2656028,
    sigpp_2=4.582221363,
    beta=0.8049660633,
    maximum_mass=300,
)

# The nine parameters the forecasts vary. ``mmin``, ``m_high``,
# ``maximum_mass`` are support, and ``m_break``/``delta_m`` are pinned; see the
# shipped configs.
MASS_FREE_PARAMETERS = [
    "alpha_1",
    "alpha_2",
    "lam_0",
    "lam_1",
    "mpp_1",
    "sigpp_1",
    "mpp_2",
    "sigpp_2",
    "beta",
]


@pytest.fixture(scope="module")
def mass_model():
    return BrokenPowerLawTwoPeakMass(
        mmin=MASS_PARAMETERS["mmin"],
        m_high=MASS_PARAMETERS["m_high"],
        maximum_mass=MASS_PARAMETERS["maximum_mass"],
        mmin_2=MASS_PARAMETERS["mmin_2"],
        delta_m_2=MASS_PARAMETERS["delta_m_2"],
        fiducial={
            name: value
            for name, value in MASS_PARAMETERS.items()
            if name
            not in ("mmin", "m_high", "maximum_mass", "mmin_2", "delta_m_2")
        },
    )


@pytest.fixture(scope="module")
def redshift_model():
    return MadauDickinsonRedshift(maximum_redshift=MAXIMUM_REDSHIFT)


@pytest.fixture(scope="module")
def spin_model():
    return DefaultSpin()


@pytest.fixture(scope="module")
def events():
    """A synthetic catalogue spanning each model's support."""
    generator = numpy.random.default_rng(20260822)
    count = 400
    low = MASS_PARAMETERS["mmin"] + 0.5
    secondary_low = MASS_PARAMETERS["mmin_2"] + 0.5
    primary = generator.uniform(low + 1.0, 150.0, count)
    secondary = secondary_low + generator.uniform(0.02, 0.95, count) * (
        primary - secondary_low
    )
    return {
        "mass_1_source": primary,
        "mass_2_source": secondary,
        "mass_ratio": secondary / primary,
        "redshift": generator.uniform(0.02, 9.5, count),
        "a_1": generator.uniform(0.02, 0.95, count),
        "a_2": generator.uniform(0.02, 0.95, count),
        "tilt_1": generator.uniform(0.0, numpy.pi, count),
        "tilt_2": generator.uniform(0.0, numpy.pi, count),
    }


def assert_score_matches_finite_difference(model, events, parameters=None, rtol=1e-4):
    """Compare an analytic score with a finite difference, column by column.

    Scaled by each column's root-mean-square: a score column passes through zero
    for individual events, and a per-entry relative error there is meaningless
    even when the derivative is exact to fourteen digits.

    The tolerance is set by the *finite difference*, not by the analytic score.
    ``m_break`` is the clearest case: it is a kink in p(m_1), so for an event
    near the break the two stencil points straddle it and the difference is a
    chord across a corner. At ``h_rel = 1e-4`` that shows up as 4e-5 of the
    column scale; shrink the step to 2.5e-5, so that no event's stencil crosses
    the kink, and it collapses to 4e-10.
    """
    analytic = model.score(events, parameters, method="analytic")
    numerical = model.score(events, parameters, method="finite-difference")
    scale = numpy.sqrt(numpy.mean(analytic**2, axis=0))
    deviation = numpy.nanmax(numpy.abs(analytic - numerical), axis=0) / scale
    worst = int(numpy.argmax(deviation))
    assert deviation[worst] < rtol, "{} column disagrees by {:.2e} of its scale".format(
        model.parameter_names[worst], deviation[worst]
    )


# ---------------------------------------------------------------------------
# The densities are the ones GWForge samples from
# ---------------------------------------------------------------------------


def test_primary_mass_matches_smoothed_mass(mass_model):
    """p(m_1) against ``_smoothed_mass``, which is what ``Mass`` samples from.

    The comparison is of *shape*, then of normalisation, separately. Both codes
    normalise by a trapezoidal integral, but on different grids: this module
    clusters its nodes towards ``mmin`` so as to resolve the taper, while
    ``_smoothed_mass`` uses a uniform one. With ``delta_m`` this small the
    integrand has an effectively hard edge, where the trapezoid converges only
    linearly, so the two constants differ at the 1e-4 level however fine either
    grid is. Folding that into a single ``assert_allclose`` would either hide a
    real shape difference behind a loose tolerance or fail on a quadrature
    detail. ``test_mass_density_normalised`` is what pins the normalisation.
    """
    reference = BrokenPowerLawTwoPeakSmoothedMassDistribution(
        mmax=200, normalization_shape=(20000, 1000)
    )
    keywords = {
        name: value
        for name, value in MASS_PARAMETERS.items()
        if name not in ("beta", "maximum_mass", "mmin_2", "delta_m_2")
    }
    grid = numpy.linspace(MASS_PARAMETERS["mmin"] + 0.5, 150.0, 500)
    expected = reference.p_m1({"mass_1": grid}, **keywords)

    normalisation = numpy.trapezoid(
        mass_model._primary_numerator(mass_model.grid, mass_model.fiducial),
        mass_model.grid,
    )
    ours = mass_model._primary_numerator(grid, mass_model.fiducial) / normalisation

    keep = expected > 0
    ratio = ours[keep] / expected[keep]
    # Same function: the ratio is one constant, to machine precision.
    assert ratio.std() / ratio.mean() < 1e-12
    # Same normalisation, to the accuracy the two quadratures share.
    assert ratio.mean() == pytest.approx(1.0, rel=1e-3)


@pytest.mark.parametrize("primary", (20.0, 60.0, 120.0))
def test_mass_ratio_matches_smoothed_mass(mass_model, primary):
    """p(q | m_1) against ``_smoothed_mass``.

    Ours is a density in ``m_2`` and ``_smoothed_mass``'s is in ``q``, so they
    differ by the Jacobian ``m_1``. That factor carries no hyper-parameter, so
    it cancels from every score, but it has to be right for ``log_prob`` -- and
    for the maximum-likelihood fits, which compare log densities.
    """
    reference = BrokenPowerLawTwoPeakSmoothedMassDistribution(
        mmax=200, normalization_shape=(20000, 2000)
    )
    ratio = numpy.linspace(1.05 * MASS_PARAMETERS["mmin_2"] / primary, 1.0, 400)
    masses = numpy.full_like(ratio, primary)
    expected = reference.p_q(
        {"mass_1": masses, "mass_ratio": ratio},
        beta=MASS_PARAMETERS["beta"],
        mmin=MASS_PARAMETERS["mmin"],
        delta_m=MASS_PARAMETERS["delta_m"],
        mmin_2=MASS_PARAMETERS["mmin_2"],
        delta_m_2=MASS_PARAMETERS["delta_m_2"],
    )
    normalisation = numpy.trapezoid(
        mass_model._primary_numerator(mass_model.grid, mass_model.fiducial),
        mass_model.grid,
    )
    joint = numpy.exp(
        mass_model.log_prob(
            {"mass_1_source": masses, "mass_ratio": ratio}, mass_model.fiducial
        )
    )
    marginal = mass_model._primary_numerator(masses, mass_model.fiducial) / normalisation
    ours = joint / marginal * primary

    keep = expected > 1e-10
    fraction = ours[keep] / expected[keep]
    # Same conditional, up to the two codes' different quadratures for its
    # normalisation -- see test_primary_mass_matches_smoothed_mass.
    assert fraction.std() / fraction.mean() < 1e-10
    assert fraction.mean() == pytest.approx(1.0, rel=1e-3)


def test_mass_density_normalised(mass_model):
    """p(m_1, m_2) integrates to one over the triangle m_min <= m_2 <= m_1."""
    # The secondary's support starts at its own taper edge, which is below the
    # primary's; integrating from mmin would miss part of the distribution.
    mmin = MASS_PARAMETERS["mmin"]
    mmin_2 = MASS_PARAMETERS["mmin_2"]
    primary_grid = numpy.linspace(mmin + 1e-3, 299.9, 900)
    fraction = numpy.linspace(0.0, 1.0, 800)
    inner = []
    for primary in primary_grid:
        secondary = mmin_2 + fraction * (primary - mmin_2)
        values = mass_model.log_prob(
            {
                "mass_1_source": numpy.full_like(secondary, primary),
                "mass_2_source": secondary,
            },
            mass_model.fiducial,
        )
        density = numpy.where(numpy.isfinite(values), numpy.exp(values), 0.0)
        inner.append(numpy.trapezoid(density, secondary))
    assert numpy.trapezoid(numpy.array(inner), primary_grid) == pytest.approx(
        1.0, rel=1e-3
    )


def test_redshift_density_normalised(redshift_model):
    grid = numpy.linspace(1e-6, MAXIMUM_REDSHIFT, 20000)
    density = numpy.exp(redshift_model.log_prob({"redshift": grid}, redshift_model.fiducial))
    assert numpy.trapezoid(density, grid) == pytest.approx(1.0, rel=1e-6)


def test_spin_density_matches_scipy(spin_model):
    """Eqs. B15-B16 against scipy's own pdfs.

    ``Spin`` only samples the Default model, so unlike the mass and redshift
    blocks this density is written out rather than reused -- which makes an
    independent check of it worth more, not less.
    """
    parameters = spin_model.fiducial

    # B15: truncated Gaussian magnitudes on [0, 1].
    magnitudes = numpy.linspace(1e-6, 1.0 - 1e-6, 5000)
    terms = truncated_normal(
        magnitudes, parameters["mu_chi"], parameters["sigma_chi"], 0.0, 1.0
    )
    numpy.testing.assert_allclose(
        terms["density"],
        stats.truncnorm.pdf(
            magnitudes,
            (0.0 - parameters["mu_chi"]) / parameters["sigma_chi"],
            (1.0 - parameters["mu_chi"]) / parameters["sigma_chi"],
            loc=parameters["mu_chi"],
            scale=parameters["sigma_chi"],
        ),
        rtol=1e-10,
    )

    # B16 marginalised over the partner: the mixture, with a *free* mean.
    cosines = numpy.linspace(-1.0, 1.0, 5001)
    gaussian = truncated_normal(
        cosines, parameters["mu_t"], parameters["sigma_t"], -1.0, 1.0
    )["density"]
    expected = stats.truncnorm.pdf(
        cosines,
        (-1.0 - parameters["mu_t"]) / parameters["sigma_t"],
        (1.0 - parameters["mu_t"]) / parameters["sigma_t"],
        loc=parameters["mu_t"],
        scale=parameters["sigma_t"],
    )
    numpy.testing.assert_allclose(gaussian, expected, rtol=1e-10)

    # The joint tilt density integrates to one over the square.
    grid = numpy.linspace(-1.0, 1.0, 1500)
    first, second = numpy.meshgrid(grid, grid, indexing="ij")
    joint = parameters["xi_spin"] * truncated_normal(
        first, parameters["mu_t"], parameters["sigma_t"], -1.0, 1.0
    )["density"] * truncated_normal(
        second, parameters["mu_t"], parameters["sigma_t"], -1.0, 1.0
    )["density"] + (1.0 - parameters["xi_spin"]) * 0.25
    assert numpy.trapezoid(
        numpy.trapezoid(joint, grid, axis=1), grid
    ) == pytest.approx(1.0, rel=1e-4)


# ---------------------------------------------------------------------------
# The scores are analytic and exact
# ---------------------------------------------------------------------------


def test_mass_score_matches_finite_difference(mass_model, events):
    assert_score_matches_finite_difference(mass_model, events)


def test_redshift_score_matches_finite_difference(redshift_model, events):
    assert_score_matches_finite_difference(redshift_model, events)


def test_spin_score_matches_finite_difference(spin_model, events):
    assert_score_matches_finite_difference(spin_model, events)


def test_scores_match_off_the_fiducial(mass_model, redshift_model, events):
    """Exactness must not be an accident of evaluating at the default point."""
    shifted = dict(mass_model.fiducial)
    shifted.update(alpha_1=1.4, m_break=31.0, lam_1=0.42, sigpp_2=5.1, delta_m=1.7)
    assert_score_matches_finite_difference(mass_model, events, shifted)
    assert_score_matches_finite_difference(
        redshift_model, events, dict(gamma=3.4, kappa=4.1, z_peak=2.4)
    )


def test_joint_score_is_the_blocks_side_by_side(
    mass_model, redshift_model, spin_model, events
):
    joint = JointPopulationModel([mass_model, redshift_model, spin_model])
    expected = numpy.hstack(
        [
            mass_model.score(events, method="analytic"),
            redshift_model.score(events, method="analytic"),
            spin_model.score(events, method="analytic"),
        ]
    )
    numpy.testing.assert_allclose(joint.score(events, method="analytic"), expected)
    numpy.testing.assert_allclose(
        joint.log_prob(events, joint.fiducial),
        mass_model.log_prob(events, mass_model.fiducial)
        + redshift_model.log_prob(events, redshift_model.fiducial)
        + spin_model.log_prob(events, spin_model.fiducial),
    )


def test_joint_rejects_duplicate_parameters(mass_model):
    with pytest.raises(ValueError, match="more than one sub-model"):
        JointPopulationModel([mass_model, mass_model])


def test_joint_fiducial_follows_a_refit(redshift_model, spin_model):
    """A block refitted after the joint model is built must be seen by it.

    The redshift block is *always* refitted (the generator convolves a time
    delay), so a fiducial snapshotted in ``__init__`` would quietly send the
    Fisher back to the pre-fit point.
    """
    joint = JointPopulationModel([redshift_model, spin_model])
    original = redshift_model.fiducial["gamma"]
    try:
        redshift_model.fiducial["gamma"] = 1.234
        assert joint.fiducial["gamma"] == 1.234
    finally:
        redshift_model.fiducial["gamma"] = original


def test_mass_derivatives_with_respect_to_the_masses(mass_model, events):
    """d ln p / d m, the chain the spectral siren rides on."""
    step = 1e-5
    primary, secondary = mass_model.dlog_prob_dmasses(events, mass_model.fiducial)
    for key, analytic in (("mass_1_source", primary), ("mass_2_source", secondary)):
        above = dict(events, **{key: events[key] + step})
        below = dict(events, **{key: events[key] - step})
        numerical = (
            mass_model.log_prob(above, mass_model.fiducial)
            - mass_model.log_prob(below, mass_model.fiducial)
        ) / (2.0 * step)
        scale = numpy.sqrt(numpy.mean(analytic**2))
        assert numpy.nanmax(numpy.abs(analytic - numerical)) / scale < 1e-3


def test_spectral_siren_score_matches_finite_difference(
    mass_model, redshift_model, spin_model
):
    """All 21 columns, including the cosmology ones, which are the point.

    The previous generation of this code finite-differenced the ``Om0`` and
    ``w0`` columns against a spline rebuilt per trial cosmology, and could not
    shrink the step past 1e-4 without the answer getting worse. These are closed
    form.
    """
    cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=-1.0)
    redshifts = MadauDickinsonRedshift(
        maximum_redshift=MAXIMUM_REDSHIFT, cosmology=cosmology
    )
    model = SpectralSirenModel(mass_model, redshifts, spin_model, cosmology=cosmology)

    generator = numpy.random.default_rng(4)
    count = 200
    redshift = generator.uniform(0.05, 8.0, count)
    primary = generator.uniform(6.0, 90.0, count)
    secondary = 5.5 + generator.uniform(0.05, 0.9, count) * (primary - 5.5)
    detector_events = {
        "mass_1": primary * (1.0 + redshift),
        "mass_2": secondary * (1.0 + redshift),
        "luminosity_distance": cosmology.luminosity_distance(redshift),
        "a_1": generator.uniform(0.02, 0.95, count),
        "a_2": generator.uniform(0.02, 0.95, count),
        "tilt_1": generator.uniform(0.0, numpy.pi, count),
        "tilt_2": generator.uniform(0.0, numpy.pi, count),
    }
    original = base.FD_RELATIVE_STEP
    try:
        base.FD_RELATIVE_STEP = 1e-5
        assert_score_matches_finite_difference(model, detector_events, rtol=1e-4)
    finally:
        base.FD_RELATIVE_STEP = original


class _PlainGaussian(base.PopulationModel):
    """A stand-in for the kind of model the finite-difference tier exists to serve.

    Deliberately given no ``analytic_score``, so ``score`` has to fall back. Its
    true score is elementary, which is what makes it a usable check on the
    machinery -- the tier has to be checked before a model that needs it exists.
    """

    parameter_names = ["mu", "sigma"]
    fiducial = {"mu": 1.0, "sigma": 2.0}
    event_keys = ["x"]

    def log_prob(self, events, parameters):
        scaled = (events["x"] - parameters["mu"]) / parameters["sigma"]
        return (
            -0.5 * scaled**2
            - numpy.log(parameters["sigma"])
            - 0.5 * numpy.log(2.0 * numpy.pi)
        )


def test_finite_difference_fallback_matches_the_exact_score():
    """The tier-2 fallback, on a model whose derivatives are known exactly."""
    model = _PlainGaussian()
    sample = {"x": numpy.linspace(-4.0, 6.0, 50)}
    scaled = (sample["x"] - model.fiducial["mu"]) / model.fiducial["sigma"]
    expected = numpy.column_stack(
        [
            scaled / model.fiducial["sigma"],
            (scaled**2 - 1.0) / model.fiducial["sigma"],
        ]
    )
    numpy.testing.assert_allclose(
        model.score(sample, method="finite-difference"), expected, rtol=1e-6
    )
    # And "auto" reaches it, since the model has no analytic score.
    numpy.testing.assert_allclose(model.score(sample), expected, rtol=1e-6)


def test_auto_falls_back_when_there_is_no_analytic_score(redshift_model, events):
    """A missing analytic score is not an error under "auto", only under an explicit ask."""

    class _NoAnalyticScore(MadauDickinsonRedshift):
        def analytic_score(self, events, parameters):
            return None

    model = _NoAnalyticScore(maximum_redshift=MAXIMUM_REDSHIFT)
    fallback = model.score(events)
    numpy.testing.assert_allclose(
        fallback, redshift_model.score(events, method="analytic"), rtol=1e-4
    )
    with pytest.raises(NotImplementedError):
        model.score(events, method="analytic")


def test_unknown_score_method_is_rejected(redshift_model, events):
    """``autodiff`` used to be a tier; asking for it now names the two that exist."""
    with pytest.raises(ValueError, match="method must be one of"):
        redshift_model.score(events, method="autodiff")


# ---------------------------------------------------------------------------
# The block combinations the shipped configs offer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "blocks",
    [
        ["redshift"],
        ["mass", "redshift"],
        ["mass", "redshift", "spin"],
        ["mass", "redshift", "cosmology"],
    ],
)
def test_every_documented_block_combination_runs(blocks, events):
    """The four combinations the ini files offer, end to end.

    Individual models were well covered; the *combinations* were not, which is
    how ``build_model`` came to pass an argument ``DefaultSpin`` never took and
    break every config naming the spin block.
    """
    import configparser

    from GWForge.population_fisher.config import build_model

    config = configparser.ConfigParser()
    config.add_section("Model")
    config.set("Model", "blocks", repr(blocks))
    config.set("Model", "maximum-redshift", str(MAXIMUM_REDSHIFT))
    model, built, sub_models, _ = build_model(config)
    assert built == blocks
    for name in blocks:
        if name != "cosmology":
            assert name in sub_models

    sample = dict(events)
    if "spin" in blocks:
        count = len(sample["redshift"])
        sample["a_1"] = numpy.full(count, 0.3)
        sample["a_2"] = numpy.full(count, 0.4)
        sample["cos_tilt_1"] = numpy.full(count, 0.5)
        sample["cos_tilt_2"] = numpy.full(count, -0.2)
    if "cosmology" in blocks:
        # The spectral-siren model is a density over detector-frame observables.
        one_plus_z = 1.0 + sample["redshift"]
        sample = {
            "mass_1": sample["mass_1_source"] * one_plus_z,
            "mass_2": sample["mass_2_source"] * one_plus_z,
            "luminosity_distance": 1e3 * one_plus_z,
        }

    scores = model.score(sample, method="analytic")
    assert scores.shape == (len(next(iter(sample.values()))), len(model.parameter_names))
    assert numpy.isfinite(scores).any()


def test_cosmology_alone_is_refused():
    """``cosmology`` is not a block of its own; it needs mass and redshift."""
    import configparser

    from GWForge.population_fisher.config import build_model

    config = configparser.ConfigParser()
    config.add_section("Model")
    config.set("Model", "blocks", "['cosmology']")
    with pytest.raises(ValueError):
        build_model(config)


# ---------------------------------------------------------------------------
# The estimator itself
# ---------------------------------------------------------------------------


def test_fisher_matches_maximum_likelihood_scatter(redshift_model):
    """Term I against the actual spread of maximum-likelihood estimates.

    With no selection every event is "detected", so this is the plain asymptotic
    MLE covariance and the Fisher must reproduce it. This is the one test that
    would fail if the formalism -- rather than its implementation -- were wrong.
    """
    grid = numpy.linspace(1e-4, MAXIMUM_REDSHIFT, 200000)
    density = numpy.exp(redshift_model.log_prob({"redshift": grid}, redshift_model.fiducial))
    cumulative = numpy.cumsum(density)
    cumulative /= cumulative[-1]
    generator = numpy.random.default_rng(7)
    names = redshift_model.parameter_names

    def draw(count):
        return numpy.interp(generator.random(count), cumulative, grid)

    count = 4000
    forecast = population_fisher(redshift_model, {"redshift": draw(count)})

    def fit(redshifts):
        sample = {"redshift": redshifts}
        return minimize(
            lambda vector: -numpy.sum(
                redshift_model.log_prob(sample, dict(zip(names, vector)))
            ),
            redshift_model.to_vector(redshift_model.fiducial),
            jac=lambda vector: -numpy.sum(
                redshift_model.score(sample, dict(zip(names, vector)), method="analytic"),
                axis=0,
            ),
            method="L-BFGS-B",
            bounds=[(0.1, 8.0), (0.5, 15.0), (0.3, 6.0)],
        ).x

    estimates = numpy.array([fit(draw(count)) for _ in range(120)])
    scatter = estimates.std(axis=0, ddof=1)
    for index, name in enumerate(names):
        # 120 realisations pin a standard deviation to about 6%, so 15% leaves
        # room for Monte-Carlo noise without leaving room for a real error.
        assert forecast.sigma[name] == pytest.approx(scatter[index], rel=0.15)


def test_fisher_is_symmetric_and_positive_semidefinite(redshift_model, events):
    result = population_fisher(redshift_model, events)
    numpy.testing.assert_allclose(result.fisher, result.fisher.T)
    assert numpy.linalg.eigvalsh(result.fisher).min() >= -1e-8 * numpy.abs(
        result.fisher
    ).max()


def test_fisher_scales_with_the_number_of_events(redshift_model, events):
    """Gamma is a sum over events, so duplicating the catalogue doubles it."""
    single = population_fisher(redshift_model, events)
    doubled = population_fisher(
        redshift_model, {"redshift": numpy.concatenate([events["redshift"]] * 2)}
    )
    numpy.testing.assert_allclose(doubled.fisher, 2.0 * single.fisher, rtol=1e-10)
    for name in single.parameter_names:
        assert doubled.sigma[name] == pytest.approx(
            single.sigma[name] / numpy.sqrt(2.0), rel=1e-8
        )
    rescaled = single.scaled_to(2 * single.n_events)
    numpy.testing.assert_allclose(rescaled.fisher, doubled.fisher, rtol=1e-10)


def test_pinned_parameters_leave_the_matrix(redshift_model, events):
    result = population_fisher(
        redshift_model, events, fixed_parameters={"z_peak": 1.9}
    )
    assert result.parameter_names == ["gamma", "kappa"]
    assert result.fisher.shape == (2, 2)
    assert result.fiducial["z_peak"] == 1.9


def test_free_and_fixed_cannot_overlap(redshift_model, events):
    with pytest.raises(ValueError, match="both free and fixed"):
        population_fisher(
            redshift_model,
            events,
            free_parameters=["gamma"],
            fixed_parameters={"gamma": 2.7},
        )


def test_unknown_free_parameter_is_rejected(redshift_model, events):
    with pytest.raises(ValueError, match="not a hyper-parameter"):
        population_fisher(redshift_model, events, free_parameters=["H0"])


def test_events_outside_the_support_are_dropped_not_zeroed(redshift_model, events):
    """An event the model says cannot exist must not enter the sum with weight one."""
    polluted = {
        "redshift": numpy.concatenate(
            [events["redshift"], [MAXIMUM_REDSHIFT + 5.0, -1.0]]
        )
    }
    result = population_fisher(redshift_model, polluted)
    assert result.n_events == len(events["redshift"])
    assert result.n_total == len(polluted["redshift"])


def test_degenerate_direction_is_named(mass_model, events, caplog):
    """A parameter with an event-independent score is unmeasurable in Term I.

    Turning the taper off makes ``m_1``'s low-mass edge a hard cutoff, and the
    score of a hard cutoff is the same number for every event in the interior,
    so its centred score vanishes identically. That is a fact about the model,
    and the right response is to name it rather than return a plausible sigma.
    """
    hard = dict(mass_model.fiducial, delta_m=0.0)
    result = population_fisher(
        mass_model,
        events,
        free_parameters=["alpha_1", "delta_m"],
        parameters=hard,
    )
    assert "delta_m" in result.degenerate


# ---------------------------------------------------------------------------
# Fits, catalogues, and the sampler they have to agree with
# ---------------------------------------------------------------------------


def test_fit_recovers_the_generating_parameters(mass_model):
    """A round trip through ``Mass``'s sampler and back out through the fit.

    Not a tautology: ``Mass`` samples by inverting a CDF built on
    ``_smoothed_mass``'s grids, and the fit maximises this package's own log
    density on a different grid entirely. Agreement means the two describe the
    same distribution, which is the whole basis for taking a fitted value as the
    fiducial of a forecast.

    This is the sharpest test in the file. The fiducial has a narrow (sigma =
    0.89 M_sun) peak sitting almost on top of a hard low-mass edge, so it
    catches a primary grid too coarse to resolve the peak, a mass-ratio grid
    that does not start at the support edge, and a peak mapped to the wrong
    index -- all of which were real bugs.
    """
    from bilby.core.utils import random as bilby_random

    # Both RNGs: the primary mass comes from a bilby Interped prior and the
    # conditional mass ratio from numpy, so seeding only one leaves the result
    # dependent on whatever test ran before this one.
    numpy.random.seed(3)
    bilby_random.seed(3)
    samples = Mass(
        mass_model="BGP", number_of_samples=40000, parameters=dict(MASS_PARAMETERS)
    ).sample()
    primary = samples["mass_1_source"]
    events = {
        "mass_1_source": primary,
        "mass_2_source": primary * samples["mass_ratio"],
    }
    fitted = fit_model(mass_model, events, free_parameters=MASS_FREE_PARAMETERS)
    # Not one sample outside the support: the sampler and the density agree on
    # where the population lives, edges included.
    assert fitted["_n_dropped"] == 0
    for name in MASS_FREE_PARAMETERS:
        assert fitted[name] == pytest.approx(MASS_PARAMETERS[name], rel=0.05), name


def test_sampled_secondaries_respect_the_minimum_mass():
    """Regression: ``mass_ratio`` must be drawn conditionally on ``mass_1``.

    ``Mass`` used to build the mass-ratio prior from ``p_q`` evaluated with the
    ``m1s`` and ``qs`` grids paired *by index* -- the diagonal of the conditional
    -- and then sample ``q`` independently of ``m_1``. Because the conditional's
    support is ``q >= mmin / m_1``, that produced secondaries below ``mmin``,
    outside the model's own support, and got the shape wrong besides. Every
    log density in this package is then ``-inf`` for those events.
    """
    numpy.random.seed(11)
    samples = Mass(
        mass_model="BGP", number_of_samples=20000, parameters=dict(MASS_PARAMETERS)
    ).sample()
    secondary = samples["mass_1_source"] * samples["mass_ratio"]
    # The secondary's support starts at *its* taper edge, which sits below the
    # primary's -- so this is mmin_2, not mmin.
    assert secondary.min() >= MASS_PARAMETERS["mmin_2"]
    assert numpy.all(samples["mass_ratio"] <= 1.0)


def test_fit_moves_off_a_start_outside_a_box_constraint(spin_model):
    """Regression: the optimiser must not stop dead at its starting point.

    A trial step that leaves a model's feasible region makes the objective
    non-finite. Returning ``inf`` there aborts L-BFGS-B's line search and it
    reports *convergence* -- the fit hands back its starting values, which reads
    as a perfect recovery rather than as a failure. A large finite penalty lets
    it backtrack instead.

    Started deliberately 15% away, the fit has to come back.
    """
    from bilby.core.utils import random as bilby_random
    from GWForge.population.spin import Spin

    # Seed both RNGs: the magnitudes come from a bilby prior and the mixture
    # split from numpy, so leaving either to the ambient state makes this test
    # pass or fail depending on what ran before it.
    numpy.random.seed(17)
    bilby_random.seed(17)
    samples = Spin(
        spin_model="Default", number_of_samples=20000, parameters=dict(spin_model.fiducial)
    ).sample()
    events = {key: samples[key] for key in ("a_1", "a_2", "tilt_1", "tilt_2")}
    truth = dict(spin_model.fiducial)
    displaced = {name: 1.15 * value for name, value in truth.items()}
    from_truth = fit_model(spin_model, events, start=truth)
    from_displaced = fit_model(spin_model, events, start=displaced)

    # The property being tested is that the optimiser *finds the maximum*, not
    # that it finds the true value. Two different starts must land in the same
    # place; a stuck optimiser returns each start unchanged and so gives two
    # different answers.
    for name in truth:
        assert from_displaced[name] != displaced[name], "{} never moved".format(name)
        scale = max(abs(from_truth[name]), 1e-3)
        assert abs(from_displaced[name] - from_truth[name]) < 0.02 * scale, name

    # It is also unbiased, but only in expectation: at 20,000 samples the
    # scatter on xi_spin alone is 0.095, comparable to the 15% displacement, so
    # a single realisation says nothing tighter than this.
    assert from_truth["sigma_chi"] == pytest.approx(truth["sigma_chi"], rel=0.05)


def test_conditional_sampling_stays_reproducible():
    """One seed must still govern the whole draw.

    The conditional mass-ratio sampler draws its uniforms from numpy's global
    RNG while the primary mass comes from a bilby prior. ``gwforge_population``
    seeds both, so a seeded run has to be repeatable -- which is worth pinning,
    because it is the kind of thing a refactor breaks silently.
    """
    from bilby.core.utils import random as bilby_random

    draws = []
    for _ in range(2):
        numpy.random.seed(42)
        bilby_random.seed(42)
        samples = Mass(
            mass_model="BGP", number_of_samples=500, parameters=dict(MASS_PARAMETERS)
        ).sample()
        draws.append((samples["mass_1_source"], samples["mass_ratio"]))
    numpy.testing.assert_array_equal(draws[0][0], draws[1][0])
    numpy.testing.assert_array_equal(draws[0][1], draws[1][1])


def test_result_round_trips_through_npz(redshift_model, events, tmp_path):
    """A saved forecast must be re-plottable without recomputing it.

    Comparing two networks means loading two archives; if this drifts, the
    comparison silently plots stale or wrong numbers.
    """
    result = population_fisher(redshift_model, events)
    path = tmp_path / "forecast.npz"
    numpy.savez(
        path,
        fisher=result.fisher,
        covariance=result.covariance,
        parameter_names=numpy.array(result.parameter_names, dtype=object),
        fiducial=numpy.array(
            [result.fiducial[name] for name in result.parameter_names]
        ),
        snr_threshold=10.0,
        n_events=result.n_events,
        n_total=result.n_total,
        condition_number=result.condition_number,
    )
    reloaded = PopulationFisherResult.from_npz(str(path))

    assert reloaded.parameter_names == result.parameter_names
    numpy.testing.assert_allclose(reloaded.fisher, result.fisher)
    numpy.testing.assert_allclose(reloaded.covariance, result.covariance)
    assert reloaded.n_events == result.n_events
    assert reloaded.condition_number == pytest.approx(result.condition_number)
    for name in result.parameter_names:
        assert reloaded.sigma[name] == pytest.approx(result.sigma[name])
        assert reloaded.fiducial[name] == pytest.approx(result.fiducial[name])
    # The per-event scores are deliberately not saved -- one row per detected
    # event would dominate the archive -- so they come back empty.
    assert reloaded.scores.shape == (0, len(result.parameter_names))


def test_expand_detectors_names_the_triangle_arms():
    assert expand_detectors(["CE40", "ET"]) == ["CE40", "ET1", "ET2", "ET3"]


def test_network_snr_adds_in_quadrature():
    data = {
        "CE40_optimal_snr": numpy.array([3.0, 6.0]),
        "CE20_optimal_snr": numpy.array([4.0, 8.0]),
    }
    snr, detectors = network_snr_from(data)
    numpy.testing.assert_allclose(snr, [5.0, 10.0])
    assert detectors == ["CE20", "CE40"]
    with pytest.raises(KeyError, match="No SNR column"):
        network_snr_from(data, ["ET1"])


def _write_catalogue(directory, rows=8, shift=0.0):
    """A minimal population/SNR file pair, as the two CLIs would write them."""
    import h5py

    population = directory / "population.h5"
    snr = directory / "snr.h5"
    primary = numpy.linspace(10.0, 60.0, rows)
    with h5py.File(population, "w") as handle:
        handle["mass_1_source"] = primary
        handle["mass_2_source"] = 0.8 * primary
        handle["mass_1"] = primary * 1.5
        handle["redshift"] = numpy.linspace(0.1, 2.0, rows)
        handle["luminosity_distance"] = numpy.linspace(500.0, 15000.0, rows)
    with h5py.File(snr, "w") as handle:
        handle["index"] = numpy.arange(rows)
        handle["mass_1"] = primary * 1.5 + shift
        handle["luminosity_distance"] = numpy.linspace(500.0, 15000.0, rows)
        handle["CE40_optimal_snr"] = numpy.linspace(5.0, 40.0, rows)
    return population, snr


def test_catalogue_joins_and_thresholds(tmp_path):
    population, snr = _write_catalogue(tmp_path)
    catalogue = load_catalogue(str(population), str(snr), snr_threshold=20.0)
    assert catalogue.n_total == 8
    assert catalogue.n_detected == int((numpy.linspace(5.0, 40.0, 8) >= 20.0).sum())
    assert catalogue.network_snr.min() >= 20.0
    assert "mass_ratio" in catalogue.events


def test_catalogue_refuses_a_mismatched_snr_file(tmp_path):
    """The join is checked, not assumed: forecasting for the wrong sources is silent."""
    population, snr = _write_catalogue(tmp_path, shift=5.0)
    with pytest.raises(ValueError, match="does not belong to this population"):
        load_catalogue(str(population), str(snr), snr_threshold=1.0)


# ---------------------------------------------------------------------------
# Defaults, shipped configs and the CLI: the forecast side of what
# ``test_population_defaults.py`` pins on the generation side. These live here
# rather than there so the committed suite does not import this package.
# ---------------------------------------------------------------------------

FIDUCIAL = DEFAULT_BBH_SPIN_PARAMETERS
PACKAGE_DIRECTORY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "GWForge"
)
CONFIGURATION_DIRECTORY = os.path.join(
    PACKAGE_DIRECTORY, "population_fisher", "configuration_files"
)
BIN = Path(__file__).resolve().parent.parent / "bin"


def parse(path):
    config = configparser.ConfigParser()
    config.read(path)
    return config


def test_fisher_mass_defaults_are_the_medians():
    model = BrokenPowerLawTwoPeakMass()
    assert model.fiducial == {
        name: BGP_PARAMETERS[name] for name in model.parameter_names
    }
    # Support parameters, which are fixed by construction rather than fitted.
    assert model.mmin == BGP_PARAMETERS["mmin"]
    assert model.mmin_2 == BGP_PARAMETERS["mmin_2"]
    assert model.delta_m_2 == BGP_PARAMETERS["delta_m_2"]
    assert model.m_high == BGP_PARAMETERS["m_high"]
    assert model.maximum_mass == BGP_PARAMETERS["maximum_mass"]


def test_fisher_spin_defaults_are_the_medians():
    model = DefaultSpin()
    assert model.fiducial == {
        name: DEFAULT_BBH_SPIN_PARAMETERS[name] for name in model.parameter_names
    }

# Every shipped Fisher analysis config, by filename. Read from the directory
# rather than listed, so a config added later is covered without anyone
# remembering to add it here.
ANALYSIS_CONFIGURATIONS = sorted(
    name for name in os.listdir(CONFIGURATION_DIRECTORY) if name.endswith(".ini")
)


@pytest.mark.parametrize("name", ANALYSIS_CONFIGURATIONS)
def test_analysis_configs_use_the_current_spin_names(name):
    """``sigma_squared_chi`` is a *variance* and belongs to the Beta models.

    ``DefaultSpin`` takes ``sigma_chi``, a standard deviation, and a free
    ``mu_t``; ``_fiducial`` rejects anything else, so these configs used to
    fail to build.
    """
    path = os.path.join(CONFIGURATION_DIRECTORY, name)
    assert "sigma_squared_chi" not in open(path).read()
    config = parse(path)
    spin = json.loads(config.get("Model", "spin-parameters").replace("'", '"'))
    assert spin == DEFAULT_BBH_SPIN_PARAMETERS

# ---------------------------------------------------------------------------
# Regressions in the Fisher config path
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ANALYSIS_CONFIGURATIONS)
def test_every_shipped_analysis_config_builds(name):
    """Each shipped example assembles a model whose free parameters exist.

    These are what the docs point users at, so a config naming a parameter the
    model does not have is a broken first experience rather than a subtle bug.
    """
    config = parse(os.path.join(CONFIGURATION_DIRECTORY, name))
    model, blocks, _, _ = build_model(config)
    assert blocks
    free = json.loads(
        config.get("Fisher", "free-parameters", fallback="[]").replace("'", '"')
    )
    unknown = [name for name in free if name not in model.parameter_names]
    assert not unknown, "{} frees {}, which the model does not have".format(
        name, unknown
    )


def test_a_config_with_a_spin_block_builds():
    """``build_model`` passed ``minimum_spin=``, which ``DefaultSpin`` never took.

    Every config naming the ``spin`` block raised ``TypeError`` at construction.
    """
    config = parse(os.path.join(CONFIGURATION_DIRECTORY, "mass_redshift_spin.ini"))
    _, blocks, sub_models, _ = build_model(config)
    assert "spin" in blocks
    assert sub_models["spin"].maximum_spin == 1.0
    assert sub_models["spin"].minimum_cosine_tilt == -1.0
    assert sub_models["spin"].fiducial == DEFAULT_BBH_SPIN_PARAMETERS


def test_the_built_mass_model_keeps_the_secondary_taper_independent():
    """The config used to override the class defaults and drop mmin_2/delta_m_2.

    That silently collapsed the secondary's taper onto the primary's, undoing
    the model fix on every Fisher run.
    """
    config = configparser.ConfigParser()
    config.add_section("Model")
    config.set("Model", "blocks", "['mass']")
    _, _, sub_models, _ = build_model(config)
    mass = sub_models["mass"]
    assert mass.mmin == BGP_PARAMETERS["mmin"]
    assert mass.mmin_2 == BGP_PARAMETERS["mmin_2"] != mass.mmin
    assert mass.delta_m_2 == BGP_PARAMETERS["delta_m_2"] != mass.fiducial["delta_m"]
    assert mass.m_high == BGP_PARAMETERS["m_high"]

def test_every_hyper_parameter_has_a_label():
    """A new parameter without a label prints a raw key with a mangled underscore."""
    models = [
        BrokenPowerLawTwoPeakMass(),
        MadauDickinsonRedshift(maximum_redshift=10.0),
        DefaultSpin(),
    ]
    missing = [
        name
        for model in models
        for name in model.parameter_names
        if name not in POPULATION_LATEX_LABELS
    ]
    assert not missing, "no LaTeX label for {}".format(missing)

def test_spin_densities_agree_with_the_fisher_model():
    """The forecast and the density checks must describe one distribution."""
    from GWForge.population.spin import (
        default_spin_magnitude_density,
        default_spin_tilt_density,
        truncated_normal,
    )
    from GWForge.population_fisher import DefaultSpin

    model = DefaultSpin()
    magnitudes = numpy.linspace(0.05, 0.95, 40)
    cosines = numpy.linspace(-0.9, 0.9, 40)
    events = {
        "a_1": magnitudes,
        "a_2": magnitudes[::-1],
        "cos_tilt_1": cosines,
        "cos_tilt_2": cosines[::-1],
    }
    gaussian = [
        truncated_normal(cosines, FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], -1.0, 1.0)[
            "density"
        ],
        truncated_normal(
            cosines[::-1], FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], -1.0, 1.0
        )["density"],
    ]
    expected = (
        numpy.log(
            default_spin_magnitude_density(
                magnitudes, FIDUCIAL["mu_chi"], FIDUCIAL["sigma_chi"]
            )
        )
        + numpy.log(
            default_spin_magnitude_density(
                magnitudes[::-1], FIDUCIAL["mu_chi"], FIDUCIAL["sigma_chi"]
            )
        )
        + numpy.log(
            FIDUCIAL["xi_spin"] * gaussian[0] * gaussian[1]
            + (1.0 - FIDUCIAL["xi_spin"]) / 4.0
        )
    )
    numpy.testing.assert_allclose(
        model.log_prob(events, model.fiducial), expected, rtol=1e-12
    )
    # And the marginal really is the partner integrated out.
    assert default_spin_tilt_density(
        cosines, FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], FIDUCIAL["xi_spin"]
    ).shape == cosines.shape

def test_population_fisher_compare_writes_a_figure(tmp_path):
    """``--compare`` overlays saved forecasts without recomputing them.

    Exercised end to end because it is the path that produces the
    network-versus-network figure, and it has no other test: the plotting is
    covered in ``test_plotting.py`` and the Fisher in
    ``test_population_fisher.py``, but nothing else runs the two together
    through the CLI.
    """
    import numpy

    names = numpy.array(["gamma", "kappa"], dtype=object)
    paths = []
    for index, scale in enumerate((1.0, 4.0)):
        path = tmp_path / "forecast_{}.npz".format(index)
        numpy.savez(
            path,
            fisher=numpy.eye(2) / scale**2,
            covariance=numpy.eye(2) * scale**2,
            parameter_names=names,
            fiducial=numpy.array([2.7, 5.6]),
            snr_threshold=10.0,
            n_events=100,
            n_total=1000,
            condition_number=1.0,
        )
        paths.append(path)

    output = tmp_path / "comparison.pdf"
    result = subprocess.run(
        [
            sys.executable,
            str(BIN / "gwforge_population_fisher"),
            "--compare",
            "A={}".format(paths[0]),
            "B={}".format(paths[1]),
            "--output-file",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output.exists()
    # The ratio column is what the figure is for; B is four times wider.
    assert "4" in result.stdout
