"""Tests for the GWTC-5.0 Default BBH spin model.

``GWForge.population.spin.Spin`` with ``spin-model = Default`` implements
Eqs. B15-B16 of arXiv:2605.27226. It used to be an alias for the GWTC-3-era
Beta model, and these tests pin the three ways that differed, each of which is
invisible to a check that only looks at one component's marginal:

* the magnitude is a **truncated Gaussian**, not a Beta, and ``sigma_chi`` is a
  standard deviation rather than a variance;
* the tilt Gaussian's mean :math:`\\mu_t` is a **free parameter**, not pinned at
  ``cos t = 1``;
* the tilt mixture is over the **binary** -- one Bernoulli decides the component
  for both tilts -- so the two cosines are correlated. Factorising it leaves
  both marginals exactly right and sets that correlation to zero, which is why
  ``test_tilts_are_correlated_across_the_binary`` compares against the analytic
  correlation rather than eyeballing a histogram.
"""

import numpy
import pytest
from scipy import stats

from GWForge.population.spin import (
    DEFAULT_BBH_SPIN_PARAMETERS as FIDUCIAL,
    Spin,
)

# Parameters chosen to make the tilt correlation large and unmistakable. At the
# fiducial the correlation is ~0.005, below the sampling noise of any catalogue
# we generate, so a discriminating test has to move mu_t towards alignment.
STRONG_CORRELATION = dict(mu_chi=0.2, sigma_chi=0.3, mu_t=1.0, sigma_t=0.3, xi_spin=0.5)


def draw(parameters, count, seed):
    from bilby.core.utils import random as bilby_random

    # Seed both RNGs: bilby's priors and numpy's shuffles both feed the
    # sampler, so seeding one leaves the draw dependent on test order.
    numpy.random.seed(seed)
    bilby_random.seed(seed)
    return Spin(
        spin_model="Default", number_of_samples=count, parameters=dict(parameters)
    ).sample()


def truncated(mu, sigma, low, high):
    return stats.truncnorm((low - mu) / sigma, (high - mu) / sigma, loc=mu, scale=sigma)


def chi_squared(samples, distribution, low, high, bins=40):
    """Reduced chi-squared of a histogram against a distribution's CDF."""
    edges = numpy.linspace(low, high, bins + 1)
    counts, _ = numpy.histogram(samples, bins=edges)
    expected = samples.size * numpy.diff(distribution.cdf(edges))
    usable = expected >= 25
    return numpy.sum(
        (counts[usable] - expected[usable]) ** 2 / expected[usable]
    ) / usable.sum()


def test_magnitudes_are_truncated_gaussian_not_beta():
    """Eq. B15. A Beta with the same mean and variance is a different shape.

    The distinction is sharpest at the endpoints: the truncated Gaussian has
    finite non-zero density at chi = 0, while a Beta either vanishes there or
    diverges.
    """
    samples = draw(FIDUCIAL, 200000, 11)
    model = truncated(FIDUCIAL["mu_chi"], FIDUCIAL["sigma_chi"], 0.0, 1.0)
    for name in ("a_1", "a_2"):
        assert chi_squared(samples[name], model, 0.0, 1.0) < 1.5
        assert samples[name].min() >= 0.0
        assert samples[name].max() <= 1.0

    # And it is genuinely not the Beta this used to draw: match the Beta's mean
    # and variance to the truncated Gaussian's and it still disagrees.
    mean, variance = model.mean(), model.var()
    alpha = mean * (mean - mean**2 - variance) / variance
    beta = alpha * (1.0 / mean - 1.0)
    assert chi_squared(samples["a_1"], stats.beta(alpha, beta), 0.0, 1.0) > 10.0


def test_tilt_marginal_is_the_mixture_with_a_free_mean():
    """Eq. B16 marginalised over the partner, with mu_t away from 1."""
    samples = draw(FIDUCIAL, 200000, 12)
    gaussian = truncated(FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], -1.0, 1.0)
    edges = numpy.linspace(-1.0, 1.0, 41)
    expected_cdf = FIDUCIAL["xi_spin"] * gaussian.cdf(edges) + (
        1.0 - FIDUCIAL["xi_spin"]
    ) * (edges + 1.0) / 2.0
    for name in ("tilt_1", "tilt_2"):
        cosines = numpy.cos(samples[name])
        counts, _ = numpy.histogram(cosines, bins=edges)
        expected = cosines.size * numpy.diff(expected_cdf)
        reduced = numpy.sum((counts - expected) ** 2 / expected) / len(counts)
        assert reduced < 1.5, name


def test_tilt_mean_is_not_pinned_at_one():
    """Regression: mu_t was a literal 1 in the model this replaced."""
    low = numpy.cos(draw(dict(FIDUCIAL, mu_t=-0.8), 40000, 13)["tilt_1"])
    high = numpy.cos(draw(dict(FIDUCIAL, mu_t=0.9), 40000, 13)["tilt_1"])
    assert low.mean() < high.mean() - 0.2


def test_tilts_are_correlated_across_the_binary():
    """Eq. B16 is joint over the pair: one draw decides the component for both.

    The old sampler drew the split per component and shuffled the two arrays
    independently, which reproduces both marginals exactly and destroys this
    correlation. Only a joint statistic sees the difference.
    """
    parameters = STRONG_CORRELATION
    samples = draw(parameters, 200000, 14)
    first = numpy.cos(samples["tilt_1"])
    second = numpy.cos(samples["tilt_2"])

    gaussian = truncated(parameters["mu_t"], parameters["sigma_t"], -1.0, 1.0)
    mean, variance = gaussian.mean(), gaussian.var()
    xi = parameters["xi_spin"]
    # Isotropic on [-1, 1] has zero mean and variance 1/3.
    expectation = xi * mean
    second_moment = xi * (variance + mean**2) + (1.0 - xi) / 3.0
    predicted = (xi * mean**2 - expectation**2) / (second_moment - expectation**2)

    measured = numpy.corrcoef(first, second)[0, 1]
    uncertainty = 1.0 / numpy.sqrt(first.size)
    assert predicted > 0.3, "test parameters must make the two models separable"
    assert abs(measured - predicted) < 5 * uncertainty
    # ... and decisively not the factorised model, whose correlation is zero.
    assert measured > 10 * uncertainty


def test_magnitudes_stay_independent():
    """B15 is iid: only the *tilts* are correlated across the binary."""
    samples = draw(STRONG_CORRELATION, 200000, 15)
    correlation = numpy.corrcoef(samples["a_1"], samples["a_2"])[0, 1]
    assert abs(correlation) < 5 / numpy.sqrt(samples["a_1"].size)


def test_missing_parameters_name_the_rename():
    """`sigma_squared_chi` belonged to the Beta model; say so rather than KeyError."""
    with pytest.raises(ValueError, match="sigma_chi"):
        Spin(
            spin_model="Default",
            number_of_samples=10,
            parameters=dict(mu_chi=0.1, sigma_squared_chi=0.02, sigma_t=1.0, xi_spin=0.5),
        ).sample()


def test_beta_model_still_available_under_its_own_name():
    """Making `Default` the GWTC-5 model must not delete the old one."""
    samples = Spin(
        spin_model="Isotropic-Beta_Gaussian_Uniform",
        number_of_samples=5000,
        parameters=dict(mu_chi=0.26, sigma_squared_chi=0.02, sigma_t=0.87, xi_spin=0.76),
    ).sample()
    assert set(("a_1", "a_2", "tilt_1", "tilt_2")).issubset(samples)
    assert samples["a_1"].max() <= 0.99


# ---------------------------------------------------------------------------
# The densities the Fisher model and the validation script share
# ---------------------------------------------------------------------------


def test_spin_densities_are_normalised():
    """Both marginals integrate to one over their support."""
    from scipy.integrate import trapezoid

    from GWForge.population.spin import (
        default_spin_magnitude_density,
        default_spin_tilt_density,
    )

    magnitudes = numpy.linspace(0.0, 1.0, 20000)
    density = default_spin_magnitude_density(
        magnitudes, FIDUCIAL["mu_chi"], FIDUCIAL["sigma_chi"]
    )
    assert trapezoid(density, magnitudes) == pytest.approx(1.0, abs=1e-6)

    cosines = numpy.linspace(-1.0, 1.0, 20000)
    density = default_spin_tilt_density(
        cosines, FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], FIDUCIAL["xi_spin"]
    )
    assert trapezoid(density, cosines) == pytest.approx(1.0, abs=1e-6)


def test_tilt_marginal_is_not_the_joint():
    """The mixture is over the *binary*, so the pair does not factorise.

    Getting this wrong leaves both marginals untouched and silently sets
    ``corr(cos t_1, cos t_2)`` to zero, which is the bug the catalogue check in
    ``validation/`` exists to catch.
    """
    from GWForge.population.spin import (
        default_spin_tilt_density,
        truncated_normal,
    )

    cosines = numpy.array([-0.5, 0.0, 0.5, 0.9])
    marginal = default_spin_tilt_density(
        cosines, FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], FIDUCIAL["xi_spin"]
    )
    gaussian = truncated_normal(
        cosines, FIDUCIAL["mu_t"], FIDUCIAL["sigma_t"], -1.0, 1.0
    )["density"]
    joint = FIDUCIAL["xi_spin"] * gaussian**2 + (1.0 - FIDUCIAL["xi_spin"]) / 4.0
    assert not numpy.allclose(joint, marginal**2)


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


# ---------------------------------------------------------------------------
# Construction and reproducibility
# ---------------------------------------------------------------------------


def test_spin_mutable_default_does_not_leak():
    # Two instances constructed without an explicit parameters dict must not
    # share mutated state (the old mutable default argument leaked alpha/beta_chi).
    Spin(spin_model="nonspinning", number_of_samples=5)
    s2 = Spin(spin_model="nonspinning", number_of_samples=5)
    assert "minimum_primary_spin" in s2.parameters
    # A fresh instance should not carry keys from a previous one's computation.
    assert (
        "alpha_chi"
        not in Spin(spin_model="nonspinning", number_of_samples=5).parameters
    )


def test_spin_sampling_reproducible_with_seed():
    # Seeding both numpy (for the tilt shuffle) and bilby (for prior draws) must
    # make spin sampling reproducible.
    a = draw(FIDUCIAL, 200, 1)
    b = draw(FIDUCIAL, 200, 1)
    assert numpy.array_equal(a["a_1"], b["a_1"])
    assert numpy.array_equal(a["tilt_1"], b["tilt_1"])
