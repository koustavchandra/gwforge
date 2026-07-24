import json
import numpy
import pytest

from GWForge.population.mass import Mass
from GWForge.population import _smoothed_mass as sm
from GWForge.population.user_defined import (
    get_xx_yy_from_population_priors,
    interped_prior_from_file,
)


# --------------------------------------------------------------------------- #
# UserDefined (tabulated xx/yy) model
# --------------------------------------------------------------------------- #
@pytest.fixture
def priors_file(tmp_path):
    m1 = numpy.linspace(5, 80, 400)
    m2 = numpy.linspace(3, 60, 300)
    q = numpy.linspace(0.05, 1, 200)
    data = {
        "mass_1_source": {"xx": m1.tolist(), "yy": (m1**-2.3).tolist()},
        "mass_2_source": {
            "xx": m2.tolist(),
            "yy": numpy.exp(-0.5 * ((m2 - 20) / 6) ** 2).tolist(),
        },
        "mass_ratio": {"xx": q.tolist(), "yy": q.tolist()},
    }
    path = tmp_path / "priors.json"
    path.write_text(json.dumps(data))
    return str(path)


def test_get_xx_yy_roundtrip(priors_file):
    xx, yy = get_xx_yy_from_population_priors(priors_file, "mass_1_source")
    assert xx.shape == yy.shape and xx[0] == 5 and xx[-1] == 80
    with pytest.raises(KeyError):
        get_xx_yy_from_population_priors(priors_file, "not_a_parameter")


def test_interped_prior_rejects_negative_density(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"x": {"xx": [1, 2, 3], "yy": [0.1, -0.2, 0.3]}}))
    with pytest.raises(ValueError):
        interped_prior_from_file(str(path), "x")


def test_userdefined_primary_and_mass_ratio(priors_file):
    s = Mass(
        "UserDefined",
        4000,
        {
            "file": priors_file,
            "primary_parameter": "mass_1_source",
            "mass_ratio_parameter": "mass_ratio",
        },
    ).sample()
    assert {"mass_1_source", "mass_2_source", "mass_ratio"} <= set(s)
    assert (s["mass_1_source"] >= s["mass_2_source"]).all()
    assert (
        s["mass_1_source"].min() >= 5 - 1e-6 and s["mass_1_source"].max() <= 80 + 1e-6
    )


def test_userdefined_primary_and_independent_secondary(priors_file):
    s = Mass(
        "UserDefined",
        4000,
        {
            "file": priors_file,
            "primary_parameter": "mass_1_source",
            "secondary_parameter": "mass_2_source",
        },
    ).sample()
    assert (s["mass_1_source"] >= s["mass_2_source"]).all()


def test_userdefined_requires_secondary_spec(priors_file):
    with pytest.raises(ValueError):
        Mass(
            "UserDefined",
            10,
            {"file": priors_file, "primary_parameter": "mass_1_source"},
        ).sample()


# --------------------------------------------------------------------------- #
# BGP (Broken Power Law + 2 Peaks) model
# --------------------------------------------------------------------------- #
BGP_PARAMS = dict(
    alpha_1=1.6,
    alpha_2=5.0,
    m_break=38.0,
    mmin=5.0,
    m_high=100.0,
    lam_0=0.9,
    lam_1=0.05,
    mpp_1=33.0,
    sigpp_1=4.0,
    mpp_2=10.0,
    sigpp_2=1.5,
    delta_m=4.8,
    beta=1.1,
)


def test_bgp_primary_is_normalisable_density():
    model = sm.BrokenPowerLawTwoPeakSmoothedMassDistribution(
        mmax=200, normalization_shape=(1000, 1000)
    )
    p = model.p_m1(
        {"mass_1": model.m1s}, **{k: v for k, v in BGP_PARAMS.items() if k != "beta"}
    )
    assert numpy.isfinite(p).all() and (p >= 0).all()
    assert numpy.isclose(numpy.trapezoid(p, model.m1s), 1.0, rtol=1e-3)


def test_bgp_broken_power_law_slopes():
    # Below the break the (declining) slope is alpha_1; above it is alpha_2 > alpha_1,
    # so the log-slope must steepen across m_break.
    m = numpy.array([10.0, 20.0, 50.0, 70.0])
    p = sm.broken_power_law(
        m, alpha_1=1.6, alpha_2=5.0, m_break=38.0, mmin=5.0, m_high=100.0
    )
    log_slope_low = numpy.log(p[1] / p[0]) / numpy.log(m[1] / m[0])
    log_slope_high = numpy.log(p[3] / p[2]) / numpy.log(m[3] / m[2])
    assert numpy.isclose(log_slope_low, -1.6, atol=1e-6)
    assert numpy.isclose(log_slope_high, -5.0, atol=1e-6)


def test_bgp_end_to_end_sampling():
    from bilby.core.utils import random as bilby_random

    numpy.random.seed(5)
    bilby_random.seed(5)
    s = Mass("BGP", 20000, BGP_PARAMS).sample()
    m1, m2 = s["mass_1_source"], s["mass_2_source"]
    assert (m1 >= m2).all()
    assert m1.max() <= 200
    # The high-mass Gaussian peak (mu=33) should enhance the density near 33 above
    # a locally smooth power-law expectation.
    near_peak = numpy.mean((m1 > 31) & (m1 < 35))
    off_peak = numpy.mean((m1 > 44) & (m1 < 48))
    assert near_peak > off_peak
