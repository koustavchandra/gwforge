"""The population mass models: BGP, the smoothed family, and UserDefined.

Three claims are made here, in widening order of how much machinery they touch:

1. **Agreement with** ``gwpopulation``. GWForge reimplements the ``gwpopulation``
   smoothed mass distributions in pure numpy (``GWForge/population/_smoothed_mass.py``)
   so that ``gwpopulation`` is not a runtime dependency. It stays a *test-only*
   dependency here, and the reimplementations are asserted to match it to
   numerical precision, so the dependency cannot be dropped by silently changing
   the physics. These three tests skip without it; everything else still runs.
2. **Self-contained properties.** Normalisation, the smoothing window's
   endpoints, and the broken power law's slopes depend on no reference
   implementation at all.
3. **End-to-end sampling.** ``Mass(...).sample()`` for BGP and UserDefined,
   including the secondary's independent taper.
"""

import json

import numpy
import pytest

from GWForge.population import _smoothed_mass as sm
from GWForge.population.mass import BGP_PARAMETERS, Mass
from GWForge.population.user_defined import (
    get_xx_yy_from_population_priors,
    interped_prior_from_file,
)

SHAPE = (1000, 1000)

# The GWTC-5.0 medians, imported rather than copied so these exercise the shape
# the package actually ships. The pre-O4b set these replaced had the peaks the
# other way round -- peak 1 the broad 33 Msun one -- which quietly contradicted
# every config in the tree.
BGP_PARAMS = dict(BGP_PARAMETERS)

# The keys that describe the mass ratio or size the grid, not the primary mass.
# Mass.sample strips exactly these before calling p_m1.
NOT_PRIMARY = ("beta", "maximum_mass", "mmin_2", "delta_m_2")


def _bgp_model(**overrides):
    """The BGP distribution at the shipped medians, plus the parameters used."""
    parameters = dict(BGP_PARAMS, **overrides)
    model = sm.BrokenPowerLawTwoPeakSmoothedMassDistribution(
        mmin=parameters["mmin"],
        mmax=parameters["maximum_mass"],
        normalization_shape=(4000, 1000),
    )
    return model, parameters


# ---------------------------------------------------------------------------
# 1. Agreement with gwpopulation
# ---------------------------------------------------------------------------

MASS_CASES = [
    (
        "PowerLawPeak",
        "SinglePeakSmoothedMassDistribution",
        dict(
            alpha=3.37,
            delta_m=5.23,
            mmin=4.89,
            mmax=88.81,
            lam=0.04,
            mpp=33.60,
            sigpp=4.59,
        ),
        dict(beta=2.0, mmin=4.89, delta_m=5.23),
    ),
    (
        "MultiPeak",
        "MultiPeakSmoothedMassDistribution",
        dict(
            alpha=2.9,
            delta_m=4.8,
            mmin=5.0,
            mmax=90.0,
            lam=0.1,
            lam_1=0.5,
            mpp_1=10.0,
            sigpp_1=2.0,
            mpp_2=35.0,
            sigpp_2=5.0,
        ),
        dict(beta=1.0, mmin=5.0, delta_m=4.8),
    ),
    (
        "BrokenPowerLaw",
        "BrokenPowerLawSmoothedMassDistribution",
        dict(
            alpha_1=1.6,
            alpha_2=5.6,
            delta_m=4.8,
            mmin=5.0,
            mmax=90.0,
            break_fraction=0.4,
        ),
        dict(beta=1.0, mmin=5.0, delta_m=4.8),
    ),
]


@pytest.mark.parametrize("name,class_name,m1_params,q_params", MASS_CASES)
def test_smoothed_mass_matches_gwpopulation(name, class_name, m1_params, q_params):
    gwpop_mass = pytest.importorskip("gwpopulation.models.mass")

    new = getattr(sm, class_name)(normalization_shape=SHAPE)
    gw = getattr(gwpop_mass, class_name)(normalization_shape=SHAPE)

    assert numpy.array_equal(new.m1s, gw.m1s)
    assert numpy.array_equal(new.qs, gw.qs)

    new_pm1 = new.p_m1({"mass_1": new.m1s}, **dict(m1_params))
    gw_pm1 = gw.p_m1({"mass_1": gw.m1s}, **dict(m1_params))
    rel_pm1 = numpy.max(numpy.abs(new_pm1 - gw_pm1)) / numpy.max(numpy.abs(gw_pm1))
    assert rel_pm1 < 1e-6, "{} p_m1 mismatch: {}".format(name, rel_pm1)

    new_pq = numpy.asarray(
        new.p_q({"mass_ratio": new.qs, "mass_1": new.m1s}, **q_params)
    )
    gw_pq = numpy.asarray(gw.p_q({"mass_ratio": gw.qs, "mass_1": gw.m1s}, **q_params))
    denom = numpy.max(numpy.abs(gw_pq))
    rel_pq = numpy.max(numpy.abs(new_pq - gw_pq)) / denom if denom > 0 else 0.0
    assert rel_pq < 1e-6, "{} p_q mismatch: {}".format(name, rel_pq)


# ---------------------------------------------------------------------------
# 2. Self-contained properties
# ---------------------------------------------------------------------------


def test_smoothing_window_endpoints():
    # Below mmin -> 0; well above mmin+delta_m and below mmax -> 1.
    m = numpy.array([1.0, 5.0, 50.0, 200.0])
    w = sm.smoothing(m, mmin=5.0, mmax=100.0, delta_m=5.0)
    assert w[0] == 0.0  # below mmin
    assert numpy.isclose(w[2], 1.0)  # in the flat region
    assert w[3] == 0.0  # above mmax
    # delta_m == 0 disables the taper entirely.
    assert numpy.array_equal(
        sm.smoothing(m, mmin=5.0, mmax=100.0, delta_m=0.0), numpy.ones(4)
    )


def test_powerlaw_and_truncnorm_normalised():
    x = numpy.linspace(2, 100, 200001)
    p = sm.powerlaw(x, alpha=-2.0, high=100, low=2)
    assert numpy.isclose(numpy.trapezoid(p, x), 1.0, rtol=1e-4)
    g = sm.truncnorm(x, mu=35.0, sigma=5.0, high=100, low=2)
    assert numpy.isclose(numpy.trapezoid(g, x), 1.0, rtol=1e-4)


def test_bgp_primary_is_normalisable_density():
    model, parameters = _bgp_model()
    p = model.p_m1(
        {"mass_1": model.m1s},
        **{k: v for k, v in parameters.items() if k not in NOT_PRIMARY},
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


# ---------------------------------------------------------------------------
# The secondary's independent taper (GWTC-5.0 Tab. 5)
# ---------------------------------------------------------------------------


def test_secondary_taper_is_independent_of_the_primary():
    """GWTC-5.0 Tab. 5 gives m2 its own ``(m2,low, delta_m,2)``.

    The paper's prior is ``m2,low ~ U(3, m1,low)``, so the secondary generally
    turns on *below* the primary, and the O4b medians separate the two pairs
    clearly (4.49/3.12 against 3.49/5.60). Sharing one pair -- which is what
    ``p_q`` used to do -- is a different model, and it shows up as a secondary
    distribution that switches on in the wrong place.
    """
    model, parameters = _bgp_model()
    primary = numpy.full(4000, 60.0)
    ratio = numpy.linspace(1e-3, 1.0, 4000)
    dataset = {"mass_1": primary, "mass_ratio": ratio}

    separate = model.p_q(
        dataset,
        beta=parameters["beta"],
        mmin=parameters["mmin"],
        delta_m=parameters["delta_m"],
        mmin_2=parameters["mmin_2"],
        delta_m_2=parameters["delta_m_2"],
    )
    shared = model.p_q(
        dataset,
        beta=parameters["beta"],
        mmin=parameters["mmin"],
        delta_m=parameters["delta_m"],
    )
    secondary = primary * ratio
    # The independent taper turns on at m2,low, below the primary's edge.
    assert secondary[separate > 0].min() < parameters["mmin"]
    assert secondary[shared > 0].min() >= parameters["mmin"] * 0.999
    # Both are still normalised conditionals.
    for density in (separate, shared):
        assert numpy.trapezoid(density, ratio) == pytest.approx(1.0, rel=1e-3)


def test_secondary_taper_defaults_to_the_primary():
    """Omitting the new arguments must reproduce the old single-taper model.

    This is what keeps ``PowerLaw+Peak``, ``MultiPeak`` and ``BrokenPowerLaw``
    -- which have no separate secondary edge -- bit-for-bit unchanged.
    """
    model, parameters = _bgp_model()
    dataset = {
        "mass_1": numpy.full(2000, 40.0),
        "mass_ratio": numpy.linspace(1e-3, 1.0, 2000),
    }
    common = dict(
        beta=parameters["beta"], mmin=parameters["mmin"], delta_m=parameters["delta_m"]
    )
    numpy.testing.assert_allclose(
        model.p_q(dataset, **common),
        model.p_q(
            dataset,
            mmin_2=parameters["mmin"],
            delta_m_2=parameters["delta_m"],
            **common,
        ),
        rtol=0,
        atol=0,
    )


# ---------------------------------------------------------------------------
# 3. End-to-end sampling
# ---------------------------------------------------------------------------


def test_bgp_end_to_end_sampling(seed):
    seed(5)
    s = Mass("BGP", 20000, dict(BGP_PARAMS)).sample()
    m1, m2 = s["mass_1_source"], s["mass_2_source"]
    assert (m1 >= m2).all()
    assert m1.max() <= BGP_PARAMS["maximum_mass"]
    assert m1.min() >= BGP_PARAMS["mmin"]
    assert m2.min() >= BGP_PARAMS["mmin_2"]
    # The broad Gaussian peak (mu = 33.3) should enhance the density near 33
    # above a locally smooth power-law expectation.
    near_peak = numpy.mean((m1 > 31) & (m1 < 35))
    off_peak = numpy.mean((m1 > 44) & (m1 < 48))
    assert near_peak > off_peak


def test_sampled_secondaries_follow_their_own_taper(seed):
    """End to end: the drawn ``m2`` turns on at ``mmin_2`` over ``delta_m_2``."""
    seed(4)
    samples = Mass(
        mass_model="BGP", number_of_samples=40000, parameters=dict(BGP_PARAMS)
    ).sample()
    secondary = samples["mass_1_source"] * samples["mass_ratio"]
    assert secondary.min() >= BGP_PARAMS["mmin_2"]
    # Secondaries below the *primary's* edge exist, and would not if the taper
    # were shared.
    assert (secondary < BGP_PARAMS["mmin"]).sum() > 0


def test_mass_fixed_model():
    # "fixed" previously raised KeyError via an unconditional re-sample.
    m = Mass(
        mass_model="fixed",
        number_of_samples=20,
        parameters={"primary_mass": 30.0, "mass_ratio": 0.8},
    )
    samples = m.sample()
    assert numpy.allclose(samples["mass_1_source"], 30.0)
    assert len(samples["mass_2_source"]) == 20


def test_mass_dipbreak_highmass_branch_nonzero():
    # The high-mass power law used maximum=gamma_high, making its probability
    # identically zero above gamma_high. After the fix, samples must appear above it.
    params = {
        "mmin": 1.0,
        "mmax": 100.0,
        "alpha_1": 1.5,
        "alpha_2": 2.0,
        "gamma_low": 3.0,
        "gamma_high": 5.0,
        "eta_low": 50.0,
        "eta_high": 50.0,
        "A": 0.5,
        "n": 50.0,
    }
    m = Mass(mass_model="powerlawdipbreak", number_of_samples=2000, parameters=params)
    samples = m.sample()
    assert numpy.max(samples["mass_1_source"]) > params["gamma_high"]


# ---------------------------------------------------------------------------
# UserDefined (tabulated xx/yy) model
# ---------------------------------------------------------------------------


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
