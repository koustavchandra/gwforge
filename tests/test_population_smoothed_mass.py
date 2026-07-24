"""Regression tests locking the in-house smoothed mass/redshift models against
``gwpopulation``.

GWForge reimplements the ``gwpopulation`` smoothed mass distributions and the
Madau-Dickinson / power-law redshift evolutions in pure numpy (see
``GWForge/population/_smoothed_mass.py`` and the ``psi_of_z`` helpers in
``GWForge/population/redshift.py``) so that ``gwpopulation`` is no longer a
runtime dependency. These tests keep ``gwpopulation`` as a *test-only*
dependency and assert the reimplementations match it to numerical precision, so
the dependency can never be dropped by silently changing the physics.

If ``gwpopulation`` is not installed the comparison tests skip (the pure-numpy
smoke tests still run).
"""

import numpy
import pytest

from GWForge.population import _smoothed_mass as sm
from GWForge.population import redshift as rz

gwpop_mass = pytest.importorskip("gwpopulation.models.mass")
gwpop_redshift = pytest.importorskip("gwpopulation.models.redshift")

SHAPE = (1000, 1000)

MASS_CASES = [
    (
        "PowerLawPeak",
        sm.SinglePeakSmoothedMassDistribution,
        gwpop_mass.SinglePeakSmoothedMassDistribution,
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
        sm.MultiPeakSmoothedMassDistribution,
        gwpop_mass.MultiPeakSmoothedMassDistribution,
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
        sm.BrokenPowerLawSmoothedMassDistribution,
        gwpop_mass.BrokenPowerLawSmoothedMassDistribution,
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


@pytest.mark.parametrize("name,new_cls,gw_cls,m1_params,q_params", MASS_CASES)
def test_smoothed_mass_matches_gwpopulation(name, new_cls, gw_cls, m1_params, q_params):
    new = new_cls(normalization_shape=SHAPE)
    gw = gw_cls(normalization_shape=SHAPE)

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


def test_madau_dickinson_psi_matches_gwpopulation():
    z = numpy.linspace(1e-6, 20, 2500)
    params = dict(gamma=2.7, kappa=5.6, z_peak=1.9)
    new = rz.madau_dickinson_psi_of_z(z, **params)
    gw = gwpop_redshift.MadauDickinsonRedshift(z_max=20).psi_of_z(redshift=z, **params)
    assert numpy.allclose(new, gw, rtol=1e-12, atol=0)


def test_power_law_psi_matches_gwpopulation():
    z = numpy.linspace(1e-6, 20, 2500)
    new = rz.power_law_psi_of_z(z, lamb=3.0)
    gw = gwpop_redshift.PowerLawRedshift(z_max=20).psi_of_z(redshift=z, lamb=3.0)
    assert numpy.allclose(new, gw, rtol=1e-12, atol=0)


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
