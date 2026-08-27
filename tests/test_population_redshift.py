"""Tests for the in-house (pycbc-free, sympy-free) redshift model.

``redshift.py`` previously convolved the star-formation rate with the
formation->merger time-delay distribution using ``pycbc.population`` plus
symbolic ``sympy`` integration. It now does this with direct vectorised
numerical integration over a log-spaced time-delay grid. These tests lock the
result against an *independent* high-accuracy adaptive-quadrature computation of
the same convolution (so they need neither pycbc nor gwpopulation), and check
the sampling contract.
"""

import numpy
import pytest
from scipy.integrate import quad
from lal import YRJUL_SI, PC_SI

from GWForge.population.redshift import (
    Redshift,
    madau_dickinson_psi_of_z,
    power_law_psi_of_z,
)

MD_PARAMS = {"gamma": 2.7, "kappa": 5.6, "z_peak": 1.9}


def _reference_rate_density(z_values, maxz=5.0):
    """High-accuracy R(z_m)/R(0) for the Madau-Dickinson + 1/tau model.

    Independent of GWForge's implementation: builds its own lookback table and
    integrates the delay convolution in log-tau with adaptive quadrature.
    """
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=30.0,
        maximum_redshift=maxz,
        gps_start_time=0,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    cosmo = r.import_cosmology()
    h0 = cosmo.H0.value / (PC_SI * 1e3) * (YRJUL_SI / 1e-9)  # Gyr^-1

    def dtdz(z):
        return 1.0 / (h0 * (1 + z) * numpy.sqrt(cosmo.Ode0 + cosmo.Om0 * (1 + z) ** 3))

    zt = numpy.linspace(0, 1000, 400000)
    lt = numpy.concatenate(
        [[0.0], numpy.cumsum((dtdz(zt[1:]) + dtdz(zt[:-1])) / 2 * numpy.diff(zt))]
    )
    age = cosmo.age(0).to("Gyr").value
    td_min = 0.02

    def sfr(z):
        return madau_dickinson_psi_of_z(z, **MD_PARAMS)

    def rate(zm):
        tm = numpy.interp(zm, zt, lt)
        return quad(
            lambda u: sfr(numpy.interp(tm + numpy.exp(u), lt, zt)),
            numpy.log(td_min),
            numpy.log(age),
            epsabs=1e-12,
            epsrel=1e-10,
            limit=500,
        )[0]

    r0 = rate(0.0)
    return numpy.array([rate(z) / r0 for z in z_values])


def test_rate_density_matches_high_accuracy_reference():
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=30.0,
        maximum_redshift=5.0,
        gps_start_time=0,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    rd = r.rate_density(elements=800)
    z = numpy.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
    mine = numpy.array([rd(zi) / rd(0.0) for zi in z])
    ref = _reference_rate_density(z)
    rel = numpy.max(numpy.abs(mine - ref) / ref)
    assert (
        rel < 1e-3
    ), "rate_density deviates from high-accuracy reference by {}".format(rel)


def test_rate_density_normalised_to_local_rate():
    local = 30.0
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=local,
        maximum_redshift=5.0,
        gps_start_time=0,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    rd = r.rate_density(elements=200)
    # R(0) must equal the local rate density (Gpc^-3 yr^-1 -> Mpc^-3 yr^-1).
    assert numpy.isclose(rd(0.0), local * 1e-9, rtol=1e-6)


def test_inverse_is_slope_one_and_unknown_model_rejected():
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=30.0,
        maximum_redshift=2.0,
        gps_start_time=0,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    assert r.time_delay_slope == 1.0
    with pytest.raises(ValueError):
        Redshift(
            redshift_model="MadauDickinson",
            local_merger_rate_density=30.0,
            maximum_redshift=2.0,
            gps_start_time=0,
            parameters=MD_PARAMS,
            time_delay_model="log_normal",
        )


def test_time_delay_probability_support_and_shape():
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=30.0,
        maximum_redshift=2.0,
        gps_start_time=0,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    age = r._maximum_time_delay()
    tau = numpy.array([0.001, 0.02, 0.1, 1.0, age * 0.99, age * 1.5])
    p = r.time_delay_probability(tau)
    assert p[0] == 0.0 and p[-1] == 0.0  # outside [td_min, age]
    # 1/tau shape inside the support
    assert numpy.allclose(p[2] * tau[2], p[3] * tau[3], rtol=1e-12)


def test_sample_contract():
    maxz = 3.0
    r = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=30.0,
        maximum_redshift=maxz,
        gps_start_time=1893024018,
        analysis_time=YRJUL_SI,
        parameters=MD_PARAMS,
        time_delay_model="inverse",
    )
    s = r.sample()
    assert set(s) == {"redshift", "time_interval", "geocent_time"}
    assert (s["redshift"] >= 0).all() and (s["redshift"] <= maxz).all()
    assert (s["time_interval"] > 0).all()
    assert (numpy.diff(s["geocent_time"]) >= 0).all()  # cumulative -> monotonic


def test_psi_of_z_closed_forms():
    assert numpy.isclose(
        madau_dickinson_psi_of_z(0.0, **MD_PARAMS), 1.0
    )  # normalised to 1 at z=0
    assert numpy.isclose(power_law_psi_of_z(0.0, lamb=3.0), 1.0)
    assert numpy.isclose(power_law_psi_of_z(1.0, lamb=3.0), 2.0**3)


def test_madau_dickinson_psi_matches_gwpopulation():
    gwpop_redshift = pytest.importorskip("gwpopulation.models.redshift")
    z = numpy.linspace(1e-6, 20, 2500)
    new = madau_dickinson_psi_of_z(z, **MD_PARAMS)
    gw = gwpop_redshift.MadauDickinsonRedshift(z_max=20).psi_of_z(
        redshift=z, **MD_PARAMS
    )
    assert numpy.allclose(new, gw, rtol=1e-12, atol=0)


def test_power_law_psi_matches_gwpopulation():
    gwpop_redshift = pytest.importorskip("gwpopulation.models.redshift")
    z = numpy.linspace(1e-6, 20, 2500)
    new = power_law_psi_of_z(z, lamb=3.0)
    gw = gwpop_redshift.PowerLawRedshift(z_max=20).psi_of_z(redshift=z, lamb=3.0)
    assert numpy.allclose(new, gw, rtol=1e-12, atol=0)


# ---------------------------------------------------------------------------
# Extrinsic parameters
# ---------------------------------------------------------------------------


def test_schutz_inclination_uses_sin_jacobian():
    # With the sin(theta) Jacobian, samples concentrate away from the poles
    # (theta = 0, pi). Without it, the poles would be over-weighted.
    from GWForge.population.extrinsic import Extrinsic

    e = Extrinsic(number_of_samples=8000, inclination_distribution="schutz")
    samples = e.sample()
    theta = numpy.asarray(samples["theta_jn"])
    assert (theta > 0).all() and (theta < numpy.pi).all()
    equatorial_fraction = numpy.mean(
        (theta > numpy.pi / 4) & (theta < 3 * numpy.pi / 4)
    )
    assert equatorial_fraction > 0.4
