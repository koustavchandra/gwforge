"""Unit tests for the population models that had crash/silent-physics bugs."""

import numpy

from GWForge.population.spin import Spin
from GWForge.population.mass import Mass
from GWForge.population.extrinsic import Extrinsic


def test_spin_default_model_runs():
    # "default" previously mapped to a string that never matched a valid choice,
    # so it always raised ValueError. It should now sample the isotropic
    # beta+gaussian+uniform model.
    s = Spin(
        spin_model="default",
        number_of_samples=100,
        parameters={
            "xi_spin": 0.5,
            "sigma_t": 1.0,
            "mu_chi": 0.2,
            "sigma_squared_chi": 0.02,
        },
    )
    samples = s.sample()
    assert len(samples["tilt_1"]) == 100
    assert len(samples["a_1"]) == 100


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
    from bilby.core.utils import random as bilby_random

    def run():
        numpy.random.seed(1)
        bilby_random.seed(1)
        return Spin(
            spin_model="default",
            number_of_samples=200,
            parameters={
                "xi_spin": 0.5,
                "sigma_t": 1.0,
                "mu_chi": 0.2,
                "sigma_squared_chi": 0.02,
            },
        ).sample()

    a = run()
    b = run()
    assert numpy.array_equal(a["a_1"], b["a_1"])
    assert numpy.array_equal(a["tilt_1"], b["tilt_1"])


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


def test_schutz_inclination_uses_sin_jacobian():
    # With the sin(theta) Jacobian, samples concentrate away from the poles
    # (theta = 0, pi). Without it, the poles would be over-weighted.
    e = Extrinsic(number_of_samples=8000, inclination_distribution="schutz")
    samples = e.sample()
    theta = numpy.asarray(samples["theta_jn"])
    assert (theta > 0).all() and (theta < numpy.pi).all()
    equatorial_fraction = numpy.mean(
        (theta > numpy.pi / 4) & (theta < 3 * numpy.pi / 4)
    )
    assert equatorial_fraction > 0.4
