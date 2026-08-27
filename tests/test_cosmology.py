"""Tests for :mod:`GWForge.cosmology`.

Two things are being pinned here.

The **background** is checked against ``astropy``, which is the reference
implementation everything in GWForge used before this module existed. Agreement
has to be at the 1e-8 level, not the 1e-3 level, because the whole point of the
module is that a spectral-siren forecast measures :math:`H_0` to a fraction of a
percent and cannot afford a systematic in its own distance-redshift relation.

The **derivatives** are checked against finite differences. They are analytic --
obtained by differentiating under the comoving-distance integral rather than
finite-differencing the integral -- and that is the claim these tests exist to
verify.

The last test is a regression: ``bin/gwforge_population`` used to sample
redshifts with the configured cosmology and then label them with luminosity
distances from ``bilby``'s default, Planck15. That is a 0.10-0.13% error in
exactly the observable a spectral siren measures.
"""

import numpy
import pytest
from astropy import units
from astropy.cosmology import FlatwCDM as AstropyFlatwCDM
from astropy.cosmology import LambdaCDM, Planck15, Planck18

from GWForge.cosmology import (
    COSMOLOGY_PARAMETERS,
    FlatwCDM,
    astropy_cosmology,
    differential_comoving_volume,
    luminosity_distance,
)

REDSHIFTS = numpy.array([1e-3, 0.05, 0.1, 0.5, 1.0, 3.0, 7.0, 10.0, 20.0])
EQUATIONS_OF_STATE = (-1.2, -1.0, -0.8)
# Finite-difference steps, chosen well inside the O(h^2) regime for each scale.
STEPS = {"H0": 1e-4, "Om0": 1e-6, "w0": 1e-6}


def reference(w0):
    """The astropy equivalent of ``FlatwCDM(67.66, 0.3111, w0)``."""
    return AstropyFlatwCDM(H0=67.66, Om0=0.3111, w0=w0, Tcmb0=0.0)


@pytest.mark.parametrize("w0", EQUATIONS_OF_STATE)
def test_distances_match_astropy(w0):
    """Comoving and luminosity distance against astropy's own integrator."""
    cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=w0)
    expected = reference(w0).luminosity_distance(REDSHIFTS).to(units.Mpc).value
    numpy.testing.assert_allclose(
        cosmology.luminosity_distance(REDSHIFTS), expected, rtol=1e-8
    )


@pytest.mark.parametrize("w0", EQUATIONS_OF_STATE)
def test_differential_comoving_volume_matches_astropy(w0):
    """dV_c/dz, remembering that astropy's is per steradian and ours is full sky."""
    cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=w0)
    expected = (
        reference(w0).differential_comoving_volume(REDSHIFTS)
        * 4.0
        * numpy.pi
        * units.sr
    ).to(units.Mpc**3).value
    numpy.testing.assert_allclose(
        cosmology.differential_comoving_volume(REDSHIFTS), expected, rtol=1e-8
    )


@pytest.mark.parametrize("w0", EQUATIONS_OF_STATE)
def test_ddL_dz_and_second_derivative(w0):
    """dd_L/dz and d2d_L/dz2 against differences of the closed-form level above."""
    cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=w0)
    step = 1e-5
    first = (
        cosmology.luminosity_distance(REDSHIFTS + step)
        - cosmology.luminosity_distance(REDSHIFTS - step)
    ) / (2.0 * step)
    numpy.testing.assert_allclose(cosmology.ddL_dz(REDSHIFTS), first, rtol=1e-7)

    second = (
        cosmology.ddL_dz(REDSHIFTS + step) - cosmology.ddL_dz(REDSHIFTS - step)
    ) / (2.0 * step)
    numpy.testing.assert_allclose(cosmology.d2dL_dz2(REDSHIFTS), second, rtol=1e-6)


@pytest.mark.parametrize("w0", EQUATIONS_OF_STATE)
def test_redshift_of_distance_round_trips(w0):
    """Newton inversion returns the redshift the forward model started from.

    Tight because the inversion is iterated against the forward model itself,
    not read off an interpolating table -- so there is no grid error to inherit.
    """
    cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=w0)
    distance = cosmology.luminosity_distance(REDSHIFTS)
    numpy.testing.assert_allclose(
        cosmology.redshift_of_distance(distance), REDSHIFTS, rtol=1e-10, atol=1e-12
    )


@pytest.mark.parametrize("w0", (-1.0, -0.85))
@pytest.mark.parametrize("name", COSMOLOGY_PARAMETERS)
@pytest.mark.parametrize(
    "quantity",
    ("comoving_distance", "luminosity_distance", "ddL_dz", "differential_comoving_volume"),
)
def test_derivatives_match_finite_differences(w0, name, quantity):
    """The analytic derivatives, which is the reason this module exists.

    They come from differentiating *under* the comoving-distance integral, so
    the integrand stays closed form and the accuracy is the quadrature's rather
    than a step size's. The previous generation of this code finite-differenced
    the integral instead and had an error floor that got *worse* as the step
    shrank.
    """
    settings = {"H0": 67.66, "Om0": 0.3111, "w0": w0}
    cosmology = FlatwCDM(**settings)
    analytic = cosmology.derivatives(REDSHIFTS, parameters=(name,))[name][quantity]

    step = STEPS[name]
    above = FlatwCDM(**dict(settings, **{name: settings[name] + step}))
    below = FlatwCDM(**dict(settings, **{name: settings[name] - step}))
    expected = (
        getattr(above, quantity)(REDSHIFTS) - getattr(below, quantity)(REDSHIFTS)
    ) / (2.0 * step)
    numpy.testing.assert_allclose(analytic, expected, rtol=1e-6)


@pytest.mark.parametrize("name", COSMOLOGY_PARAMETERS)
def test_dredshift_dparameter_matches_finite_differences(name):
    """dz/dLambda at *fixed distance*, the chain the spectral siren rides on."""
    settings = {"H0": 67.66, "Om0": 0.3111, "w0": -1.0}
    cosmology = FlatwCDM(**settings)
    distance = cosmology.luminosity_distance(REDSHIFTS)
    analytic = cosmology.dredshift_dparameter(REDSHIFTS, parameters=(name,))[name]

    step = STEPS[name]
    above = FlatwCDM(**dict(settings, **{name: settings[name] + step}))
    below = FlatwCDM(**dict(settings, **{name: settings[name] - step}))
    expected = (
        above.redshift_of_distance(distance) - below.redshift_of_distance(distance)
    ) / (2.0 * step)
    numpy.testing.assert_allclose(analytic, expected, rtol=1e-6)


def test_derivatives_rejects_unknown_parameter():
    with pytest.raises(ValueError, match="Unknown cosmology parameter"):
        FlatwCDM().derivatives(1.0, parameters=("Ode0",))


def test_constructor_rejects_unphysical_density():
    with pytest.raises(ValueError, match="Om0 must lie"):
        FlatwCDM(Om0=1.5)


def test_astropy_cosmology_resolves_realizations_and_custom():
    """A realization name wins; anything else needs explicit densities."""
    assert astropy_cosmology("Planck18") is Planck18
    custom = astropy_cosmology("custom", H0=70.0, Om0=0.3, Ode0=0.7, Tcmb0=0.0)
    assert isinstance(custom, LambdaCDM)
    assert custom.H0.value == pytest.approx(70.0)
    with pytest.raises(ValueError, match="not an astropy realization"):
        astropy_cosmology("not_a_realization")


def test_to_astropy_round_trips():
    """Tcmb0 = 0 and no neutrinos makes the astropy twin exact, not approximate."""
    cosmology = FlatwCDM(H0=70.0, Om0=0.28, w0=-1.0)
    twin = cosmology.to_astropy()
    numpy.testing.assert_allclose(
        cosmology.luminosity_distance(REDSHIFTS),
        twin.luminosity_distance(REDSHIFTS).to(units.Mpc).value,
        rtol=1e-8,
    )


def test_from_astropy_folds_neutrinos():
    cosmology = FlatwCDM.from_astropy(Planck18)
    assert cosmology.Om0 == pytest.approx(Planck18.Om0 + Planck18.Onu0)
    assert cosmology.H0 == pytest.approx(Planck18.H0.value)


def test_helpers_accept_either_cosmology_flavour():
    """The population code passes astropy objects; the Fisher passes FlatwCDM."""
    ours = FlatwCDM.from_astropy(Planck18)
    for function, reference_value in (
        (luminosity_distance, ours.luminosity_distance(REDSHIFTS)),
        (differential_comoving_volume, ours.differential_comoving_volume(REDSHIFTS)),
    ):
        numpy.testing.assert_allclose(
            function(ours, REDSHIFTS), reference_value, rtol=1e-12
        )
        # The astropy branch returns plain floats in Mpc, not Quantities.
        values = function(Planck18, REDSHIFTS)
        assert isinstance(values, numpy.ndarray)
        numpy.testing.assert_allclose(values, reference_value, rtol=5e-3)


def test_population_distances_use_the_configured_cosmology():
    """Regression: distances must not silently come from bilby's Planck15.

    ``bin/gwforge_population`` calls
    :func:`GWForge.cosmology.luminosity_distance` with the cosmology from the
    ini file. It used to call ``bilby.gw.conversion``'s version, which falls back
    to bilby's *default* cosmology whatever the config said. The two differ by
    0.10% at z = 0.1 and 0.13% at z = 10 -- immaterial for an injection
    campaign, and comparable to the signal for a spectral siren.
    """
    redshifts = numpy.array([0.1, 1.0, 3.0, 10.0])
    configured = luminosity_distance(astropy_cosmology("Planck18"), redshifts)
    numpy.testing.assert_allclose(
        configured, Planck18.luminosity_distance(redshifts).to(units.Mpc).value,
        rtol=1e-12,
    )
    # And the two cosmologies really are distinguishable at this precision, so
    # the assertion above has teeth.
    wrong = Planck15.luminosity_distance(redshifts).to(units.Mpc).value
    assert numpy.max(numpy.abs(configured / wrong - 1.0)) > 1e-4
