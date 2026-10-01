"""Cross-validate the population Fisher against the seminumeric implementation.

A fast version of ``validation/population_fisher_vs_seminumeric.py``: the
mass + redshift forecast at one SNR threshold, which is the comparison that runs
without the cosmology sector's finite-difference noise and can therefore be held
to a tight tolerance.

The reference is an absolute path outside the repository rather than an
installable package, so this skips when it is not there. Point at it with
``GWFORGE_SEMINUMERIC_DIRECTORY``.
"""

import os
import sys

import pytest

VALIDATION = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "validation"
)
sys.path.insert(0, VALIDATION)

SNR_THRESHOLD = 10.0


@pytest.fixture(scope="module")
def comparison():
    """The validation script's own machinery, or a skip."""
    import population_fisher_vs_seminumeric as script

    directory = script.DEFAULT_REFERENCE_DIRECTORY
    if not os.path.isdir(directory):
        pytest.skip("no seminumeric reference at {}".format(directory))
    if not os.path.exists(os.path.join(directory, script.REFERENCE_CATALOGUE)):
        pytest.skip("seminumeric reference has no catalogue checked out")
    load, forecast_mass_redshift, spectral_sirens_forecast = script.reference_modules(
        directory
    )
    fiducial = script.reference_fiducial(
        directory, forecast_mass_redshift, spectral_sirens_forecast
    )
    return script, directory, load, fiducial


def test_mass_and_redshift_reproduce_the_published_forecast(comparison):
    """sigma, parameter by parameter, against the checked-in .npz.

    Both codes run on the *same* detected events, so this is not a statistical
    agreement -- any difference is a difference in the estimator or the density.
    """
    script, directory, load, fiducial = comparison
    assert script.compare_mass_and_redshift(directory, load, fiducial, SNR_THRESHOLD)
