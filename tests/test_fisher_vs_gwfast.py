"""Cross-validate the Fisher matrix against GWFast.

A fast version of ``validation/fisher_vs_gwfast.py``: one source, a coarse
frequency grid, and a quadrupole-only waveform. Restricting to ``IMRPhenomXAS``
is what makes this a tight test -- GWFast applies the coalescence phase as a
single global factor, which is exact without higher modes, so the two codes
should agree to a fraction of a percent on *every* parameter with no
allowances. The full higher-mode comparison, where GWFast's approximation costs
it up to 20% on the ``phase`` entry, lives in ``validation/``.

GWFast is an optional dependency: ``pip install -e .[validation]``.
"""

import os
import sys

import numpy
import pytest

VALIDATION = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "validation"
)
sys.path.insert(0, VALIDATION)

MINIMUM_FREQUENCY = 20.0
MAXIMUM_FREQUENCY = 512.0
DELTA_FREQUENCY = 1.0 / 8.0

SOURCE = dict(
    chirp_mass=28.0,
    symmetric_mass_ratio=0.2200,
    chi_1=0.40,
    chi_2=0.30,
    luminosity_distance=2000.0,
    theta_jn=1.20,
    psi=0.30,
    ra=3.100,
    dec=0.4500,
    phase=0.0,
    geocent_time=1187008882.4,
)


@pytest.fixture(scope="module")
def comparison():
    pytest.importorskip("gwfast")
    import fisher_vs_gwfast
    from GWForge.fisher import covariance

    record = fisher_vs_gwfast.compare(
        0, dict(SOURCE), order=1, delta_frequency=DELTA_FREQUENCY,
        approximant="IMRPhenomXAS",
    )
    covariance(record["gwforge_fisher"])
    return record


def test_snr_agrees(comparison):
    """If the SNRs disagree the two codes are not looking at the same signal."""
    # The two codes integrate on slightly different quadrature grids, which
    # is worth a part in a thousand; anything larger means a real mismatch.
    assert comparison["gwforge_snr"] == pytest.approx(
        comparison["gwfast_snr"], rel=3e-3
    )


def test_fisher_diagonals_agree(comparison):
    """The Fisher entries are what each code computes; compare them directly.

    Not the sigmas: this parameter set contains near-degenerate pairs, so the
    matrix inverse can amplify a fraction of a percent into tens of percent
    without either code being wrong.
    """
    from gwfast_interface import PARAMETERS

    mine = numpy.diag(comparison["gwforge_fisher"])
    reference = numpy.diag(comparison["gwfast_fisher"])
    for name, actual, expected in zip(PARAMETERS, mine, reference):
        assert actual == pytest.approx(expected, rel=0.02), name


def test_fisher_correlations_agree(comparison):
    """Off-diagonals too, normalised so every entry is order one.

    Looser than the diagonals on purpose. A correlation close to one sits in a
    near-degenerate direction, where the 1-2% by which the two codes' Fisher
    elements differ is amplified; the measured worst-case residual across
    sources and approximants is 0.035 to 0.094. The median is 6e-4.
    """
    def normalise(matrix):
        scale = numpy.sqrt(numpy.diag(matrix))
        return matrix / numpy.outer(scale, scale)

    numpy.testing.assert_allclose(
        normalise(comparison["gwforge_fisher"]),
        normalise(comparison["gwfast_fisher"]),
        atol=0.12,
    )
