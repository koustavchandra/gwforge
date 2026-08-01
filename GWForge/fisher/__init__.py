"""Fisher-matrix and DALI forecasting for GWForge.

The Fisher matrix is the Vallisneri (`arXiv:gr-qc/0703086
<https://arxiv.org/abs/gr-qc/0703086>`_) linearised-signal approximation,

.. math::

   \\Gamma_{ab} = \\langle \\partial_a h \\mid \\partial_b h \\rangle,
   \\qquad \\Sigma = \\Gamma^{-1},

summed over the detectors of a network. ``order = 2`` adds the doublet-DALI
correction of Sellentin, Quartin & Amendola (`arXiv:1401.6892
<https://arxiv.org/abs/1401.6892>`_), as applied to gravitational waves by
GWDALI (`arXiv:2307.10154 <https://arxiv.org/abs/2307.10154>`_).

Unlike GWFast, which differentiates its own JAX-native waveforms, GWForge works
with any waveform :class:`bilby.gw.WaveformGenerator` can produce. The waveform
is therefore opaque and most derivatives are numerical -- except for the five
parameters the waveform never sees, which are exact. See
:mod:`GWForge.fisher.derivatives`.
"""

from .antenna_derivatives import AntennaDerivatives
from .derivatives import WaveformDerivatives, analytic_derivatives
from .inner_product import inner_product, inner_product_matrix, optimal_snr_squared
from .matrix import FisherMatrix, check_order, covariance
from .parameters import (
    ANALYTIC_PARAMETERS,
    DEFAULT_PARAMETERS,
    ECCENTRIC_PRECESSING_PARAMETERS,
    PRECESSING_PARAMETERS,
)
from .stepsize import estimate_noise

__all__ = [
    "ANALYTIC_PARAMETERS",
    "DEFAULT_PARAMETERS",
    "ECCENTRIC_PRECESSING_PARAMETERS",
    "PRECESSING_PARAMETERS",
    "AntennaDerivatives",
    "FisherMatrix",
    "WaveformDerivatives",
    "analytic_derivatives",
    "check_order",
    "covariance",
    "estimate_noise",
    "inner_product",
    "inner_product_matrix",
    "optimal_snr_squared",
]
