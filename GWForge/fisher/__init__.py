"""
Fisher-matrix and DALI forecasting for GWForge.
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
    strip_shadowed_parameters,
)
from .stepsize import estimate_noise

__all__ = [
    "ANALYTIC_PARAMETERS",
    "DEFAULT_PARAMETERS",
    "ECCENTRIC_PRECESSING_PARAMETERS",
    "PRECESSING_PARAMETERS",
    "strip_shadowed_parameters",
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
