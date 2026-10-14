r"""Population (hyper-parameter) Fisher forecasting for GWForge."""

from .base import JointPopulationModel, PopulationModel
from .catalogue import Catalogue, load_catalogue, load_injections, network_snr_from
from .fit import compare_to, fit_model
from .model import (
    BrokenPowerLawTwoPeakMass,
    DefaultSpin,
    MadauDickinsonRedshift,
    SpectralSirenModel,
)
from .population_fisher_term_I import (
    PopulationFisherResult,
    covariance_from_fisher,
    population_fisher,
)

__all__ = [
    "BrokenPowerLawTwoPeakMass",
    "Catalogue",
    "DefaultSpin",
    "JointPopulationModel",
    "MadauDickinsonRedshift",
    "PopulationFisherResult",
    "PopulationModel",
    "SpectralSirenModel",
    "compare_to",
    "covariance_from_fisher",
    "fit_model",
    "load_catalogue",
    "load_injections",
    "network_snr_from",
    "population_fisher",
]
