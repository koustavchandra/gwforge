r"""Population (hyper-parameter) Fisher forecasting for GWForge.

Where :mod:`GWForge.fisher` asks how well a network measures *one source*, this
asks how well it measures the *population*: the shape of the mass function, the
history of the merger rate, the spin distribution, and -- through spectral
sirens -- the cosmology.

The estimator is the first term of the hyper-parameter Fisher expansion of Gair,
Ghosh et al. (`arXiv:2205.07893 <https://arxiv.org/abs/2205.07893>`_ Eq. 21),

.. math::

   \Gamma^{ij} = \sum_{k\ \rm detected}
     \big(s_k^i - \langle s\rangle^i\big)\big(s_k^j - \langle s\rangle^j\big),
   \qquad
   s_k^i = \frac{\partial \ln p(\theta_k \mid \Lambda)}{\partial \Lambda^i},

which is the small-measurement-error limit. The mean subtraction is the
selection term, so there is no separate detection integral to evaluate and no
rate normalisation to get right -- see :mod:`GWForge.population_fisher.population_fisher_term_I`.

Pipeline
--------

.. code-block:: bash

   gwforge_population        --config-file bgp-gwtc5.ini --output-file population.h5
   gwforge_optimal_snr       --injection-file population.h5 --output-file snr.h5 --ifos CE40 CE20 ET
   gwforge_population_fisher --config-file mass_redshift.ini

Models
------

All four live in :mod:`GWForge.population_fisher.model`. The densities are
GWForge's own, called directly rather than reimplemented:
:mod:`GWForge.population._smoothed_mass` for the BGP mass model,
:func:`GWForge.population.redshift.madau_dickinson_psi_of_z` for the merger
rate, and :func:`GWForge.population.spin.truncated_normal` for both halves of
the ``Default`` spin model -- so a forecast and the catalogue it runs on
describe the same population by construction.

Derivatives
-----------

Analytic, everywhere, cosmology included -- see
:mod:`GWForge.population_fisher.base` for the two-tier policy and
:mod:`GWForge.cosmology` for why differentiating under the integral beats
finite-differencing the integral. Finite differences are the test oracle and the
fallback for a model added later, not a production path; there is no autodiff
tier and no JAX dependency.
"""

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
