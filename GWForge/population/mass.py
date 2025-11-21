import numpy
import logging
import bilby
from .. import utils
from .. import conversion

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)


def notch_filter(val, parameters):
    return 1.0 - parameters["A"] / (
        (1 + (parameters["gamma_low"] / val) ** parameters["eta_low"])
        * (1 + (val / parameters["gamma_high"]) ** parameters["eta_high"])
    )


def low_pass_filter(val, parameters):
    return 1.0 / (1 + (val / parameters["mmax"]) ** parameters["n"])


choices = [
    "PowerLaw+Peak",
    "MultiPeak",
    "BrokenPowerLaw",
    "UniformSecondary",
    "DoubleGaussian",
    "LogNormal",
    "PowerLawDipBreak",
    "PowerLaw",
    "Uniform_components",
]


choices = ['PowerLaw+Peak', 'MultiPeak', 'BrokenPowerLaw', 'UniformSecondary', 
           'DoubleGaussian', 'LogNormal', 'PowerLawDipBreak', 'PowerLaw', 'Uniform_components', 'Uniform_M_q']
class Mass:
    def __init__(
        self,
        mass_model,
        number_of_samples,
        parameters={
            "alpha": 3.37,
            "beta": 0.76,
            "delta_m": 5.23,
            "mmin": 4.89,
            "mmax": 88.81,
            "lam": 0.04,
            "mpp": 33.60,
            "sigpp": 4.59,
        },
    ):
        """
        Parameters:
        ----------
        mass_model : str
            The parameterized mass model. [Options: {}]
        number_of_samples : (int)
            The number of samples to generate. [Ideal: Exactly same as redshift samples]
        parameters: (dict, optional)
            A dictionary of model parameters. Default is provided assuming PowerLawPeak
        """.format(
            choices
        )
        self.mass_model = utils.remove_special_characters(mass_model.lower())
        self.number_of_samples = number_of_samples
        self.parameters = parameters

    def sample(self):
        """
        Generate mass distribution samples based on the chosen parameterised model and its parameters.

        Returns:
        --------
            dict: A dictionary containing source frame mass distribution samples.
        """
        samples = {}
        if self.mass_model in [
            "powerlawpeak",
            "multipeak",
            "brokenpowerlaw",
        ]:  # Implemented GWPopulation Models
            if "powerlawpeak" in self.mass_model:
                from gwpopulation.models.mass import SinglePeakSmoothedMassDistribution

                model = SinglePeakSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
            elif "multipeak" in self.mass_model:
                from gwpopulation.models.mass import MultiPeakSmoothedMassDistribution

                model = MultiPeakSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
            elif "brokenpowerlaw" in self.mass_model:
                from gwpopulation.models.mass import (
                    BrokenPowerLawSmoothedMassDistribution,
                )

                model = BrokenPowerLawSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )

            mass1, mass_ratio = model.m1s, model.qs

            # Create dictionaries for supported parameters
            mass_parameters = {
                param: self.parameters[param]
                for param in self.parameters
                if param not in ("beta")
            }
            mass_ratio_parameters = {
                param: self.parameters[param]
                for param in self.parameters
                if param in ("beta", "mmin", "delta_m")
            }

            prob_mass_1 = model.p_m1({"mass_1": mass1}, **mass_parameters)
            prob_mass_ratio = model.p_q(
                {"mass_ratio": mass_ratio, "mass_1": mass1}, **mass_ratio_parameters
            )

            primary_mass_prior = bilby.core.prior.Interped(
                mass1,
                prob_mass_1,
                minimum=numpy.min(mass1),
                maximum=numpy.max(mass1),
                name="mass_1_source",
            )

            mass_ratio_prior = bilby.core.prior.Interped(
                mass_ratio,
                prob_mass_ratio,
                minimum=numpy.min(mass_ratio),
                maximum=numpy.max(mass_ratio),
                name="mass_ratio",
            )
            mass_prior = bilby.gw.prior.BBHPriorDict(
                dictionary=utils.reference_prior_dict
            )
            mass_prior["mass_ratio"] = mass_ratio_prior
            mass_prior["mass_1_source"] = primary_mass_prior
            prior_samples = mass_prior.sample(self.number_of_samples)
            samples["mass_1_source"] = prior_samples["mass_1_source"]
            samples["mass_ratio"] = prior_samples["mass_ratio"]

        else:
            logging.warning("Parameterised mass model does not exist in gwpopulation")
            logging.info("Generating samples using {} model".format(self.mass_model))

            mass_prior = bilby.gw.prior.BBHPriorDict(
                dictionary=utils.reference_prior_dict
            )
            if "uniformsecondary" in self.mass_model:
                from gwpopulation.models.mass import SinglePeakSmoothedMassDistribution

                model = SinglePeakSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
                mass_parameters = {
                    param: self.parameters[param]
                    for param in self.parameters
                    if param
                    not in (
                        "beta",
                        "minimum_secondary_mass",
                        "maximum_secondary_mass",
                        "minimum_mass_ratio",
                    )
                }
                mass1 = model.m1s
                prob_mass_1 = model.p_m1({"mass_1": mass1}, **mass_parameters)
                primary_mass_prior = bilby.core.prior.Interped(
                    mass1,
                    prob_mass_1,
                    minimum=numpy.min(mass1),
                    maximum=numpy.max(mass1),
                    name="mass_1_source",
                )
                secondary_mass_prior = bilby.core.prior.analytical.Uniform(
                    minimum=self.parameters["minimum_secondary_mass"],
                    maximum=self.parameters["maximum_secondary_mass"],
                    name="mass_2_source",
                )
                # Waveform limitations + Population synthesis limitations: https://arxiv.org/abs/2009.06655
                minimum_mass_ratio = self.parameters.get("minimum_mass_ratio", 0.02)
                mass_prior["mass_ratio"] = bilby.gw.prior.Constraint(
                    minimum=minimum_mass_ratio, maximum=1, name="mass_ratio"
                )
                mass_prior["mass_1_source"] = primary_mass_prior
                mass_prior["mass_2_source"] = secondary_mass_prior

            elif "doublegaussian" in self.mass_model:
                """
                Consider checking https://arxiv.org/pdf/2005.00032.pdf
                """
                mass = numpy.linspace(
                    self.parameters["mmin"], self.parameters["mmax"], 5001
                )
                prior_1 = bilby.core.prior.analytical.TruncatedGaussian(
                    mu=self.parameters["mu_1"],
                    sigma=self.parameters["sigma_1"],
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["mmax"],
                )
                prob_1 = prior_1.prob(mass) * self.parameters["breaking_fraction"]
                prior_2 = bilby.core.prior.analytical.TruncatedGaussian(
                    mu=self.parameters["mu_2"],
                    sigma=self.parameters["sigma_2"],
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["mmax"],
                )
                prob_2 = prior_2.prob(mass) * (1 - self.parameters["breaking_fraction"])
                prob = prob_1 + prob_2
                prior = bilby.core.prior.Interped(
                    mass, prob, minimum=numpy.min(mass), maximum=numpy.max(mass)
                )
                mass_prior["mass_1_source"] = prior
                mass_prior["mass_2_source"] = prior
                mass_prior["mass_ratio"] = bilby.gw.prior.Constraint(
                    minimum=0.5, maximum=1, name="mass_ratio"
                )

            elif "lognormal" in self.mass_model or "loggaussian" in self.mass_model:
                prior = bilby.core.prior.analytical.LogNormal(
                    mu=self.parameters["mu"], sigma=self.parameters["sigma"]
                )
                mass_prior["mass_1_source"] = prior
                mass_prior["mass_2_source"] = prior
                # Chosen based on waveform limitation
                mass_prior["mass_ratio"] = bilby.gw.prior.Constraint(
                    minimum=self.parameters.get("minimum_mass_ratio", 0.056),
                    maximum=1,
                    name="mass_ratio",
                )

            elif "dip" in self.mass_model:
                """
                Consider checking Eq.1 of https://arxiv.org/pdf/2111.03498.pdf
                """
                mass = numpy.linspace(
                    self.parameters["mmin"], self.parameters["mmax"], 5001
                )
                prob = numpy.zeros_like(mass)
                prior_1 = bilby.core.prior.analytical.PowerLaw(
                    alpha=self.parameters["alpha_1"],
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["gamma_high"],
                )
                prob_1 = prior_1.prob(mass[mass <= self.parameters["gamma_high"]])
                prob[mass <= self.parameters["gamma_high"]] = (
                    prob_1
                    * notch_filter(
                        val=mass[mass <= self.parameters["gamma_high"]],
                        parameters=self.parameters,
                    )
                    * low_pass_filter(
                        val=mass[mass <= self.parameters["gamma_high"]],
                        parameters=self.parameters,
                    )
                )

                prior_2 = bilby.core.prior.analytical.PowerLaw(
                    alpha=self.parameters["alpha_2"],
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["gamma_high"],
                )
                prob_2 = prior_2.prob(mass[mass > self.parameters["gamma_high"]])
                prob[mass > self.parameters["gamma_high"]] = (
                    prob_2
                    * notch_filter(
                        val=mass[mass > self.parameters["gamma_high"]],
                        parameters=self.parameters,
                    )
                    * low_pass_filter(
                        val=mass[mass > self.parameters["gamma_high"]],
                        parameters=self.parameters,
                    )
                )
                prior = bilby.core.prior.Interped(
                    mass, prob, minimum=numpy.min(mass), maximum=numpy.max(mass)
                )
                mass_prior["mass_1_source"] = prior
                mass_prior["mass_2_source"] = prior

            elif "powerlaw" in self.mass_model:
                prior = bilby.core.prior.analytical.PowerLaw(
                    alpha=self.parameters["alpha"],
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["mmax"],
                )

                minimum_mass_ratio = self.parameters.get("minimum_mass_ratio", 0.056)
                mass_prior["mass_ratio"] = bilby.gw.prior.Constraint(
                    minimum=minimum_mass_ratio, maximum=1, name="mass_ratio"
                )
                mass_prior["mass_1_source"] = prior
                mass_prior["mass_2_source"] = prior

                prior_samples = mass_prior.sample(self.number_of_samples)
                samples['mass_1_source'] = prior_samples['mass_1_source']
                samples['mass_2_source'] = prior_samples['mass_2_source']
                
            elif 'fixed' in self.mass_model:
                samples['mass_1_source'] = numpy.ones(self.number_of_samples) * self.parameters['primary_mass']
                samples['mass_2_source'] = samples['mass_1_source'] * (self.parameters['mass_ratio'] if self.parameters['mass_ratio'] < 1 else 1 / self.parameters['mass_ratio'])
            elif self.mass_model == 'uniformcomponents':
                mass_prior['mass_1_source'] = bilby.core.prior.analytical.Uniform(minimum=self.parameters['mmin'], maximum=self.parameters['mmax'], name='mass_1_source')
                mass_prior['mass_2_source'] = bilby.core.prior.analytical.Uniform(minimum=self.parameters['mmin'], maximum=self.parameters['mmax'], name='mass_2_source')
            elif self.mass_model == 'uniformmq':
                mass_prior['total_mass_source'] = bilby.core.prior.analytical.Uniform(minimum=self.parameters['minimum_total_mass'], maximum=self.parameters['maximum_total_mass'], name='total_mass_source')
                mass_prior['mass_ratio'] = bilby.core.prior.analytical.Uniform(minimum=self.parameters['minimum_mass_ratio'], maximum=self.parameters['maximum_mass_ratio'], name='mass_ratio')
            else:
                raise ValueError(
                    "{} is not implemented in gwpopulation. Please choose from {}".format(
                        self.mass_model, choices
                    )
                )

            prior_samples = mass_prior.sample(self.number_of_samples)
            samples["mass_1_source"] = prior_samples["mass_1_source"]
            samples["mass_2_source"] = prior_samples["mass_2_source"]

            if self.mass_model == "uniformcomponents":
                m1_tmp = samples["mass_1_source"]
                m2_tmp = samples["mass_2_source"]
                samples["mass_1_source"] = numpy.where(m1_tmp > m2_tmp, m1_tmp, m2_tmp)
                samples["mass_2_source"] = numpy.where(m1_tmp > m2_tmp, m2_tmp, m1_tmp)
        # Generate all source frame mass parameters from samples
        samples = conversion.generate_mass_parameters(samples, source=True)
        return samples
