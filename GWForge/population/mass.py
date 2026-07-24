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
    "BGP",
    "UniformSecondary",
    "DoubleGaussian",
    "LogNormal",
    "PowerLawDipBreak",
    "PowerLaw",
    "Uniform_components",
    "Uniform_M_q",
    "FullPop_GWTC4",
    "UserDefined",
]

sampler_choices = ["importance_m1_m2", "importance_m1_q", "lint"]


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
        full_pop_sampler="importance_m1_m2",
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
        full_pop_sampler : str
            Sampler to be used for full pop gwtc-4 model. [Options: {}] [Default: importance_m1_m2]
        """.format(choices, sampler_choices)
        self.mass_model = utils.remove_special_characters(mass_model.lower())
        self.number_of_samples = number_of_samples
        self.parameters = parameters
        self.full_pop_sampler = full_pop_sampler

    def sample(self, **sampler_kwargs):
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
            "bgp",
        ]:  # Smoothed mass models (in-house, formerly from GWPopulation)
            # Parameters not passed to the primary-mass model (beta drives the mass
            # ratio; maximum_mass only sizes the BGP evaluation grid).
            excluded_from_primary = ("beta",)
            if "powerlawpeak" in self.mass_model:
                from ._smoothed_mass import SinglePeakSmoothedMassDistribution

                model = SinglePeakSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
            elif "multipeak" in self.mass_model:
                from ._smoothed_mass import MultiPeakSmoothedMassDistribution

                model = MultiPeakSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
            elif "brokenpowerlaw" in self.mass_model:
                from ._smoothed_mass import BrokenPowerLawSmoothedMassDistribution

                model = BrokenPowerLawSmoothedMassDistribution(
                    normalization_shape=(1000, 1000)
                )
            elif "bgp" in self.mass_model:
                from ._smoothed_mass import (
                    BrokenPowerLawTwoPeakSmoothedMassDistribution,
                )

                # BGP's Gaussian peaks are only left-truncated, so allow a wider
                # evaluation grid than the default 100 Msun upper bound.
                maximum_mass = self.parameters.get("maximum_mass", 200)
                model = BrokenPowerLawTwoPeakSmoothedMassDistribution(
                    mmax=maximum_mass, normalization_shape=(1000, 1000)
                )
                excluded_from_primary = ("beta", "maximum_mass")

            mass1, mass_ratio = model.m1s, model.qs

            # Create dictionaries for supported parameters
            mass_parameters = {
                param: self.parameters[param]
                for param in self.parameters
                if param not in excluded_from_primary
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

        elif self.mass_model == "fullpopgwtc4":
            logging.info("Generating samples using {} model".format(self.mass_model))
            mass_prior = bilby.gw.prior.BBHPriorDict(
                dictionary=utils.reference_prior_dict
            )

            if self.full_pop_sampler == "importance_m1_m2":
                from .pdb_mass_sampler import importance_sampling_m1_m2_prop

                m1_samples, m2_samples, ess = importance_sampling_m1_m2_prop(
                    n_samples=self.number_of_samples,
                    verbose=True,
                    **{
                        param: self.parameters[param]
                        for param in self.parameters
                        if param
                        in (
                            "A",
                            "A2",
                            "NSmin",
                            "NSmax",
                            "BHmin",
                            "BHmax",
                            "UPPERmin",
                            "UPPERmax",
                            "n0",
                            "n1",
                            "n2",
                            "n3",
                            "n4",
                            "n5",
                            "alpha_1",
                            "alpha_2",
                            "alpha_dip",
                            "mu1",
                            "sig1",
                            "mix1",
                            "mu2",
                            "sig2",
                            "mix2",
                            "beta_pair_1",
                            "beta_pair_2",
                            "mbreak",
                            "mmin",
                            "mmax",
                        )
                    },
                    **sampler_kwargs,
                )
            elif self.full_pop_sampler == "importance_m1_q":
                from .pdb_mass_sampler import importance_sampling_m1_q_prop

                m1_samples, m2_samples, ess = importance_sampling_m1_q_prop(
                    n_samples=self.number_of_samples,
                    verbose=True,
                    **{
                        param: self.parameters[param]
                        for param in self.parameters
                        if param
                        in (
                            "A",
                            "A2",
                            "NSmin",
                            "NSmax",
                            "BHmin",
                            "BHmax",
                            "UPPERmin",
                            "UPPERmax",
                            "n0",
                            "n1",
                            "n2",
                            "n3",
                            "n4",
                            "n5",
                            "alpha_1",
                            "alpha_2",
                            "alpha_dip",
                            "mu1",
                            "sig1",
                            "mix1",
                            "mu2",
                            "sig2",
                            "mix2",
                            "beta_pair_1",
                            "beta_pair_2",
                            "mbreak",
                            "mmin",
                            "mmax",
                        )
                    },
                    **sampler_kwargs,
                )
            elif self.full_pop_sampler == "lint":
                from .pdb_mass_sampler import lintsampling

                m1_samples, m2_samples = lintsampling(
                    n_samples=self.number_of_samples,
                    verbose=True,
                    **{
                        param: self.parameters[param]
                        for param in self.parameters
                        if param
                        in (
                            "A",
                            "A2",
                            "NSmin",
                            "NSmax",
                            "BHmin",
                            "BHmax",
                            "UPPERmin",
                            "UPPERmax",
                            "n0",
                            "n1",
                            "n2",
                            "n3",
                            "n4",
                            "n5",
                            "alpha_1",
                            "alpha_2",
                            "alpha_dip",
                            "mu1",
                            "sig1",
                            "mix1",
                            "mu2",
                            "sig2",
                            "mix2",
                            "beta_pair_1",
                            "beta_pair_2",
                            "mbreak",
                            "mmin",
                            "mmax",
                        )
                    },
                    **sampler_kwargs,
                )
            else:
                raise ValueError(
                    f"{self.full_pop_sampler} sampler not available. Please choose from {sampler_choices}."
                )
            samples["mass_1_source"] = m1_samples
            samples["mass_2_source"] = m2_samples
        else:
            logging.info("Using an in-house parameterised mass model")
            mass_prior = bilby.gw.prior.BBHPriorDict(
                dictionary=utils.reference_prior_dict
            )

            if "uniformsecondary" in self.mass_model:
                from ._smoothed_mass import SinglePeakSmoothedMassDistribution

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
                    minimum=self.parameters["gamma_high"],
                    maximum=self.parameters["mmax"],
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
                samples["mass_1_source"] = prior_samples["mass_1_source"]
                samples["mass_2_source"] = prior_samples["mass_2_source"]

            elif "fixed" in self.mass_model:
                samples["mass_1_source"] = (
                    numpy.ones(self.number_of_samples) * self.parameters["primary_mass"]
                )
                samples["mass_2_source"] = samples["mass_1_source"] * (
                    self.parameters["mass_ratio"]
                    if self.parameters["mass_ratio"] < 1
                    else 1 / self.parameters["mass_ratio"]
                )
            elif "userdefined" in self.mass_model:
                # Tabulated (xx, yy) distributions from an external population, e.g.
                # from a population-synthesis run. The user provides a JSON file and
                # the primary-mass key, plus EITHER a mass-ratio key OR a
                # secondary-mass key. See GWForge.population.user_defined.
                from .user_defined import interped_prior_from_file

                prior_file = self.parameters["file"]
                primary_parameter = self.parameters.get(
                    "primary_parameter", "mass_1_source"
                )
                primary_prior = interped_prior_from_file(
                    prior_file, primary_parameter, name="mass_1_source"
                )
                samples["mass_1_source"] = primary_prior.sample(self.number_of_samples)

                if "mass_ratio_parameter" in self.parameters:
                    ratio_prior = interped_prior_from_file(
                        prior_file,
                        self.parameters["mass_ratio_parameter"],
                        name="mass_ratio",
                    )
                    samples["mass_ratio"] = ratio_prior.sample(self.number_of_samples)
                elif "secondary_parameter" in self.parameters:
                    secondary_prior = interped_prior_from_file(
                        prior_file,
                        self.parameters["secondary_parameter"],
                        name="mass_2_source",
                    )
                    secondary = secondary_prior.sample(self.number_of_samples)
                    # Independent marginals do not preserve the source pairing; enforce
                    # m1 >= m2 by ordering each pair (logged so the user is aware).
                    logging.info(
                        "UserDefined: primary and secondary sampled independently; enforcing mass_1 >= mass_2 by ordering."
                    )
                    primary = samples["mass_1_source"]
                    samples["mass_1_source"] = numpy.maximum(primary, secondary)
                    samples["mass_2_source"] = numpy.minimum(primary, secondary)
                else:
                    raise ValueError(
                        "UserDefined mass model needs either 'mass_ratio_parameter' or 'secondary_parameter' in mass-parameters."
                    )
            elif self.mass_model == "uniformcomponents":
                mass_prior["mass_1_source"] = bilby.core.prior.analytical.Uniform(
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["mmax"],
                    name="mass_1_source",
                )
                mass_prior["mass_2_source"] = bilby.core.prior.analytical.Uniform(
                    minimum=self.parameters["mmin"],
                    maximum=self.parameters["mmax"],
                    name="mass_2_source",
                )
            elif self.mass_model == "uniformmq":
                mass_prior["total_mass_source"] = bilby.core.prior.analytical.Uniform(
                    minimum=self.parameters["minimum_total_mass"],
                    maximum=self.parameters["maximum_total_mass"],
                    name="total_mass_source",
                )
                mass_prior["mass_ratio"] = bilby.core.prior.analytical.Uniform(
                    minimum=self.parameters["minimum_mass_ratio"],
                    maximum=self.parameters["maximum_mass_ratio"],
                    name="mass_ratio",
                )
            else:
                raise ValueError(
                    "{} is not an implemented mass model. Please choose from {}".format(
                        self.mass_model, choices
                    )
                )

            # Only sample from mass_prior for models that set priors on it but did
            # not already populate `samples` directly (e.g. "powerlaw" and "fixed"
            # fill `samples` themselves). This avoids re-sampling/overwriting the
            # "powerlaw" draw and avoids a KeyError for "fixed" (whose mass_prior
            # has no mass_1_source key).
            already_sampled = {
                "mass_1_source",
                "mass_2_source",
                "total_mass_source",
            } & samples.keys()
            if not already_sampled:
                prior_samples = mass_prior.sample(self.number_of_samples)
                if self.mass_model == "uniformmq":
                    samples["total_mass_source"] = prior_samples["total_mass_source"]
                    samples["mass_ratio"] = prior_samples["mass_ratio"]
                else:
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
