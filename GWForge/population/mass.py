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

# BGP fiducials: the posterior medians of the O4b ``Default`` analysis released
# with `arXiv:2605.27226 <https://arxiv.org/abs/2605.27226>`_, from
# ``o4b_default_mass_TwoPeakBrokenPowerLawSmoothedMassDistribution_redshift_``
# ``PowerLawRedshift_magnitude_iid_spin_magnitude_gaussian_tilt_iid_spin_``
# ``orientation_result.hdf5``. The release spells several of them differently:
# ``mlow_1`` is ``mmin``, ``delta_m_1`` is ``delta_m``, ``mlow_2`` is
# ``mmin_2``, ``mmax`` is ``m_high`` and ``break_mass`` is ``m_break``.
BGP_PARAMETERS = {
    "alpha_1": 1.456442737,
    "alpha_2": 5.100400428,
    "m_break": 37.9806382,
    "mmin": 4.489078237,
    "delta_m": 3.123416382,
    # The secondary gets its own Planck taper, independent of the primary's.
    "mmin_2": 3.487749145,
    "delta_m_2": 5.60056759,
    # A delta-function prior in the analysis, not a measurement.
    "m_high": 300.0,
    "lam_0": 0.4206178692,
    "lam_1": 0.5228817702,
    # Peak 1 is the *narrow low-mass* one and carries the larger weight; peak 2
    # is the broad ~33 Msun feature. The pre-O4b defaults had them the other way
    # round, which is worth knowing when reading an old config.
    "mpp_1": 9.989384122,
    "sigpp_1": 0.6601839103,
    "mpp_2": 33.2656028,
    "sigpp_2": 4.582221363,
    "beta": 0.8049660633,
    # Upper edge of the evaluation grid, not a model parameter.
    "maximum_mass": 300.0,
}

CONDITIONAL_SAMPLING_CHUNK = 4096

NODES_PER_PEAK_WIDTH = 40

# Bounds on the primary-mass grid size. The floor is the historical default;
# the ceiling keeps the (m1, q) mesh in ``_norm_p_q`` to a few hundred MB.
MINIMUM_MASS_GRID_NODES = 1000
MAXIMUM_MASS_GRID_NODES = 40000


def primary_mass_grid_nodes(parameters, minimum_mass, maximum_mass):
    """
    Number of primary-mass grid nodes needed to resolve the narrowest peak.

    The smoothed mass models are sampled by wrapping the density in a
    :class:`bilby.core.prior.Interped`, which builds a **linear interpolant** on
    this grid and inverse-CDF samples from that. The interpolant, not the
    density, is what the catalogue is drawn from, so the grid has to resolve
    every feature the density has.

    Parameters
    ----------
    parameters : dict
        Mass-model parameters. Any key beginning with ``sigpp`` is treated as a
        peak width; models without peaks keep the floor.
    minimum_mass, maximum_mass : float
        Range the grid spans.

    Returns
    -------
    int
        Clipped to ``[MINIMUM_MASS_GRID_NODES, MAXIMUM_MASS_GRID_NODES]``.
    """
    widths = [
        float(value)
        for name, value in parameters.items()
        if name.startswith("sigpp") and float(value) > 0
    ]
    if not widths:
        return MINIMUM_MASS_GRID_NODES
    needed = NODES_PER_PEAK_WIDTH * (maximum_mass - minimum_mass) / min(widths)
    nodes = int(numpy.clip(needed, MINIMUM_MASS_GRID_NODES, MAXIMUM_MASS_GRID_NODES))
    if needed > MAXIMUM_MASS_GRID_NODES:
        logging.warning(
            "Resolving a peak of width %.3g Msun over [%g, %g] would need %d "
            "grid nodes; capping at %d. The sampled distribution will be a "
            "slightly smoothed version of the requested one.",
            min(widths),
            minimum_mass,
            maximum_mass,
            int(needed),
            MAXIMUM_MASS_GRID_NODES,
        )
    return nodes


def sample_mass_ratio_given_primary(
    model, primary_samples, mass_ratio_parameters, random_state=None
):
    r"""Draw ``mass_ratio`` conditionally on ``mass_1``, by inverse CDF.

    Parameters
    ----------
    model : BaseSmoothedMassDistribution
        Supplies the ``m1s``/``qs`` grids and ``p_q``.
    primary_samples : array_like
        Sampled ``mass_1_source``.
    mass_ratio_parameters : dict
        ``beta``, ``mmin`` and ``delta_m``, and optionally ``mmin_2`` and
        ``delta_m_2`` for models whose secondary has its own taper.
    random_state : numpy.random.Generator or None
        Defaults to numpy's global RNG, so a single seed governs it.

    Returns
    -------
    numpy.ndarray
        Mass-ratio samples, the same length as ``primary_samples``.
    """
    primary_samples = numpy.asarray(primary_samples, dtype=float)
    m1s = model.m1s
    # The secondary's own taper edge, which is where p(q | m_1) turns on. 
    minimum_mass = mass_ratio_parameters.get(
        "mmin_2", mass_ratio_parameters["mmin"]
    )

    fraction = numpy.linspace(0.0, 1.0, len(model.qs))
    lowest_ratio = numpy.clip(minimum_mass / m1s, 0.0, 1.0)[:, numpy.newaxis]
    ratio_grid = lowest_ratio + (1.0 - lowest_ratio) * fraction[numpy.newaxis, :]
    primary_grid = numpy.broadcast_to(
        m1s[:, numpy.newaxis], ratio_grid.shape
    )
    conditional = model.p_q(
        {"mass_1": primary_grid, "mass_ratio": ratio_grid}, **mass_ratio_parameters
    )

    # Row-wise cumulative integral: the grids differ per row, so this is a
    # trapezoid over each row's own spacing rather than a shared one.
    widths = numpy.diff(ratio_grid, axis=1)
    increments = 0.5 * (conditional[:, 1:] + conditional[:, :-1]) * widths
    cumulative = numpy.concatenate(
        [numpy.zeros((len(m1s), 1)), numpy.cumsum(increments, axis=1)], axis=1
    )
    totals = cumulative[:, -1:]
    # Rows where the conditional has no support at all (m1 below mmin, so no q
    # gives a valid secondary) would divide by zero; they are never selected,
    # because p_m1 vanishes there too.
    cumulative = numpy.divide(
        cumulative, totals, out=numpy.zeros_like(cumulative), where=totals > 0
    )

    uniform = (
        numpy.random.uniform(size=primary_samples.size)
        if random_state is None
        else random_state.uniform(size=primary_samples.size)
    )
    upper = numpy.clip(numpy.searchsorted(m1s, primary_samples), 1, len(m1s) - 1)
    lower = upper - 1
    weight = (primary_samples - m1s[lower]) / (m1s[upper] - m1s[lower])

    samples = numpy.empty(primary_samples.size)
    for start in range(0, primary_samples.size, CONDITIONAL_SAMPLING_CHUNK):
        stop = start + CONDITIONAL_SAMPLING_CHUNK
        piece = slice(start, stop)
        draws = uniform[piece]
        below = _invert_rows(
            cumulative[lower[piece]], ratio_grid[lower[piece]], draws
        )
        above = _invert_rows(
            cumulative[upper[piece]], ratio_grid[upper[piece]], draws
        )
        samples[piece] = (1.0 - weight[piece]) * below + weight[piece] * above
    # The bracketing rows have slightly different support edges
    return numpy.clip(samples, minimum_mass / primary_samples, 1.0)


def _invert_rows(cumulative, support, draws):
    """Invert one CDF per row at one draw per row, by linear interpolation.

    Parameters
    ----------
    cumulative : numpy.ndarray
        ``(n_rows, n_nodes)`` normalised CDFs.
    support : numpy.ndarray
        ``(n_rows, n_nodes)`` nodes -- one grid per row, because each primary
        mass has its own mass-ratio support edge.
    draws : numpy.ndarray
        ``(n_rows,)`` uniform draws.

    Returns
    -------
    numpy.ndarray
    """
    index = numpy.clip(
        numpy.sum(cumulative < draws[:, numpy.newaxis], axis=1),
        1,
        cumulative.shape[1] - 1,
    )
    rows = numpy.arange(len(draws))
    low = cumulative[rows, index - 1]
    high = cumulative[rows, index]
    span = high - low
    fraction = numpy.divide(
        draws - low, span, out=numpy.zeros_like(draws), where=span > 0
    )
    left = support[rows, index - 1]
    return left + fraction * (support[rows, index] - left)


class Mass:
    def __init__(
        self,
        mass_model,
        number_of_samples,
        parameters=None,
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
            A dictionary of model parameters. Omit it to get the model's
            fiducial values from :data:`BGP_PARAMETERS` -- the GWTC-5.0
            medians for BGP if mass_model is "bgp".
        full_pop_sampler : str
            Sampler to be used for full pop gwtc-4 model. [Options: {}] [Default: importance_m1_m2]
        """.format(choices, sampler_choices)
        self.mass_model = utils.remove_special_characters(mass_model.lower())
        self.number_of_samples = number_of_samples
        if parameters is None:
            parameters = BGP_PARAMETERS if self.mass_model == "bgp" else {}
        # Copy so the module-level defaults above cannot be mutated through an
        # instance, and so a caller's dictionary is left alone.
        self.parameters = dict(parameters)
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
            # BGP's Gaussian peaks are only left-truncated, so it needs a wider
            # evaluation grid than the 100 Msun default the others use.
            maximum_mass = (
                self.parameters.get("maximum_mass", 200)
                if "bgp" in self.mass_model
                else 100
            )
            minimum_mass = self.parameters.get("mmin", 2)
            # Enough nodes to resolve the narrowest peak, for the same reason.
            nodes = primary_mass_grid_nodes(
                self.parameters, minimum_mass, maximum_mass
            )
            logging.info(
                "Primary-mass grid: %d nodes over [%g, %g] Msun",
                nodes,
                minimum_mass,
                maximum_mass,
            )
            if "powerlawpeak" in self.mass_model:
                from ._smoothed_mass import SinglePeakSmoothedMassDistribution

                model = SinglePeakSmoothedMassDistribution(
                    mmin=minimum_mass, normalization_shape=(nodes, 1000)
                )
            elif "multipeak" in self.mass_model:
                from ._smoothed_mass import MultiPeakSmoothedMassDistribution

                model = MultiPeakSmoothedMassDistribution(
                    mmin=minimum_mass, normalization_shape=(nodes, 1000)
                )
            elif "brokenpowerlaw" in self.mass_model:
                from ._smoothed_mass import BrokenPowerLawSmoothedMassDistribution

                model = BrokenPowerLawSmoothedMassDistribution(
                    mmin=minimum_mass, normalization_shape=(nodes, 1000)
                )
            elif "bgp" in self.mass_model:
                from ._smoothed_mass import (
                    BrokenPowerLawTwoPeakSmoothedMassDistribution,
                )

                model = BrokenPowerLawTwoPeakSmoothedMassDistribution(
                    mmin=minimum_mass,
                    mmax=maximum_mass,
                    normalization_shape=(nodes, 1000),
                )
                # These describe the mass ratio, not the primary mass, and
                # maximum_mass only sizes the evaluation grid.
                excluded_from_primary = (
                    "beta",
                    "maximum_mass",
                    "mmin_2",
                    "delta_m_2",
                )

            mass1 = model.m1s

            # Create dictionaries for supported parameters
            mass_parameters = {
                param: self.parameters[param]
                for param in self.parameters
                if param not in excluded_from_primary
            }
            mass_ratio_parameters = {
                param: self.parameters[param]
                for param in self.parameters
                if param in ("beta", "mmin", "delta_m", "mmin_2", "delta_m_2")
            }

            prob_mass_1 = model.p_m1({"mass_1": mass1}, **mass_parameters)

            primary_mass_prior = bilby.core.prior.Interped(
                mass1,
                prob_mass_1,
                minimum=numpy.min(mass1),
                maximum=numpy.max(mass1),
                name="mass_1_source",
            )
            samples["mass_1_source"] = primary_mass_prior.sample(
                self.number_of_samples
            )
            # p(q | m_1), not a marginal p(q): the conditional's support is
            # q >= mmin / m_1, so drawing q independently of m_1 puts secondaries
            # below mmin, outside the model's own support. See
            # sample_mass_ratio_given_primary.
            samples["mass_ratio"] = sample_mass_ratio_given_primary(
                model, samples["mass_1_source"], mass_ratio_parameters
            )

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
