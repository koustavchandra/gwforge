import numpy
import logging
import bilby
from scipy.special import ndtr
from .. import utils
from ..conversion import *

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

choices = [
    "Non-spinning",
    "Aligned",
    "Aligned-Bilby",
    "Aligned-Uniform",
    "Beta-Aligned",
    "Aligned-Gaussian-Uniform",
    "Gaussian-Non_spinning",
    "Isotropic-Bilby",
    "Isotropic-Beta",
    "Isotropic-Beta_Gaussian",
    "Isotropic-Beta_Gaussian_Uniform",
    "Default",
]

# Default BBH spin fiducials: the posterior medians of the O4b ``Default``
# analysis released with `arXiv:2605.27226 <https://arxiv.org/abs/2605.27226>`_.
# The release calls ``mu_t`` and ``sigma_t`` ``mu_spin`` and ``sigma_spin``.
# ``sigma_chi`` is a standard deviation, not the variance the Beta models take.
#
# This is the *only* copy of these numbers; the Fisher model
# (:class:`GWForge.population_fisher.DefaultSpin`) and the shipped configs read
# from here.
DEFAULT_BBH_SPIN_PARAMETERS = {
    "mu_chi": 0.0751318233,
    "sigma_chi": 0.3667469243,
    "mu_t": 0.2788662234,
    "sigma_t": 0.9661688072,
    "xi_spin": 0.6650168221,
}

# 1 / sqrt(2 pi), for the standard normal written out below.
_INVERSE_ROOT_TWO_PI = 1.0 / numpy.sqrt(2.0 * numpy.pi)


def _standard_normal(value):
    """Standard normal density."""
    return _INVERSE_ROOT_TWO_PI * numpy.exp(-0.5 * numpy.asarray(value, float) ** 2)


def truncated_normal(value, mu, sigma, low, high):
    r"""Truncated-normal density and its log-derivatives in ``mu`` and ``sigma``.

    The kernel both halves of the ``Default`` BBH spin model are built from: the
    magnitudes are one of these on :math:`[0, a_{\max}]`, and the cos-tilt
    mixture's Gaussian component is one on :math:`[t_{\min}, 1]`.

    With :math:`z = (x-\mu)/\sigma`, :math:`a = (\text{low}-\mu)/\sigma`,
    :math:`b = (\text{high}-\mu)/\sigma` and :math:`D = \Phi(b) - \Phi(a)`,

    .. math::

       \frac{\partial \ln N}{\partial \mu}
         = \frac{z}{\sigma} - \frac{\phi(a) - \phi(b)}{\sigma D},
       \qquad
       \frac{\partial \ln N}{\partial \sigma}
         = \frac{z^2 - 1}{\sigma} - \frac{a\phi(a) - b\phi(b)}{\sigma D}.

    Parameters
    ----------
    value : array_like
    mu, sigma, low, high : float

    Returns
    -------
    dict
        ``{"density": ..., "mu": ..., "sigma": ...}``, the latter two being
        ``d ln N / d mu`` and ``d ln N / d sigma``. ``density`` is zero outside
        ``[low, high]``.
    """
    value = numpy.asarray(value, dtype=float)
    inside = (value >= low) & (value <= high)
    scaled = (value - mu) / sigma
    lower = (low - mu) / sigma
    upper = (high - mu) / sigma
    normalisation = ndtr(upper) - ndtr(lower)
    lower_pdf = _standard_normal(lower)
    upper_pdf = _standard_normal(upper)

    density = numpy.where(
        inside, _standard_normal(scaled) / (sigma * normalisation), 0.0
    )
    edge_mu = (lower_pdf - upper_pdf) / (sigma * normalisation)
    edge_sigma = (lower * lower_pdf - upper * upper_pdf) / (sigma * normalisation)
    return {
        "density": density,
        "mu": scaled / sigma - edge_mu,
        "sigma": (scaled**2 - 1.0) / sigma - edge_sigma,
    }


def default_spin_magnitude_density(magnitude, mu_chi, sigma_chi, amax=1.0):
    r"""``p(chi)`` for the ``Default`` BBH model: Eq. B15 of arXiv:2605.27226.

    A truncated Gaussian :math:`N_{[0, a_{\max}]}(\mu_\chi, \sigma_\chi)`, the
    same for both components. ``sigma_chi`` is a **standard deviation**; the Beta
    models take a variance and are a different model.
    """
    return truncated_normal(magnitude, mu_chi, sigma_chi, 0.0, amax)["density"]


def default_spin_tilt_density(cosine, mu_t, sigma_t, xi_spin, t_min=-1.0):
    r"""The **marginal** ``p(cos t)`` of one component: Eq. B16, one tilt at a time.

    .. math::
        p(\cos t) = \xi\, N_{[t_{\min}, 1]}(\mu_t, \sigma_t)
                    + \frac{1 - \xi}{1 - t_{\min}}.

    Careful: the model is joint over the *binary* -- one draw decides whether
    both tilts come from the Gaussian or both from the isotropic component -- so
    the density of the pair is **not** the product of two of these. It is

    .. math::
        p(\cos t_1, \cos t_2) = \xi\, N(\cos t_1) N(\cos t_2)
                              + \frac{1 - \xi}{(1 - t_{\min})^2},

    which is what :class:`GWForge.population_fisher.DefaultSpin` evaluates and
    what makes ``corr(cos t_1, cos t_2)`` non-zero. This marginal is the right
    thing to compare a one-dimensional histogram against, and nothing else.
    """
    gaussian = truncated_normal(cosine, mu_t, sigma_t, t_min, 1.0)["density"]
    return xi_spin * gaussian + (1.0 - xi_spin) / (1.0 - t_min)


# Fiducial ``spin-parameters`` per model, keyed by the normalised model name.
# Only the GWTC-5.0 ``Default`` BBH model has one
DEFAULT_PARAMETERS = {
    "default": DEFAULT_BBH_SPIN_PARAMETERS,
}


class Spin:
    def __init__(self, spin_model, number_of_samples, parameters=None):
        """
        Parameters:
        ----------
        spin_model : str
            The parameterized spin model. [Options: {}]
        number_of_samples : (int)
            The number of samples to generate. [Ideal: Exactly same as redshift samples]
        parameters: (dict, optional)
            A dictionary of model parameters. Omit it to get the model's
            fiducial values from :data:`DEFAULT_PARAMETERS`, which for
            ``Default`` are the GWTC-5.0 medians. [Default: Empty]
        """.format(", ".join(choices))

        self.spin_model = utils.remove_special_characters(spin_model.lower())
        self.number_of_samples = number_of_samples
        # Copy into a fresh dict so the mutations below never leak into a caller's
        # dict or a shared default argument across instances.
        self.parameters = (
            dict(parameters)
            if parameters
            else dict(DEFAULT_PARAMETERS.get(self.spin_model, {}))
        )
        self.parameters["minimum_primary_spin"] = self.parameters.get(
            "minimum_primary_spin", 0
        )
        self.parameters["maximum_primary_spin"] = self.parameters.get(
            "maximum_primary_spin", 0.99
        )
        self.parameters["minimum_secondary_spin"] = self.parameters.get(
            "minimum_secondary_spin", 0
        )
        self.parameters["maximum_secondary_spin"] = self.parameters.get(
            "maximum_secondary_spin", 0.99
        )
        if "mu_chi" in self.parameters and "sigma_squared_chi" in self.parameters:
            self.parameters["alpha_chi"] = (
                self.parameters["mu_chi"]
                * (
                    self.parameters["mu_chi"]
                    - self.parameters["mu_chi"] ** 2
                    - self.parameters["sigma_squared_chi"]
                )
                / self.parameters["sigma_squared_chi"]
            )
            self.parameters["beta_chi"] = self.parameters["alpha_chi"] * (
                1.0 / self.parameters["mu_chi"] - 1.0
            )

        self.parameters.setdefault("amax", 1.0)
        self.parameters.setdefault("t_min", -1.0)

    def _sample_default_bbh(self):
        r"""The GWTC-4.0/5.0 **Default BBH** spin model, Eqs. B15-B16 of
        `arXiv:2605.27226 <https://arxiv.org/abs/2605.27226>`_.

        Magnitudes are independent and identically distributed truncated
        Gaussians (Eq. B15),

        .. math::

           \pi(\chi_i \mid \mu_\chi, \sigma_\chi)
             = N_{[0, a_{\max}]}(\chi_1 \mid \mu_\chi, \sigma_\chi)\,
               N_{[0, a_{\max}]}(\chi_2 \mid \mu_\chi, \sigma_\chi),

        and the cosine tilts are identically but **not independently**
        distributed (Eq. B16),

        .. math::

           \pi(\cos\theta_i \mid \mu_t, \sigma_t, \xi)
             = \xi\, N_{[t_{\min}, 1]}(\cos\theta_1 \mid \mu_t, \sigma_t)\,
                     N_{[t_{\min}, 1]}(\cos\theta_2 \mid \mu_t, \sigma_t)
             + \frac{1 - \xi}{(1 - t_{\min})^2}.

        Parameters, from ``spin-parameters``: ``mu_chi``, ``sigma_chi``,
        ``mu_t``, ``sigma_t``, ``xi_spin``, and optionally ``amax`` (default 1)
        and ``t_min`` (default -1). The posterior files call ``mu_t`` and
        ``sigma_t`` ``mu_spin`` and ``sigma_spin``.

        Returns:
        --------
            dict: ``a_1``, ``a_2``, ``tilt_1``, ``tilt_2``, ``phi_12``, ``phi_jl``.
        """
        logging.info("Generating spins from the GWTC-5.0 Default BBH model")
        missing = [
            name
            for name in ("mu_chi", "sigma_chi", "mu_t", "sigma_t", "xi_spin")
            if name not in self.parameters
        ]
        if missing:
            raise ValueError(
                "The Default spin model needs {}. Note that it takes "
                "'sigma_chi' (a standard deviation) and a free tilt mean "
                "'mu_t'; 'sigma_squared_chi' belongs to the Beta models, which "
                "are still available as 'Isotropic-Beta_Gaussian_Uniform' and "
                "friends.".format(missing)
            )

        count = self.number_of_samples
        amax = self.parameters["amax"]
        t_min = self.parameters["t_min"]
        samples = {}

        magnitude = bilby.core.prior.analytical.TruncatedNormal(
            name="a",
            mu=self.parameters["mu_chi"],
            sigma=self.parameters["sigma_chi"],
            minimum=0.0,
            maximum=amax,
        )
        samples["a_1"] = magnitude.sample(count)
        samples["a_2"] = magnitude.sample(count)

        aligned = bilby.core.prior.analytical.TruncatedNormal(
            name="cos_tilt",
            mu=self.parameters["mu_t"],
            sigma=self.parameters["sigma_t"],
            minimum=t_min,
            maximum=1.0,
        )
        isotropic = bilby.core.prior.analytical.Uniform(
            name="cos_tilt", minimum=t_min, maximum=1.0
        )
        # One draw per *binary*: this is what correlates the two tilts.
        in_gaussian = numpy.random.uniform(size=count) < self.parameters["xi_spin"]
        for index in (1, 2):
            cosines = numpy.where(
                in_gaussian, aligned.sample(count), isotropic.sample(count)
            )
            samples["tilt_{}".format(index)] = numpy.arccos(cosines)

        samples["phi_12"] = bilby.gw.prior.Uniform(
            name="phi_12", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
        ).sample(count)
        samples["phi_jl"] = bilby.gw.prior.Uniform(
            name="phi_jl", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
        ).sample(count)
        return samples

    def sample(self):
        """
        Generate spin distribution samples based on the chosen parameterised model and its parameters.

        Returns:
        --------
            dict: A dictionary containing spin distribution samples.
        """

        samples = {}
        err_msg = "Please choose from one of these spin models {}".format(choices)

        if self.spin_model in [
            utils.remove_special_characters(choice.lower()) for choice in choices
        ]:
            pass
        else:
            raise ValueError(err_msg)
        if self.spin_model == "default":
            return self._sample_default_bbh()
        if "nonspinning" in self.spin_model:
            if "gaussiannonspinning" in self.spin_model:
                mu_chi_1 = self.parameters["mu_chi_1"]
                sigma_chi_1 = self.parameters["sigma_chi_1"]
                chi_1 = bilby.core.prior.analytical.TruncatedGaussian(
                    name="chi_1",
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                    mu=mu_chi_1,
                    sigma=sigma_chi_1,
                )
                samples["chi_1"] = chi_1.sample(self.number_of_samples)
                samples["chi_2"] = numpy.zeros(self.number_of_samples)
            else:
                logging.info("You chose non-spinning distribution")
                samples["chi_1"] = numpy.zeros(self.number_of_samples)
                samples["chi_2"] = numpy.zeros(self.number_of_samples)

        elif "aligned" in self.spin_model:
            logging.info("Generating spin samples from Aligned Spin Distribution")
            if "bilby" in self.spin_model or "aligned" == self.spin_model:
                chi_1 = bilby.gw.prior.AlignedSpin(
                    name="chi_1",
                    a_prior=bilby.gw.prior.Uniform(
                        minimum=self.parameters["minimum_primary_spin"],
                        maximum=self.parameters["maximum_primary_spin"],
                    ),
                )
                chi_2 = bilby.gw.prior.AlignedSpin(
                    name="chi_2",
                    a_prior=bilby.gw.prior.Uniform(
                        minimum=self.parameters["minimum_secondary_spin"],
                        maximum=self.parameters["maximum_secondary_spin"],
                    ),
                )
            elif "gaussianuniform" in self.spin_model:
                mu_chi_1 = self.parameters["mu_chi_1"]
                sigma_chi_1 = self.parameters["sigma_chi_1"]
                chi_1 = bilby.core.prior.analytical.TruncatedGaussian(
                    name="chi_1",
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                    mu=mu_chi_1,
                    sigma=sigma_chi_1,
                )

                chi_2 = bilby.gw.prior.Uniform(
                    name="chi_2",
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                )
            elif "uniform" in self.spin_model:
                chi_1 = bilby.gw.prior.Uniform(
                    name="chi_1",
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                )
                chi_2 = bilby.gw.prior.Uniform(
                    name="chi_2",
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                )
            elif "beta" in self.spin_model:
                chi_1 = bilby.gw.prior.AlignedSpin(
                    name="chi_1",
                    a_prior=bilby.core.prior.analytical.Beta(
                        name="chi_1",
                        alpha=self.parameters["alpha_chi"],
                        beta=self.parameters["beta_chi"],
                        minimum=self.parameters["minimum_primary_spin"],
                        maximum=self.parameters["maximum_primary_spin"],
                    ),
                )
                chi_2 = bilby.gw.prior.AlignedSpin(
                    name="chi_2",
                    a_prior=bilby.core.prior.analytical.Beta(
                        name="chi_2",
                        alpha=self.parameters["alpha_chi"],
                        beta=self.parameters["beta_chi"],
                        minimum=self.parameters["minimum_secondary_spin"],
                        maximum=self.parameters["maximum_secondary_spin"],
                    ),
                )

            samples["chi_1"] = chi_1.sample(self.number_of_samples)
            samples["chi_2"] = chi_2.sample(self.number_of_samples)

        elif "isotropic" in self.spin_model:
            logging.info("Generating spin samples from Isotropic Spin Distribution")
            if "bilby" in self.spin_model:
                a_1 = bilby.gw.prior.Uniform(
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                    name="a_1",
                )
                a_2 = bilby.gw.prior.Uniform(
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                    name="a_2",
                )

                tilt_1 = bilby.core.prior.analytical.Sine(name="tilt_1")
                tilt_2 = bilby.core.prior.analytical.Sine(name="tilt_2")
                samples["tilt_1"] = tilt_1.sample(self.number_of_samples)
                samples["tilt_2"] = tilt_2.sample(self.number_of_samples)

            elif "betagaussianuniform" in self.spin_model:
                n = int(self.parameters["xi_spin"] * self.number_of_samples)
                m = self.number_of_samples - n
                a_1 = bilby.core.prior.analytical.Beta(
                    name="a_1",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                )

                a_2 = bilby.core.prior.analytical.Beta(
                    name="a_2",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                )

                cos_tilt_gaussian = bilby.core.prior.analytical.TruncatedNormal(
                    name="cos_tilt_gaussian",
                    mu=1,
                    sigma=self.parameters["sigma_t"],
                    minimum=-1,
                    maximum=1,
                )

                cos_tilt_isotropic = bilby.core.prior.analytical.Uniform(
                    name="cos_tilt_isotropic", minimum=-1, maximum=1
                )

                samples["tilt_1"] = numpy.arccos(
                    numpy.concatenate(
                        [cos_tilt_gaussian.sample(n), cos_tilt_isotropic.sample(m)]
                    )
                )
                samples["tilt_2"] = numpy.arccos(
                    numpy.concatenate(
                        [cos_tilt_gaussian.sample(n), cos_tilt_isotropic.sample(m)]
                    )
                )
                # Shuffle so the gaussian/isotropic block ordering does not correlate
                # tilt_1 with tilt_2. Use numpy's RNG so a single global seed governs it.
                numpy.random.shuffle(samples["tilt_1"])
                numpy.random.shuffle(samples["tilt_2"])
            elif "betagaussian" in self.spin_model:
                a_1 = bilby.core.prior.analytical.Beta(
                    name="a_1",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                )

                a_2 = bilby.core.prior.analytical.Beta(
                    name="a_2",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                )

                cos_tilt = bilby.core.prior.analytical.TruncatedNormal(
                    name="cos_tilt",
                    mu=1,
                    sigma=self.parameters["sigma_t"],
                    minimum=-1,
                    maximum=1,
                )

                samples["tilt_1"] = numpy.arccos(
                    cos_tilt.sample(self.number_of_samples)
                )
                samples["tilt_2"] = numpy.arccos(
                    cos_tilt.sample(self.number_of_samples)
                )

            elif "beta" in self.spin_model:
                a_1 = bilby.core.prior.analytical.Beta(
                    name="a_1",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_primary_spin"],
                    maximum=self.parameters["maximum_primary_spin"],
                )
                a_2 = bilby.core.prior.analytical.Beta(
                    name="a_2",
                    alpha=self.parameters["alpha_chi"],
                    beta=self.parameters["beta_chi"],
                    minimum=self.parameters["minimum_secondary_spin"],
                    maximum=self.parameters["maximum_secondary_spin"],
                )
                tilt_1 = bilby.core.prior.analytical.Sine(name="tilt_1")
                tilt_2 = bilby.core.prior.analytical.Sine(name="tilt_2")
                samples["tilt_1"] = tilt_1.sample(self.number_of_samples)
                samples["tilt_2"] = tilt_2.sample(self.number_of_samples)

            phi_12 = bilby.gw.prior.Uniform(
                name="phi_12", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
            )
            phi_jl = bilby.gw.prior.Uniform(
                name="phi_jl", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
            )

            samples["a_1"] = a_1.sample(self.number_of_samples)
            samples["a_2"] = a_2.sample(self.number_of_samples)
            samples["phi_12"] = phi_12.sample(self.number_of_samples)
            samples["phi_jl"] = phi_jl.sample(self.number_of_samples)
        return samples
