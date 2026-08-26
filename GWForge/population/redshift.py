import numpy
import logging
import bilby
from .. import utils
from ..cosmology import astropy_cosmology, differential_comoving_volume
from lal import YRJUL_SI, PC_SI
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Supported time-delay models.
TIME_DELAY_MODELS = ("inverse", "powerlaw")
# Regimbau et al. 2012 said use 20 Myr (0.02 Gyr)
# Phys. Rev. D 86, 122001), below which p(tau) = 0.
DEFAULT_TIME_DELAY_MINIMUM = 0.02
# Merger rate density at z = 0 in Gpc^-3 yr^-1: the posterior median of the O4b
# `Default` analysis (arXiv:2605.27226). One should get roughly 33000 BBHs if z_max = 10
GWTC5_LOCAL_MERGER_RATE_DENSITY = 23.3895396
# Lookback-time table resolution. The table must reach the age of the universe,
# so it extends to very high redshift (lookback time saturates as z -> inf).
FORMATION_GRID_POINTS = 20000
FORMATION_REDSHIFT_MAXIMUM = 1000.0
# Number of points on the log-spaced time-delay integration grid.
TIME_DELAY_GRID_POINTS = 3000
# Number of points on the coalescence-rate grid.
COALESCENCE_GRID_POINTS = 10000

# Hyperparameter names each GWPopulation redshift model's psi_of_z accepts. Used to
# pass only the relevant parameters so that, e.g., a PowerLaw model is not handed the
# Madau-Dickinson {gamma, kappa, z_peak} defaults (which would raise a TypeError).
REDSHIFT_MODEL_PARAMETERS = {
    "madaudickinson": ("gamma", "kappa", "z_peak"),
    "powerlaw": ("lamb",),
}


def madau_dickinson_psi_of_z(redshift, gamma, kappa, z_peak):
    r"""Madau-Dickinson redshift evolution :math:`\psi(z)`, normalised so that
    :math:`\psi(0) = 1`.

    .. math::
        \psi(z) = \frac{(1 + z)^\gamma}{1 + ((1 + z) / (1 + z_p))^\kappa}
                  \left[1 + (1 + z_p)^{-\kappa}\right]

    Matches ``gwpopulation.models.redshift.MadauDickinsonRedshift.psi_of_z``
    (Fishbach+ arXiv:1805.10270 Eq. 33; normalisation from Callister+
    arXiv:2003.12152 Eq. 2). Written with plain arithmetic so it also accepts a
    ``sympy`` symbol, as required by :meth:`Redshift.transform`.
    """
    psi = (1 + redshift) ** gamma / (1 + ((1 + redshift) / (1 + z_peak)) ** kappa)
    psi *= 1 + (1 + z_peak) ** (-kappa)
    return psi


def power_law_psi_of_z(redshift, lamb):
    r"""Power-law redshift evolution :math:`\psi(z) = (1 + z)^\lambda`.

    Matches ``gwpopulation.models.redshift.PowerLawRedshift.psi_of_z``
    (Fishbach+ arXiv:1805.10270).
    """
    return (1 + redshift) ** lamb


# Redshift grid used to build the merger-rate interpolant, matching the
# gwpopulation ``_Redshift.zs`` convention (linspace(1e-6, z_max, 2500)).
REDSHIFT_GRID_POINTS = 2500
REDSHIFT_GRID_MINIMUM = 1e-6


class Redshift:
    def __init__(
        self,
        redshift_model,
        local_merger_rate_density,
        maximum_redshift,
        gps_start_time,
        analysis_time=YRJUL_SI,
        cosmology="Planck18",
        parameters={"gamma": 2.7, "kappa": 5.6, "z_peak": 1.9},
        time_delay_model="inverse",
        time_delay_slope=1.0,
        time_delay_minimum=DEFAULT_TIME_DELAY_MINIMUM,
        H0=70,
        Om0=0.3,
        Ode0=0.7,
        Tcmb0=2.725,
        Ob0=None,
    ):
        """
        Initialise the Redshift class

        Parameters:
        ----------
        redshift_model : str
            Redshift evolution model to use ("MadauDickinson" or "PowerLaw")
        local_merger_rate_density : float
            Estimated local merger rate density in /Gpc^3/yr
        maximum_redshift: float
            Maximum redshift for calculation
        gps_start_time : float
            GPS start time of analysis
        analysis_time : float
            Total analysis time [Default: 1 year]
        cosmology: str, optional
            Name of astropy cosmology class to use [Default: Planck18]
        parameters: dict, optional
            Dictionary of parameters
        time_delay_model : str, optional
            Formation->merger time-delay model. One of {"inverse", "powerlaw"}.
            Both use p(tau) ∝ tau^-slope; "inverse" fixes slope = 1 (the 1/tau
            model). [Default: inverse]
        time_delay_slope : float, optional
            Power-law slope of the time-delay distribution p(tau) ∝ tau^-slope.
            Ignored (forced to 1) when time_delay_model == "inverse". [Default: 1]
        time_delay_minimum : float, optional
            Minimum formation->merger delay in Gyr below which p(tau) = 0.
            [Default: 0.02]
        H0 : float, optional
            Hubble constant value for custom LDCM cosmology [Default: 70]
        Om0 : float, optional
            Matter density parameter value for custom LDCM cosmology [Default: 0.3]
        Ode0 : float, optional
            Dark energy density parameter value for custom LDCM cosmology [Default: 0.7]
        Tcmb0: float, optional
            Temperature of the CMB z=0 [Default: 2.725]
        Ob0: float, optional
            Density of baryonic matter in units of the critical density at z=0. [Default: None]
        """

        self.redshift_model = utils.remove_special_characters(redshift_model.lower())
        self.local_merger_rate_density = local_merger_rate_density * 1e-9
        self.maximum_redshift = maximum_redshift
        self.analysis_time = analysis_time
        self.gps_start_time = gps_start_time
        self.cosmology = cosmology
        self.parameters = parameters
        self.time_delay_model = utils.remove_special_characters(
            time_delay_model.lower()
        )
        if self.time_delay_model not in TIME_DELAY_MODELS:
            raise ValueError(
                "time_delay_model must be one of {} (got '{}')".format(
                    TIME_DELAY_MODELS, time_delay_model
                )
            )
        # The "inverse" (1/tau) model is the slope == 1 special case.
        self.time_delay_slope = (
            1.0 if self.time_delay_model == "inverse" else time_delay_slope
        )
        self.time_delay_minimum = time_delay_minimum
        self.H0, self.Om0, self.Ode0, self.Tcmb0, self.Ob0 = H0, Om0, Ode0, Tcmb0, Ob0

    def import_cosmology(self):
        """Resolve this object's cosmology specification.

        Delegates to :func:`GWForge.cosmology.astropy_cosmology`, which is the
        single place a GWForge cosmology specification is resolved, so that the
        redshifts a population is sampled from and the luminosity distances it
        is labelled with cannot come from two different universes.

        Returns
        -------
        astropy.cosmology.FLRW
        """
        return astropy_cosmology(
            name=self.cosmology,
            H0=self.H0,
            Om0=self.Om0,
            Ode0=self.Ode0,
            Tcmb0=self.Tcmb0,
            Ob0=self.Ob0,
        )

    def differential_lookback_time(self, redshift):
        """
        Derivative of lookback time t(z) with respect to z, in Gyr.
        See Eq.(A3) of <arXiv:2011.02717v3>. Vectorised over ``redshift``.
        """
        cosmo = self.import_cosmology()
        H0 = cosmo.H0.value / (PC_SI * 1e3)  # H0 in 1/s (km/s/Mpc -> 1/s)
        H0 = H0 * (YRJUL_SI / 1e-9)  # H0 in Gyr^-1
        redshift = numpy.asarray(redshift, dtype=float)
        dz_dt = (
            H0
            * (1 + redshift)
            * numpy.sqrt(cosmo.Ode0 + cosmo.Om0 * (1 + redshift) ** 3)
        )
        return 1.0 / dz_dt  # returns dt/dz

    def _psi_of_z(self, redshift):
        """Star-formation-rate redshift evolution psi(z) for the chosen model.

        The overall amplitude is irrelevant here: :meth:`rate_density` rescales to
        ``local_merger_rate_density`` at z = 0, so any constant prefactor cancels.
        """
        if self.redshift_model == "madaudickinson":
            return madau_dickinson_psi_of_z(redshift, **self._model_parameters())
        elif self.redshift_model == "powerlaw":
            return power_law_psi_of_z(redshift, **self._model_parameters())
        raise ValueError(
            "Redshift model {} is not implemented".format(self.redshift_model)
        )

    def _maximum_time_delay(self):
        """Maximum formation->merger delay (Gyr): the age of the universe at z=0."""
        cosmo = self.import_cosmology()
        return cosmo.age(0).to("Gyr").value

    def time_delay_probability(self, tau):
        r"""Formation->merger time-delay distribution, :math:`p(\tau)`.

        Power-law model :math:`p(\tau) \propto \tau^{-\mathrm{slope}}` supported on
        ``[time_delay_minimum, age_of_universe]`` and zero elsewhere. ``slope = 1``
        recovers the ``1/\tau`` ("inverse") model. The normalisation is a global
        constant and cancels in :meth:`rate_density`, so it is omitted here.
        """
        tau = numpy.asarray(tau, dtype=float)
        support = (tau >= self.time_delay_minimum) & (tau <= self._maximum_time_delay())
        prob = numpy.zeros_like(tau)
        numpy.power(tau, -self.time_delay_slope, where=support, out=prob)
        return numpy.where(support, prob, 0.0)

    def _model_parameters(self):
        """
        Select only the hyperparameters the chosen redshift model's psi_of_z accepts,
        raising a clear error if a required one is missing.
        """
        expected = REDSHIFT_MODEL_PARAMETERS[self.redshift_model]
        missing = [p for p in expected if p not in self.parameters]
        if missing:
            raise ValueError(
                "Redshift model '{}' requires parameter(s) {}. Provide them via "
                "'redshift-parameters' (got {}).".format(
                    self.redshift_model, missing, list(self.parameters)
                )
            )
        return {p: self.parameters[p] for p in expected}

    def rate_density(self, elements=1000):
        r"""Merger-rate density :math:`R(z_m)` as a function of merger redshift.

        Convolves the star-formation rate with the time-delay distribution,

        .. math::
            R(z_m) \propto \int_{z_m}^{\infty} \psi(z_f)\, p\!\left(t(z_f) - t(z_m)\right)
                          \frac{dt}{dz_f}\, dz_f,

        where :math:`\psi` is the SFR, :math:`p(\tau)` the time-delay distribution
        (:meth:`time_delay_probability`), and :math:`t(z)` the lookback time. The
        result is normalised and rescaled to ``local_merger_rate_density`` at
        :math:`z_m = 0`. Replaces the previous pycbc + sympy implementation with
        direct vectorised numerical integration.

        Return:
        -------
        merger rate density : scipy.interpolate.interp1d
            Callable over ``[0, maximum_redshift]``, in Mpc^-3 yr^-1.
        """
        merger_redshift = numpy.linspace(0, self.maximum_redshift, elements)

        # Lookback-time table t(z) used to convert between formation redshift and
        # time delay. It must reach up to the age of the universe, so the grid
        # extends to very high redshift (t(z) saturates at the age as z -> inf).
        formation_redshift = numpy.linspace(
            0, FORMATION_REDSHIFT_MAXIMUM, FORMATION_GRID_POINTS
        )
        dt_dz = self.differential_lookback_time(formation_redshift)  # Gyr
        lookback_time = cumulative_trapezoid(dt_dz, formation_redshift, initial=0.0)
        lookback_at_merger = numpy.interp(
            merger_redshift, formation_redshift, lookback_time
        )

        # Integrate over the time delay tau rather than formation redshift. The
        # convolution R(z_m) = \int SFR(z_f(tau)) p(tau) dtau is evaluated on a
        # log-spaced tau grid: with p(tau) ∝ tau^-slope the integrand in
        # d ln(tau) is SFR(z_f) * tau^(1-slope), which is smooth (flat for the
        # 1/tau "inverse" model) and free of the sharp short-delay peak that a
        # uniform formation-redshift grid fails to resolve.
        tau_maximum = self._maximum_time_delay()
        ln_tau = numpy.linspace(
            numpy.log(self.time_delay_minimum),
            numpy.log(tau_maximum),
            TIME_DELAY_GRID_POINTS,
        )
        tau = numpy.exp(ln_tau)

        formation_lookback = lookback_at_merger[:, None] + tau[None, :]
        # Formation redshift for each (z_m, tau): invert the lookback table.
        formation_z = numpy.interp(
            formation_lookback.ravel(), lookback_time, formation_redshift
        ).reshape(formation_lookback.shape)
        # Delays that would place formation before the start of the grid (t > t(z_max))
        # contribute nothing (SFR ~ 0 there for physical models).
        within = formation_lookback <= lookback_time[-1]
        sfr = self._psi_of_z(formation_z)
        integrand = sfr * tau[None, :] ** (1.0 - self.time_delay_slope) * within
        merger_rate_density = numpy.trapezoid(integrand, ln_tau, axis=1)

        merger_rate_density = (
            merger_rate_density
            / merger_rate_density[0]
            * self.local_merger_rate_density
        )
        return interp1d(merger_redshift, merger_rate_density)

    def coalescence_rate(self):
        r"""Detector-frame merger rate per unit redshift, :math:`dR/dz` (yr^-1).

        .. math::
            \frac{dR}{dz} = R(z)\, \frac{1}{1 + z}\, \frac{dV_c}{dz},

        with :math:`dV_c/dz = 4\pi\, (dV_c/d\Omega/dz)` the comoving volume element
        and the :math:`1/(1+z)` factor converting source-frame to detector-frame
        time. Reproduces ``pycbc.population.population_models.coalescence_rate``.
        """
        rate_density = self.rate_density()
        cosmo = self.import_cosmology()
        z = numpy.linspace(0, self.maximum_redshift, COALESCENCE_GRID_POINTS)
        dr_dz = (
            rate_density(z) * differential_comoving_volume(cosmo, z) / (1 + z)
        )
        return interp1d(z, dr_dz, fill_value="extrapolate")

    def average_time_between_signals(self):
        """
        Calculates the average time interval (in seconds) between two signals of
        the same type, out to ``maximum_redshift``.
        """
        merger_rate = self.coalescence_rate()
        z = numpy.linspace(0, self.maximum_redshift, COALESCENCE_GRID_POINTS)
        # Total rate (yr^-1) out to maximum_redshift; see Eq. (7) of
        # <Phys. Rev. D 93, 024018 (2016)>.
        total_rate = numpy.trapezoid(merger_rate(z), z)
        return 1.0 / total_rate * YRJUL_SI

    def sample(self):
        """
        Return:
        -------
        parameters : dict
            dictionary of redshift, time_interval, tc
        """
        if self.redshift_model == "madaudickinson":
            logging.info("Generating samples assumming Madau-Dickinson Model")
        elif self.redshift_model == "powerlaw":
            logging.info("Generating samples Power-Law Model")
        else:
            raise ValueError(
                "Redshift model {} is not implemented".format(self.redshift_model)
            )

        rate = self.coalescence_rate()
        zs = numpy.linspace(
            REDSHIFT_GRID_MINIMUM, self.maximum_redshift, REDSHIFT_GRID_POINTS
        )
        dataset = {"redshift": zs}
        prob = rate(dataset["redshift"])
        prior = bilby.core.prior.Interped(
            xx=dataset["redshift"],
            yy=prob,
            minimum=0.0,
            maximum=self.maximum_redshift,
            name="redshift",
        )
        average_time_interval = self.average_time_between_signals()
        logging.info(
            "Average time interval between signals = {:.2f}".format(
                average_time_interval
            )
        )
        number_of_samples = int(self.analysis_time / average_time_interval)
        logging.info("Number of samples generated = {}".format(number_of_samples))
        z = prior.sample(number_of_samples)

        time_gap = bilby.core.prior.analytical.Exponential(mu=average_time_interval)
        time_interval = time_gap.sample(number_of_samples)
        tc = time_interval.cumsum() + self.gps_start_time

        return {"redshift": z, "time_interval": time_interval, "geocent_time": tc}
