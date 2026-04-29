import numpy
import sympy
import logging
import bilby
from .. import utils
from sympy import symbols, lambdify, integrate
from lal import YRJUL_SI, PC_SI
from scipy.interpolate import interp1d
from scipy.integrate import quad
from rich.progress import track

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")


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
        td_min=0.02,
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
            Implemented gwpopulation redshift model to use. Options: 'MadauDickinson', 'PowerLaw'.
            For 'PowerLaw', parameters should contain {'lamb': value} (power-law index).
            For 'MadauDickinson', parameters should contain {'gamma', 'kappa', 'z_peak'}.
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
            Dictionary of parameters for the chosen redshift model
        time_delay_model : str, optional
            Time delay model to use [Default: inverse]. Options: 'inverse', 'log_normal', 'gaussian', 'power_law'
        td_min : float, optional
            Minimum time delay in Gyr used in the inverse time delay model [Default: 0.02]
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
        self.time_delay_model = time_delay_model
        self.td_min = td_min
        self.H0, self.Om0, self.Ode0, self.Tcmb0, self.Ob0 = H0, Om0, Ode0, Tcmb0, Ob0

    def import_cosmology(self):

        try:
            from astropy import cosmology

            cosmology = getattr(cosmology, self.cosmology)
            return cosmology
        except (ImportError, AttributeError):
            if self.H0 is not None and self.Om0 is not None and self.Ode0 is not None:
                logging.info("Importing FLRW Cosmology with the provided constants")
                from astropy.cosmology import LambdaCDM

                return LambdaCDM(
                    H0=self.H0,
                    Om0=self.Om0,
                    Ode0=self.Ode0,
                    Tcmb0=self.Tcmb0,
                    Ob0=self.Ob0,
                )
            else:
                raise ValueError("Could not import cosmology")

    def differential_lookback_time(self, redshift):
        """
        Derivative of lookback time t(z)
        See Eq.(A3) of <arXiv:2011.02717v3>
        """
        cosmo = self.import_cosmology()
        H0 = cosmo.H0.value / (PC_SI * 1e3)  # H0 in seconds = km --> Mpc
        H0 = H0 * (YRJUL_SI / 1e-9)  # H0 in Gyr^-1
        dz_dt = H0 * (1 + redshift) * sympy.sqrt((cosmo.Ode0 + cosmo.Om0 * (1 + redshift) ** 3))
        return 1 / dz_dt  # returns dt_dz

    def p_tau(self, tau):
        """
        The probability distribution of the time delay.
        Adapted from pycbc.population.population_models.p_tau

        Parameters
        ----------
        tau : float, numpy.ndarray, or sympy expression
            The merger delay time in Gyr.

        Returns
        -------
        p_t : float, numpy.ndarray, or sympy expression
            The probability at time delay tau.

        Notes
        -----
            See the Appendix in <arXiv:2011.02717v3> for more details.
        """
        from sympy import sqrt, exp, log, Piecewise

        td_model = self.time_delay_model
        td_min = self.td_min
        td_max = self.import_cosmology().age(0).to("Gyr").value

        if td_model == "log_normal":
            t_ln = 2.9  # Gyr
            sigma_ln = 0.2
            p_t = exp(-(log(tau) - log(t_ln)) ** 2 / (2 * sigma_ln ** 2)) / (sqrt(2 * numpy.pi) * sigma_ln)
        elif td_model == "gaussian":
            t_g = 2  # Gyr
            sigma_g = 0.3
            p_t = exp(-(tau - t_g) ** 2 / (2 * sigma_g ** 2)) / (sqrt(2 * numpy.pi) * sigma_g)
        elif td_model == "power_law":
            alpha_t = 0.81
            p_t = tau ** (-alpha_t)
        elif td_model == "inverse":
            if isinstance(tau, (float, int)) or isinstance(tau, numpy.ndarray):
                norm_const = 1 / numpy.log(td_max / td_min)
                p_t = numpy.where((tau < td_min) | (tau > td_max), 0, norm_const * tau ** (-0.999))
            else:
                norm_const = 1 / numpy.log(td_max / td_min)
                p_t = Piecewise((0, tau < td_min), (0, tau > td_max), (norm_const * tau ** (-0.999), True))
        else:
            raise ValueError("'time_delay_model' must be one of ['log_normal', 'gaussian', 'power_law', 'inverse'].")
        return p_t

    def transform(self):
        """
        Adapted from pycbc.population.population_models

        Note:

        This function combines the star formation rate, time delay probability, and differential lookback time
        to compute the merger rate density. It uses the specified cosmological model.
        """
        from functools import partial

        z, z_peak = symbols("z"), symbols("z_0")
        diff_lookback_time = partial(self.differential_lookback_time)
        time_delay = integrate(diff_lookback_time(redshift=z), (z, z_peak, z))

        if self.redshift_model == "madaudickinson":
            from gwpopulation.models.redshift import MadauDickinsonRedshift

            # Note the GWPopulation misses out the 0.015 constant in psi_of_z
            psi_of_z = 0.015 * MadauDickinsonRedshift(z_max=self.maximum_redshift).psi_of_z(redshift=z, **self.parameters)
        elif self.redshift_model == "powerlaw":
            from gwpopulation.models.redshift import PowerLawRedshift

            psi_of_z = PowerLawRedshift(z_max=self.maximum_redshift).psi_of_z(redshift=z, **self.parameters)
        else:
            raise ValueError(f"Redshift model {self.redshift_model} is not implemented in GWPopulation")

        return psi_of_z * self.p_tau(tau=time_delay) * diff_lookback_time(z)

    def rate_density(self, elements=1000):
        """
        Adapted from pycbc.population.population_models.

        Return:
        -------
        merger rate density : scipy.interpolate.interp1d
        """
        redshift = numpy.linspace(0, self.maximum_redshift, elements)

        z, z_0 = symbols("z"), symbols("z_0")
        merger_rate_density = numpy.zeros_like(redshift)

        function = self.transform()

        for k in track(range(len(redshift))):
            function_2 = lambdify(z, function.subs(z_0, redshift[k]), "scipy")
            merger_rate_density[k] = quad(function_2, redshift[k], numpy.inf, epsabs=1e-3)[0]
        merger_rate_density = merger_rate_density / merger_rate_density[0] * self.local_merger_rate_density  # Normalize & Rescale
        return interp1d(redshift, merger_rate_density)

    def coalescence_rate(self):
        """
        Calculates the merger rate from rate_density function.
        Adapted from pycbc.population.population_models.coalescence_rate
        """
        from astropy import units

        rate_density = self.rate_density()
        cosmology = self.import_cosmology()

        z_array = numpy.linspace(0, self.maximum_redshift, 1000)
        dr_dz = []
        for z in z_array:
            dr = cosmology.differential_comoving_volume(z) / (1 + z)
            dr_dz.append((dr * 4 * numpy.pi * units.sr * rate_density(z) * (units.Mpc) ** (-3)).value)

        return interp1d(z_array, dr_dz, fill_value="extrapolate")

    def average_time_between_signals(self):
        """
        Calculates the average time interval between two signals of same type.
        Adapted from pycbc.population.population_models.total_rate_upto_redshift
        """
        merger_rate = self.coalescence_rate()
        total_rate = quad(merger_rate, 0, self.maximum_redshift, epsabs=2.00e-4, epsrel=2.00e-4, limit=1000)[0]
        return 1 / total_rate * YRJUL_SI

    def sample(self):
        """
        Return:
        -------
        parameters : dict
            dictionary of redshift, time_interval, tc
        """
        if self.redshift_model == "madaudickinson":
            logging.info("Generating samples assumming Madau-Dickinson Model")
            from gwpopulation.models.redshift import MadauDickinsonRedshift

            model = MadauDickinsonRedshift(z_max=self.maximum_redshift)
        elif self.redshift_model == "powerlaw":
            logging.info("Generating samples Power-Law Model")
            from gwpopulation.models.redshift import PowerLawRedshift

            model = PowerLawRedshift(z_max=self.maximum_redshift)
        else:
            raise ValueError("Redshift model {} is not implemented in GWPopulation")

        rate = self.coalescence_rate()
        dataset = {"redshift": model.zs}
        prob = rate(dataset["redshift"])
        prior = bilby.core.prior.Interped(
            xx=dataset["redshift"],
            yy=prob,
            minimum=0.0,
            maximum=self.maximum_redshift,
            name="redshift",
        )
        average_time_interval = self.average_time_between_signals()
        logging.info("Average time interval between signals = {:.2f}".format(average_time_interval))
        number_of_samples = int(self.analysis_time / average_time_interval)
        logging.info("Number of samples generated = {}".format(number_of_samples))
        z = prior.sample(number_of_samples)

        time_gap = bilby.prior.analytical.Exponential(mu=average_time_interval)
        time_interval = time_gap.sample(number_of_samples)
        tc = time_interval.cumsum() + self.gps_start_time

        return {"redshift": z, "time_interval": time_interval, "geocent_time": tc}
