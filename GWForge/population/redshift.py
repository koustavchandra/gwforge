import numpy, sympy, logging, bilby
from .. import utils
from sympy import symbols, lambdify, integrate
from lal import YRJUL_SI, PC_SI
from scipy.interpolate import interp1d
from scipy.integrate import quad
from pycbc.population import population_models
from rich.progress import track
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S')


class Redshift:
    def __init__(self,
                 redshift_model, 
                 local_merger_rate_density,
                 maximum_redshift, 
                 gps_start_time,
                 analysis_time = YRJUL_SI,
                 cosmology = 'Planck18', 
                 parameters = {'gamma' : 2.7, 'kappa': 5.6, 'z_peak':1.9},
                 time_delay_model = 'inverse',
                 H0=70, Om0=0.3, Ode0=0.7, Tcmb0=2.725, Ob0=None):

        '''
        Initialise the Redshift class

        Parameters:
        ----------
        redshift_model : str
            Implemented gwpopulation redshift model to use
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
            Implemented PyCBC time delay model to use [Default: inverse]
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
        '''
        
        self.redshift_model = utils.remove_special_characters(redshift_model.lower())
        self.local_merger_rate_density = local_merger_rate_density * 1e-9
        self.maximum_redshift = maximum_redshift
        self.analysis_time = analysis_time
        self.gps_start_time = gps_start_time
        self.cosmology = cosmology
        self.parameters = parameters
        self.time_delay_model = time_delay_model
        self.H0, self.Om0, self.Ode0, self.Tcmb0, self.Ob0 = H0, Om0, Ode0, Tcmb0, Ob0

    def import_cosmology(self):
        import importlib
        try:
            from astropy import cosmology
            cosmology = getattr(cosmology, self.cosmology)
            return cosmology
        except (ImportError, AttributeError):
            if self.H0 is not None and self.Om0 is not None and self.Ode0 is not None:
                logging.info('Importing FLRW Cosmology with the provided constants')
                from astropy.cosmology import LambdaCDM
                return LambdaCDM(H0=self.H0, Om0=self.Om0, Ode0=self.Ode0, Tcmb0=self.Tcmb0, Ob0=self.Ob0)
            else:
                raise ValueError("Could not import cosmology")

    def differential_lookback_time(self, redshift):
        '''
        Derivative of lookback time t(z) 
        See Eq.(A3) of <arXiv:2011.02717v3>
        '''
        cosmo = self.import_cosmology()
        H0 = cosmo.H0.value / (PC_SI * 1e3) # H0 in seconds = km --> Mpc 
        H0 = H0 * (YRJUL_SI / 1e-9)  # H0 in Gyr^-1
        dz_dt = H0 * (1+redshift) * sympy.sqrt((cosmo.Ode0 + cosmo.Om0 * (1+redshift)**3))
        return 1/dz_dt # returns dt_dz 

    def transform(self):
        '''
        Adapted from pycbc.population.population_models
        
        Note:
        
        This function combines the star formation rate, time delay probability, and differential lookback time
        to compute the merger rate density. It uses the specified cosmological model.    
        '''
        from functools import partial
        
        z, z_peak = symbols('z'), symbols('z_0')
        diff_lookback_time = partial(self.differential_lookback_time)
        time_delay = integrate(diff_lookback_time(redshift=z), 
                               (z, z_peak, z))
        
        if self.redshift_model == 'madaudickinson':
            from gwpopulation.models.redshift import MadauDickinsonRedshift
            # Note the GWPopulation misses out the 0.015 constant in psi_of_z
            psi_of_z = 0.015 * MadauDickinsonRedshift(z_max=self.maximum_redshift).psi_of_z(redshift=z, **self.parameters) 
        elif self.redshift_model == 'powerlaw':
            from gwpopulation.models.redshift import PowerLawRedshift
            psi_of_z = PowerLawRedshift(z_max=self.maximum_redshift).psi_of_z(redshift=z, **self.parameters)
        else:
            raise ValueError('Redshift model {} is not implemented in GWPopulation')
    
        return psi_of_z * population_models.p_tau(tau=time_delay, 
                                          td_model=self.time_delay_model) * diff_lookback_time(z)


    def rate_density(self, elements = 1000):
        '''
        Adapted from pycbc.population.population_models.
        
        Return:
        -------
        merger rate density : scipy.interpolate.interp1d
        '''
        redshift = numpy.linspace(0, self.maximum_redshift, elements)
    
        z, z_0 = symbols('z'), symbols('z_0')
        merger_rate_density = numpy.zeros_like(redshift)
        
        function = self.transform() 
        
        for k in track(range(len(redshift))):
            function_2 = lambdify(z, function.subs(z_0, redshift[k]), 'scipy')
            merger_rate_density[k] = quad(function_2, redshift[k], numpy.inf, epsabs=1e-3)[0]
        merger_rate_density = merger_rate_density / merger_rate_density[0] * self.local_merger_rate_density  # Normalize & Rescale
        return interp1d(redshift, merger_rate_density)

    def average_time_between_signals(self):
        '''
        Calculates the average time interval between two signals of same type
        '''
        rate_density = self.rate_density()
        # this PyCBC function returns an interpolation function ---> merger rate at maximum redshift    
        merger_rate = population_models.coalescence_rate(rate_den=rate_density, maxz=self.maximum_redshift)
        # The following returns the average time delay. See Eq. (7) of <Phys. Rev. D 93, 024018 (2016)>
        return 1 / population_models.total_rate_upto_redshift(z=self.maximum_redshift, merger_rate=merger_rate) * YRJUL_SI    

    def sample(self):
        '''
        Return:
        -------
        parameters : dict
            dictionary of redshift, time_interval, tc
        '''
        if self.redshift_model == 'madaudickinson':
            logging.info('Generating samples assumming Madau-Dickinson Model')
            from gwpopulation.models.redshift import MadauDickinsonRedshift
            model = MadauDickinsonRedshift(z_max=self.maximum_redshift)
        elif self.redshift_model == 'powerlaw':
            logging.info('Generating samples Power-Law Model')
            from gwpopulation.models.redshift import PowerLawRedshift
            model = PowerLawRedshift(z_max=self.maximum_redshift)
        else:
            raise ValueError('Redshift model {} is not implemented in GWPopulation')
        dataset = {'redshift' : model.zs}
        prob = model.probability(dataset, **self.parameters)
        prior = bilby.core.prior.Interped(xx=dataset['redshift'], yy=prob, 
                                          minimum=0.,maximum=self.maximum_redshift, name='redshift')
        average_time_interval = self.average_time_between_signals()
        logging.info('Average time interval between signals = {:.2f}'.format(average_time_interval))
        number_of_samples = int(self.analysis_time / average_time_interval)
        logging.info('Number of samples generated = {}'.format(number_of_samples))
        z = prior.sample(number_of_samples)
        
        time_gap = bilby.prior.analytical.Exponential(mu=average_time_interval)
        time_interval = time_gap.sample(number_of_samples)
        tc = time_interval.cumsum() + self.gps_start_time
        
        return {'redshift' : z, 'time_interval': time_interval, 'geocent_time': tc}