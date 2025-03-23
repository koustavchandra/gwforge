import numpy, logging, bilby, random
from .. import utils
from ..conversion import *

logging.basicConfig(level = logging.INFO,
                    format = '%(asctime)s %(message)s',
                    datefmt = '%Y-%m-%d %H:%M:%S')

choices = ['Non-spinning', 
           'Aligned', 'Aligned-Bilby','Aligned-Uniform', 'Beta-Aligned', 'Aligned-Gaussian-Uniform',  'Isotropic-Bilby', 'Isotropic-Beta', 'Isotropic-Beta_Gaussian', 'Isotropic-Beta_Gaussian_Uniform', 
           'Default']

class Spin:
    def __init__(self, 
                 spin_model,
                 number_of_samples,
                 parameters={}):
        '''
        Parameters:
        ----------
        spin_model : str 
            The parameterized spin model. [Options: {}]
        number_of_samples : (int)
            The number of samples to generate. [Ideal: Exactly same as redshift samples]
        parameters: (dict, optional)
            A dictionary of model parameters. [Default: Empty]
        '''.format(", ".join(choices))

        
        self.spin_model = utils.remove_special_characters(spin_model.lower())
        self.number_of_samples = number_of_samples
        self.parameters = parameters
        self.parameters['minimum_primary_spin'] = parameters.get('minimum_primary_spin', 0) 
        self.parameters['maximum_primary_spin'] = parameters.get('maximum_primary_spin', 0.99)
        self.parameters['minimum_secondary_spin'] = parameters.get('minimum_secondary_spin', 0)
        self.parameters['maximum_secondary_spin'] = parameters.get('maximum_secondary_spin', 0.99)
        if 'mu_chi' in self.parameters and 'sigma_squared_chi' in self.parameters:
            self.parameters['alpha_chi'] = self.parameters['mu_chi'] * (self.parameters['mu_chi'] - self.parameters['mu_chi'] **2 
                                                                        - self.parameters['sigma_squared_chi']) / self.parameters['sigma_squared_chi']
            self.parameters['beta_chi'] = self.parameters['alpha_chi'] * (1. / self.parameters['mu_chi'] - 1.)
            
        if self.spin_model == 'default':
            self.spin_model = 'betagaussianuniformisotropic'
    def sample(self):
        '''
        Generate spin distribution samples based on the chosen parameterised model and its parameters.

        Returns:
        --------
            dict: A dictionary containing spin distribution samples.
        '''

        samples={}                
        err_msg = 'Please choose from one of these spin models {}'.format(choices)
        
        try:
            if 'nonspinning' in self.spin_model:
                logging.info('You chose non-spinning distribution')
                samples['chi_1'] = numpy.zeros(self.number_of_samples)
                samples['chi_2'] = numpy.zeros(self.number_of_samples)
        
            elif 'aligned' in self.spin_model:
                logging.info('Generating spin samples from Aligned Spin Distribution')
                if 'bilby' in self.spin_model or 'aligned' == self.spin_model:
                    chi_1 = bilby.gw.prior.AlignedSpin(name='chi_1', a_prior=bilby.gw.prior.Uniform(minimum=self.parameters['minimum_primary_spin'], 
                                                                                                    maximum=self.parameters['maximum_primary_spin']))
                    chi_2 = bilby.gw.prior.AlignedSpin(name='chi_2', a_prior=bilby.gw.prior.Uniform(minimum=self.parameters['minimum_secondary_spin'], 
                                                                                                    maximum=self.parameters['maximum_secondary_spin']))
                elif 'gaussianuniform' in self.spin_model:
                    mu_chi_1 = self.parameters['mu_chi_1']
                    sigma_chi_1 = self.parameters['sigma_chi_1']
                    chi_1 = bilby.core.prior.analytical.TruncatedGaussian(name='chi_1', 
                                                                          minimum=self.parameters['minimum_primary_spin'], 
                                                                          maximum=self.parameters['maximum_primary_spin'],
                                                                          mu=mu_chi_1, sigma=sigma_chi_1)
                    
                    chi_2 = bilby.gw.prior.Uniform(name='chi_2',
                                                   minimum=self.parameters['minimum_secondary_spin'], 
                                                   maximum=self.parameters['maximum_secondary_spin'])                    
                elif 'uniform' in self.spin_model: 
                    chi_1 = bilby.gw.prior.Uniform(name='chi_1', 
                                                   minimum=self.parameters['minimum_primary_spin'], 
                                                   maximum=self.parameters['maximum_primary_spin'])
                    chi_2 = bilby.gw.prior.Uniform(name='chi_2', 
                                                   minimum=self.parameters['minimum_secondary_spin'], 
                                                   maximum=self.parameters['maximum_secondary_spin'])
                elif 'beta' in self.spin_model:
                    chi_1 = bilby.gw.prior.AlignedSpin(name='chi_1', a_prior=bilby.core.prior.analytical.Beta(name='chi_1', 
                                                             alpha=self.parameters['alpha_chi'], 
                                                             beta=self.parameters['beta_chi'], 
                                                             minimum=self.parameters['minimum_primary_spin'], 
                                                             maximum=self.parameters['maximum_primary_spin'],))
                    chi_2 = bilby.gw.prior.AlignedSpin(name='chi_2', a_prior=bilby.core.prior.analytical.Beta(name='chi_2',
                                                             alpha=self.parameters['alpha_chi'], 
                                                             beta=self.parameters['beta_chi'],                                                      
                                                             minimum=self.parameters['minimum_secondary_spin'], 
                                                             maximum=self.parameters['maximum_secondary_spin']))

                samples['chi_1'] = chi_1.sample(self.number_of_samples)
                samples['chi_2'] = chi_2.sample(self.number_of_samples)
        
            elif 'isotropic' in self.spin_model:
                logging.info('Generating spin samples from Isotropic Spin Distribution')
                if 'bilby' in self.spin_model:                
                    a_1 = bilby.gw.prior.Uniform(minimum=self.parameters['minimum_primary_spin'], 
                                                 maximum=self.parameters['maximum_primary_spin'], name = 'a_1')
                    a_2 = bilby.gw.prior.Uniform(minimum=self.parameters['minimum_secondary_spin'], 
                                                 maximum=self.parameters['maximum_secondary_spin'], name = 'a_2')
                    
                    tilt_1 = bilby.core.prior.analytical.Sine(name='tilt_1')
                    tilt_2 = bilby.core.prior.analytical.Sine(name='tilt_2')            
                    samples['tilt_1'] = tilt_1.sample(self.number_of_samples)
                    samples['tilt_2'] = tilt_2.sample(self.number_of_samples)

                elif 'betagaussianuniform' in self.spin_model:
                    n = int(self.parameters['xi_spin'] * self.number_of_samples)
                    m = self.number_of_samples - n
                    a_1 = bilby.core.prior.analytical.Beta(name='a_1',   
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_primary_spin'], 
                                                           maximum=self.parameters['maximum_primary_spin'],)
                    
                    a_2 = bilby.core.prior.analytical.Beta(name='a_2',
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_secondary_spin'], 
                                                           maximum=self.parameters['maximum_secondary_spin'])
            
                    cos_tilt_gaussian = bilby.core.prior.analytical.TruncatedNormal(name='cos_tilt_gaussian', 
                                                                                    mu=1,
                                                                                    sigma=self.parameters['sigma_t'],
                                                                                    minimum=-1, maximum=1)
                    
                    cos_tilt_isotropic = bilby.core.prior.analytical.Uniform(name='cos_tilt_isotropic', 
                                                                             minimum=-1, 
                                                                             maximum=1)
                    
                    samples['tilt_1'] = numpy.arccos(numpy.concatenate([cos_tilt_gaussian.sample(n),
                                                                        cos_tilt_isotropic.sample(m)]))
                    samples['tilt_2'] = numpy.arccos(numpy.concatenate([cos_tilt_gaussian.sample(n),
                                                                        cos_tilt_isotropic.sample(m)]))
                    random.shuffle(samples['tilt_1'])
                    random.shuffle(samples['tilt_2'])
                elif 'betagaussian' in self.spin_model:
                    a_1 = bilby.core.prior.analytical.Beta(name='a_1',   
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_primary_spin'], 
                                                           maximum=self.parameters['maximum_primary_spin'],)
                    
                    a_2 = bilby.core.prior.analytical.Beta(name='a_2',
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_secondary_spin'], 
                                                           maximum=self.parameters['maximum_secondary_spin'])
            
                    cos_tilt = bilby.core.prior.analytical.TruncatedNormal(name='cos_tilt', 
                                                                           mu=1,
                                                                           sigma=self.parameters['sigma_t'],
                                                                           minimum=-1, 
                                                                           maximum=1)
                    
                    samples['tilt_1'] = numpy.arccos(cos_tilt.sample(self.number_of_samples))
                    samples['tilt_2'] = numpy.arccos(cos_tilt.sample(self.number_of_samples))                     
                                    
                elif 'beta' in self.spin_model:
                    a_1 = bilby.core.prior.analytical.Beta(name='a_1',   
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_primary_spin'], 
                                                           maximum=self.parameters['maximum_primary_spin'],)
                    a_2 = bilby.core.prior.analytical.Beta(name='a_2',
                                                           alpha=self.parameters['alpha_chi'], 
                                                           beta=self.parameters['beta_chi'], 
                                                           minimum=self.parameters['minimum_secondary_spin'], 
                                                           maximum=self.parameters['maximum_secondary_spin'])
                    tilt_1 = bilby.core.prior.analytical.Sine(name='tilt_1')
                    tilt_2 = bilby.core.prior.analytical.Sine(name='tilt_2')
                    samples['tilt_1'] = tilt_1.sample(self.number_of_samples)
                    samples['tilt_2'] = tilt_2.sample(self.number_of_samples)
                                                
                phi_12 = bilby.gw.prior.Uniform(name='phi_12', minimum=0, maximum=2 * numpy.pi, boundary='periodic')
                phi_jl = bilby.gw.prior.Uniform(name='phi_jl', minimum=0, maximum=2 * numpy.pi, boundary='periodic')
                
                samples['a_1'] = a_1.sample(self.number_of_samples)
                samples['a_2'] = a_2.sample(self.number_of_samples)
                samples['phi_12'] = phi_12.sample(self.number_of_samples)
                samples['phi_jl'] = phi_jl.sample(self.number_of_samples)
        except:
            raise ValueError(err_msg)
        return samples