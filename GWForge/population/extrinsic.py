import numpy, bilby, logging

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S')


class Extrinsic:
    def __init__(self, 
                 number_of_samples, 
                 prior_file=None, 
                 inclination_distribution=None):
        '''
        Parameters
        ----------
        number_of_samples : (int)
            The number of samples to generate. [Ideal: Exactly same as redshift samples]
        prior_file : str
            Bilby prior file with custom definition
        '''
        self.number_of_samples = number_of_samples
        self.prior_file = prior_file
        self.inclination_distribution = inclination_distribution
    def schutz_inclination_prob(self, inclination):
        '''
        See Eq. (25) of https://arxiv.org/pdf/1102.5421.pdf
        '''
        return 1/8 * (1 +  6 * numpy.cos(inclination) ** 2 + numpy.cos(inclination) ** 4)
    def sample(self):
        samples = {}
        if self.prior_file:
            prior = bilby.gw.prior.PriorDict(filename=self.prior_file)
        else:
            logging.warn('Using default priors')
            prior = bilby.gw.prior.PriorDict()
            prior['dec'] = bilby.core.prior.analytical.Cosine(name='dec')
            prior['ra'] = bilby.core.prior.analytical.Uniform(name='ra', minimum=0, maximum=2 * numpy.pi, boundary='periodic')
            prior['theta_jn'] = bilby.core.prior.analytical.Sine(name='theta_jn')
            prior['psi'] = bilby.core.prior.analytical.Uniform(name='psi', minimum=0, maximum=numpy.pi, boundary='periodic')
            prior['phase'] = bilby.core.prior.analytical.Uniform(name='phase', minimum=0, maximum=2 * numpy.pi, boundary='periodic')
            if  self.inclination_distribution is not None and 'schutz' in self.inclination_distribution.lower():
                inclination = numpy.linspace(0, numpy.pi, 5001)
                prob = self.schutz_inclination_prob(inclination)
                prior['theta_jn'] = bilby.core.prior.Interped(inclination, prob, minimum=numpy.min(inclination), maximum=numpy.max(inclination))
            elif self.inclination_distribution is not None and 'schutz' not in self.inclination_distribution.lower():
                raise ValueError('Distribution is not implemented!')
        priors = prior.sample(self.number_of_samples)
        for key in ['ra', 'dec', 'theta_jn', 'psi', 'phase']:
            samples[key] = priors[key]
        return samples