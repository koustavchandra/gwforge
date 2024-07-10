from ..ifo.detectors import Network
import bilby, logging
from rich.progress import track
from ..conversion import *
from pycbc.detector import add_detector_on_earth, Detector
from ..utils import pycbc_labels
from pycbc.waveform import get_td_waveform
from pycbc.waveform import utils as wfutils
from pycbc.types import float64, float32
from gwpy.timeseries import TimeSeries
from lal import LIGOTimeGPS
import lalsimulation
injection_func_map = {
    numpy.dtype(float32): lambda *args: lalsimulation.SimAddInjectionREAL4TimeSeries(*args),
    numpy.dtype(float64): lambda *args: lalsimulation.SimAddInjectionREAL8TimeSeries(*args),
}
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S')

bilby.core.utils.setup_logger(log_level='warning')
class PyCBCInject:
    def __init__(self, 
                 ifos, 
                 data, 
                 injection_parameters,
                 waveform_arguments,
                 injection_type = 'bbh', ):
        '''
        Parameters:
        --------------
        ifos: list
            List of interferometers
        data: dict
            Dictionary of gwpy.timeseries.TimeSeries for each interferometer
        injection_parameters: dict
            Dictionary of injection parameters.
        injection_type: str, optional
            Type of injection (default: 'bbh').
        waveform_arguments: dict, optional
            Arguments for the waveform generator
        '''
        # Type checking
        if not isinstance(ifos, list):
            raise TypeError("ifos should be a list.")
        if not isinstance(data, dict):
            raise TypeError("data should be a dictionary.")
        if not isinstance(injection_parameters, dict):
            raise TypeError("parameters should be a dictionary.")
        self.ifos = Network(ifos=ifos).initialise_ifos()
            
        self.injection_parameters = injection_parameters
        self.injection_type = injection_type.lower()
        for ifo in self.ifos:
            ifo.strain_data.set_from_gwpy_timeseries(data[ifo.name])
            add_detector_on_earth(name=ifo.name,
                                  longitude=ifo.longitude_radians,
                                  latitude=ifo.latitude_radians,
                                  yangle=numpy.deg2rad(360 - ifo.xarm_azimuth), 
                                  xangle=numpy.deg2rad(360 - ifo.yarm_azimuth), 
                                  height=0,
                                  xlength=ifo.length * 1e3, 
                                  ylength=ifo.length * 1e3)            

        
        if self.injection_type in ['bbh', 'imbhb', 'pbh', 'imbbh', 'nsbh']:
            self.injection_type = 'bbh'
            self.frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole
            self.parameter_conversion = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters
        elif self.injection_type == 'bns':
            self.frequency_domain_source_model = bilby.gw.source.lal_binary_neutron_star
            self.parameter_conversion = bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters
        else:
            raise ValueError('Currently supports only CBC sources')
            
        self.waveform_arguments = waveform_arguments
        
    def get_pycbc_signal(self, ifo, signal_parameters):
        logging.info(f'Injecting signals with parameters: {signal_parameters}')
        hplus, hcross = get_td_waveform(**signal_parameters, 
                                        approximant= self.waveform_arguments['waveform_approximant'], 
                                        f_lower = self.waveform_arguments['minimum_frequency'], 
                                        f_ref = self.waveform_arguments['reference_frequency'],
                                        delta_t = 1./ifo.sampling_frequency)
        hplus.start_time = hplus.start_time + signal_parameters['tc']
        hcross.start_time = hcross.start_time + signal_parameters['tc']
        hplus = wfutils.taper_timeseries(hplus, 'startend')
        hcross = wfutils.taper_timeseries(hcross, 'startend')
        detector = Detector(ifo.name)
        signal = detector.project_wave(hplus, hcross, signal_parameters['ra'], signal_parameters['dec'], signal_parameters['psi'], method='lal', reference_time = signal_parameters['tc'])
        return signal
    
    def inject(self, ifo, strain, signal_parameters):
        signal = self.get_pycbc_signal(ifo=ifo, signal_parameters=signal_parameters)
        signal = signal.astype(strain.dtype)
        lal_signal = signal.lal()
        lalstrain = strain.lal()
        add_injection = injection_func_map[strain.dtype]
        add_injection(lalstrain, lal_signal, None)
        strain.data[:] = lalstrain.data.data[:]
        return strain
        
    def inject_signal_using_pycbc_method(self):
        '''
        Inject a gravitational wave signal into the ifo network using 
        pycbc methods
        Returns:
        -------
        '''
        for ifo in self.ifos:
            #convert it into pycbc time-series
            strain = ifo.strain_data.to_pycbc_timeseries()
            for k in track(range(len(self.injection_parameters['mass1'])), description='Injecting signals in {}'.format(ifo.name)):  
                signal_parameters = {key: self.injection_parameters[key][k] for key in self.injection_parameters.keys()}
                strain = self.inject(ifo=ifo,
                                     strain=strain, 
                                     signal_parameters=signal_parameters)
            strain = TimeSeries.from_pycbc(strain)
            ifo.strain_data.set_from_gwpy_timeseries(strain)
        return self.ifos
