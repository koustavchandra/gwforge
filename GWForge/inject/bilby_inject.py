from ..ifo.detectors import Network
from ..ifo.antenna import inject_signal_with_response
import bilby
from ..conversion import *

bilby.core.utils.setup_logger(log_level="warning")


class BilbyInject:
    def __init__(
        self,
        ifos,
        data,
        injection_parameters,
        waveform_arguments,
        injection_type="bbh",
        earth_rotation=True,
        finite_size=True,
    ):
        """
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
        earth_rotation: bool, optional
            Evolve the antenna pattern and arrival delay along the signal's
            time-frequency track rather than freezing them at coalescence
            (default: True).
        finite_size: bool, optional
            Apply the per-arm finite light-travel-time transfer function
            (default: True).
        waveform_arguments: dict, optional
            Arguments for the waveform generator
        """
        # Type checking
        if not isinstance(ifos, list):
            raise TypeError("ifos should be a list.")
        if not isinstance(data, dict):
            raise TypeError("data should be a dictionary.")
        if not isinstance(injection_parameters, dict):
            raise TypeError("parameters should be a dictionary.")
        self.ifos = Network(ifos=ifos).initialise_ifos()

        for ifo in self.ifos:
            ifo.strain_data.set_from_gwpy_timeseries(data[ifo.name])

        self.injection_parameters = injection_parameters
        self.injection_type = injection_type.lower()

        if self.injection_type in ["bbh", "imbhb", "pbh", "imbbh", "nsbh"]:
            self.injection_type = "bbh"
            self.frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole
            self.parameter_conversion = (
                bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters
            )
        elif self.injection_type == "bns":
            self.frequency_domain_source_model = bilby.gw.source.lal_binary_neutron_star
            self.parameter_conversion = (
                bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters
            )
        else:
            raise ValueError("Currently supports only CBC sources")

        self.waveform_arguments = waveform_arguments
        self.earth_rotation = earth_rotation
        self.finite_size = finite_size

    def inject_signal_using_bilby_method(self):
        """
        Inject a gravitational wave signal into the interferometer network.
        using bilby methods
        Returns:
        -------
        Network: bilby.gw.detectors.InterferometerList
            Injected interferometer network
        """
        waveform_generator = bilby.gw.WaveformGenerator(
            duration=self.ifos.duration,
            sampling_frequency=self.ifos.sampling_frequency,
            frequency_domain_source_model=self.frequency_domain_source_model,
            parameter_conversion=self.parameter_conversion,
            waveform_arguments=self.waveform_arguments,
        )
        # bilby's InterferometerList.inject_signal hard-codes its own
        # long-wavelength, static-pattern response, so the projection is done
        # here instead. The meta_data it fills is identical.
        inject_signal_with_response(
            self.ifos,
            waveform_generator=waveform_generator,
            parameters=self.injection_parameters,
            earth_rotation=self.earth_rotation,
            finite_size=self.finite_size,
        )
        return self.ifos
