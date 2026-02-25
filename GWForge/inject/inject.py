from ..ifo.detectors import Network
import bilby

bilby.core.utils.setup_logger(log_level="warning")


class Inject:
    def __init__(
        self,
        ifos,
        data,
        injection_parameters,
        injection_type="bbh",
        **waveform_arguments,
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

        self.waveform_arguments = {
            "waveform_approximant": waveform_arguments.get(
                "waveform_approximant", "IMRPhenomXPHM"
            ),
            "reference_frequency": waveform_arguments.get("reference_frequency", 7),
            "minimum_frequency": waveform_arguments.get("minimum_frequency", 7),
        }

    def inject(self):
        """
        Inject a gravitational wave signal into the interferometer network.
        Returns:
        --------------
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
        # Add snippet to check if the signal is in the data segment
        self.ifos.inject_signal(
            waveform_generator=waveform_generator, parameters=self.injection_parameters
        )
        return self.ifos
