import bilby
import numpy
import logging
import os
from bilby.gw.detector.psd import PowerSpectralDensity as psd
from bilby.gw.detector.networks import TriangularInterferometer

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)


def load_xg_interferometer(filename, power_spectral_density):
    """Load a GWForge XG ``.ifo`` file, attaching ``power_spectral_density``.

    This mirrors ``bilby.gw.detector.load_interferometer``'s tiny key = value
    parser but injects a valid PSD before the interferometer is constructed.
    Two things make this necessary rather than a plain ``load_interferometer``:

    - bilby's ``TriangularInterferometer`` indexes ``power_spectral_density[ii]``
      for each of the three arms, so a ``.ifo`` file carrying
      ``power_spectral_density = None`` (as ET does) raises a ``TypeError`` at
      load time.
    - even for a single detector the PSD has to reach the interferometer(s), not
      just the enclosing list, so we set it explicitly on every arm.

    Parameters
    ----------
    filename : str
        Path to the ``.ifo`` file.
    power_spectral_density : bilby.gw.detector.psd.PowerSpectralDensity
        The PSD to attach to the interferometer (all three arms for a triangle).

    Returns
    -------
    bilby.gw.detector.Interferometer or bilby.gw.detector.networks.TriangularInterferometer
    """
    parameters = {}
    with open(filename, "r") as parameter_file:
        for line in parameter_file:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            key, _, value = line.partition("=")
            parameters[key.strip()] = eval(value)
    parameters["power_spectral_density"] = power_spectral_density
    shape = str(parameters.pop("shape", "L"))
    if shape.lower() in ("triangular", "triangle"):
        ifo = TriangularInterferometer(**parameters)
        for arm in ifo:
            arm.power_spectral_density = power_spectral_density
    else:
        ifo = bilby.gw.detector.Interferometer(**parameters)
        ifo.power_spectral_density = power_spectral_density
    return ifo


class IFO:
    """
    Usage:
    -----
    >>> ifo = IFO(name='H1', asharp=True)
    >>> ifo.initialise_ifo(seed=42)

    >>> ifo = IFO(name='ET')
    >>> ifo.initialize_ifo()
    """

    def __init__(self, name, psd_file=None, asharp=None):
        """
        Parameters:
        ----------
        name : str
           Name of interferometer
        psd_file : str (optional)
           Path of PSD/ASD File [Default: None]
        asharp : bool (optional)
           A flag indicating whether to use A# PSD for H1, L1 and A1 ifo
        """
        self.name = name.upper()
        self.psd_file = psd_file
        self.asharp = asharp
        if self.name == "CEA":
            self.name = "CE40"
        elif self.name == "CEB":
            self.name = "CE20"

    def initialise_ifo(self, seed=None):
        """
        Initialise the interferometer with optional random seed.

        Parameters:
        -----------
        seed : int (optional)
           The random seed for initialisation. [Default: None]

        Returns:
        --------
          ifo: bilby.gw.detector.Intereferometer

        Raises:
        -------
          ValueError: If the interferometer is not implemented.
        """
        if seed is not None:
            numpy.random.seed(seed)
        try:
            # Execute if ifo is ET or CEs
            filename = os.path.join(
                os.path.dirname(__file__), "ifos", "{}.ifo".format(self.name)
            )
            if not os.path.isfile(filename):
                raise FileNotFoundError(filename)
            if self.name == "ET":
                file_type = "psd"
                load_psd = psd.from_power_spectral_density_file
            else:
                file_type = "asd"
                load_psd = psd.from_amplitude_spectral_density_file
            noise_file = os.path.join(
                os.path.dirname(__file__),
                "noise_curves",
                "{}-{}.txt".format(self.name, file_type),
            )
            temp = load_psd(noise_file)
            # Attach the PSD as the interferometer is built: bilby's
            # TriangularInterferometer rejects a None PSD (as ET's .ifo carries),
            # and the PSD must reach each ET arm, not just the enclosing list.
            ifo = load_xg_interferometer(filename, temp)
        except FileNotFoundError:
            # Execute if detector is LIGO, Virgo, NEMO or LISA or any other bilby ifos
            logging.info(
                "{} is not XG. Checking whether it is implemented in Bilby".format(
                    self.name
                )
            )
            ifo = bilby.gw.detector.get_empty_interferometer(self.name)
            if self.name in ["H1", "L1", "A1"] and self.asharp:
                logging.info("Assigning {} ifo A# PSD".format(self.name))
                asharp_file = os.path.join(
                    os.path.dirname(__file__), "noise_curves", "Asharp-asd.txt"
                )
                ifo.power_spectral_density = psd.from_amplitude_spectral_density_file(
                    asharp_file
                )
        except Exception as e:
            raise ValueError(
                "Interferometer {} is not implemented. Error: {}".format(self.name, e)
            )
        if self.psd_file:
            temp = numpy.loadtxt(self.psd_file)
            frequency_array, psd_array = temp[:, 0], temp[:, -1]
            if numpy.min(psd_array) > 1e-30:
                load_psd = psd.from_amplitude_spectral_density_array
            else:
                load_psd = psd.from_power_spectral_density_array
            ifo.power_spectral_density = load_psd(frequency_array, psd_array)
        return ifo

    # Define an alias for the function
    initialize_ifo = initialise_ifo


class Network(IFO):
    """
    Usage:
    ------
    >>> ifos = Network(ifos = ['ET', 'CE20', 'CE40'])
    >>> ifos.initialise_network()
    """

    def __init__(self, ifos, psd_files=None):
        """
        Parameters:
        ----------
        ifos: str list
           A list of interferometers to include in the network
        psd_files: dict or list or str, optional
           PSD files corresponding to interferometers [Default: None]
        """
        ifo_list = []
        k = 0
        for ifo in ifos:  # Fix the typo here
            if isinstance(ifo, str) and not psd_files:
                ifo_list.append(IFO(name=ifo).initialise_ifo())
            elif isinstance(ifo, str) and isinstance(psd_files, dict):
                ifo_list.append(IFO(name=ifo, psd_file=psd_files[ifo]).initialise_ifo())
            elif isinstance(ifo, str) and isinstance(psd_files, list):
                ifo_list.append(IFO(name=ifo, psd_file=psd_files[k]).initialise_ifo())
                k = k + 1
            elif (
                isinstance(ifo, str)
                and isinstance(psd_files, str)
                and psd_files.lower() == "asharp"
            ):
                ifo_list.append(IFO(name=ifo, asharp=True).initialise_ifo())
            else:
                raise ValueError("Input is not compatible")
        self.ifos = bilby.gw.detector.InterferometerList(ifo_list)

    def initialise_network(self):
        """
        Returns:
        --------
        ifos: bilby.gw.detector.InterferometerList
        """
        return self.ifos

    # Define aliases for the function
    initialize_network = initialise_network
    initialise_ifos = initialise_network
    initialize_ifos = initialise_network
