from bilby.gw.conversion import *
import numpy
from lal import MTSUN_SI
from scipy.interpolate import interp1d


def generate_detector_frame_parameters(samples):
    """
    Generate detector frame parameters
    Parameters
    ----------
    sample: dict, pandas.DataFrame

    Returns
    -------
    output_sample: dict, pandas.DataFrame
    """
    output_samples = samples.copy()
    for key in [
        "mass_1_source",
        "mass_2_source",
        "chirp_mass_source",
        "total_mass_source",
    ]:
        if key in output_samples:
            detector_frame_key = key.replace("_source", "")
            output_samples[detector_frame_key] = output_samples[key] * (
                1 + output_samples["redshift"]
            )
    return output_samples


def get_lambda(file, source_mass):
    """
    Retrieve Lambda (tidal deformability) for a given source mass from a data file.

    Parameters
    ----------
    file : str
        The path to the data file containing mass and Lambda information.
    source_mass : float
        The mass for which to retrieve the corresponding Lambda.

    Returns
    -------
    float
        The Lambda value interpolated for the given source mass.
    """
    data = numpy.genfromtxt(file, names=True, delimiter=",")
    mass = numpy.zeros(len(data))
    lambda_data = numpy.zeros_like(mass)
    for j in range(len(data)):
        mass[j] = data[j][2]
        lambda_data[j] = (2.0 / 3.0) * data[j][5] / (data[j][0] ** 5.0)

    index_maximum_mass = numpy.argmax(mass)
    mass_new = mass[:index_maximum_mass]
    lambda_new = lambda_data[:index_maximum_mass]

    # Above the table's maximum mass the star has collapsed: no tides. Below
    # its minimum, hold the last tabulated value rather than extrapolate.
    mass_to_lambda = interp1d(
        mass_new, lambda_new, kind="cubic", bounds_error=False, fill_value=(lambda_new[0], 0.0)
    )

    return mass_to_lambda(source_mass)


def get_safe_signal_durations(
    mass_1,
    mass_2,
    spin_1z,
    spin_2z,
    waveform_minimum_frequency,
    safety=1.1,
    approximant="IMRPhenomXPHM",
):
    """
    Calculate the safe signal durations for gravitational waveforms.

    Parameters
    ----------
    mass_1: float or array-like
        The mass of the primary compact object(s) in solar masses.
    mass_2: float or array-like
        The mass of the secondary compact object(s) in solar masses.
    spin_1z: float or array-like
        The z-component of the spin of the primary compact object(s).
    spin_2z: (float or array-like)
        The z-component of the spin of the secondary compact object(s).
    waveform_minimum_frequency (float):
        The minimum frequency of the gravitational waveform.
    safety (float, optional):
        A safety factor to adjust the calculated durations. [Default: 1.1]
    approximant (str, optional):
        The waveform approximant to use. [Default: 'IMRPhenomXPHM']

    Returns
    -------
    float or array-like:
        The safe signal durations in seconds.

    Raises
    ------
    ValueError:
        If any of the input arrays (mass_1, mass_2, spin_1z, spin_2z) is not a float array.
    RuntimeError:
        If the provided waveform approximant is not supported.
    """
    import numpy
    import lalsimulation
    import lal

    # Convert input parameters to appropriate units
    mass_1, mass_2, spin_1z, spin_2z = (
        numpy.atleast_1d(numpy.asarray(mass_1, dtype=float)) * lal.MSUN_SI,
        numpy.atleast_1d(numpy.asarray(mass_2, dtype=float)) * lal.MSUN_SI,
        numpy.atleast_1d(numpy.asarray(spin_1z, dtype=float)),
        numpy.atleast_1d(numpy.asarray(spin_2z, dtype=float)),
    )

    waveform_minimum_frequency = float(waveform_minimum_frequency)

    # Calculate the safe signal durations based on the waveform approximant
    if "IMRPhenom" in approximant:
        durations = [
            safety
            * lalsimulation.SimIMRPhenomXASDuration(
                m1, m2, s1z, s2z, waveform_minimum_frequency
            )
            for m1, m2, s1z, s2z in zip(mass_1, mass_2, spin_1z, spin_2z)
        ]
    elif "SEOBNR" in approximant:
        durations = [
            safety
            * lalsimulation.SimIMRSEOBNRv5ROMTimeOfFrequency(
                waveform_minimum_frequency, m1, m2, s1z, s2z
            )
            for m1, m2, s1z, s2z in zip(mass_1, mass_2, spin_1z, spin_2z)
        ]
    else:
        raise RuntimeError(
            "Failed to compute durations for approximant {}".format(approximant)
        )

    # Return the value directly if durations is a single-element array
    return durations[0] if len(durations) == 1 else durations


def time_to_coalescence(frequency, chirp_mass, symmetric_mass_ratio, mode=2):
    r"""Time remaining to merger at a given frequency, to 3.5PN.

    This is the time-frequency relation that turns the frequency-domain response
    into a time-dependent one: the bin at frequency :math:`f` is emitted
    :math:`\tau(f)` seconds before merger, so the Earth's orientation there is
    the orientation at ``geocent_time - tau(f)``.

    Uses eq. (3.8b) of `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_.

    Parameters
    ----------
    frequency : array_like
        Frequency in Hz.
    chirp_mass : float
        Detector-frame chirp mass in solar masses.
    symmetric_mass_ratio : float
        :math:`\eta = m_1 m_2 / (m_1+m_2)^2`.
    mode : int
        Harmonic mode number. Mode ``m`` reaches a given frequency when the
        dominant quadrupole is at :math:`2f/m`, so ``tau_m(f) = tau_2(2f/m)``.

    Returns
    -------
    numpy.ndarray
        Seconds before coalescence, positive and decreasing with frequency.
    """
    frequency = numpy.asarray(frequency, dtype=float)
    if mode != 2:
        frequency = 2.0 * frequency / mode

    total_mass = chirp_mass * MTSUN_SI / (symmetric_mass_ratio ** (3.0 / 5.0))
    # Guard f = 0 (the DC bin): tau diverges there, and the bin is masked out of
    # every analysis anyway.
    with numpy.errstate(divide="ignore", invalid="ignore"):
        velocity = (numpy.pi * total_mass * frequency) ** (1.0 / 3.0)

        eta = symmetric_mass_ratio
        eta2 = eta * eta
        prefactor = 5.0 / 256 * total_mass / (eta * velocity**8)

        through_2p5 = (
            1.0
            + (743.0 / 252.0 + 11.0 / 3.0 * eta) * velocity**2
            - 32.0 / 5.0 * numpy.pi * velocity**3
            + (3058673.0 / 508032.0 + 5429.0 / 504.0 * eta + 617.0 / 72.0 * eta2)
            * velocity**4
            - (7729.0 / 252.0 - 13.0 / 3.0 * eta) * numpy.pi * velocity**5
        )
        order_3 = (
            -10052469856691.0 / 23471078400.0
            + 128.0 / 3.0 * numpy.pi**2
            + 6848.0 / 105.0 * numpy.euler_gamma
            + (3147553127.0 / 3048192.0 - 451.0 / 12.0 * numpy.pi**2) * eta
            - 15211.0 / 1728.0 * eta2
            + 25565.0 / 1296.0 * eta2 * eta
            + 3424.0 / 105.0 * numpy.log(16.0 * velocity**2)
        ) * velocity**6
        order_3p5 = (
            -15419335.0 / 127008.0 - 75703.0 / 756.0 * eta + 14809.0 / 378.0 * eta2
        ) * numpy.pi * velocity**7

        tau = prefactor * (through_2p5 + order_3 + order_3p5)

    return numpy.nan_to_num(tau, nan=0.0, posinf=0.0, neginf=0.0)


def chirp_mass_and_symmetric_mass_ratio(parameters):
    """Extract ``(chirp_mass, symmetric_mass_ratio)`` from a parameter dict.

    Accepts any of the mass parameterisations GWForge and bilby pass around:
    component masses, chirp mass with mass ratio, or chirp mass with symmetric
    mass ratio.

    Returns
    -------
    tuple of float, or None
        ``None`` when the masses cannot be determined -- e.g. a custom non-CBC
        source model -- in which case the caller falls back to a static response.
    """
    if "chirp_mass" in parameters and "symmetric_mass_ratio" in parameters:
        return float(parameters["chirp_mass"]), float(
            parameters["symmetric_mass_ratio"]
        )
    if "mass_1" in parameters and "mass_2" in parameters:
        mass_1 = float(parameters["mass_1"])
        mass_2 = float(parameters["mass_2"])
        total = mass_1 + mass_2
        eta = mass_1 * mass_2 / total**2
        return total * eta ** (3.0 / 5.0), eta
    if "chirp_mass" in parameters and "mass_ratio" in parameters:
        mass_ratio = float(parameters["mass_ratio"])
        eta = mass_ratio / (1.0 + mass_ratio) ** 2
        return float(parameters["chirp_mass"]), eta
    return None
