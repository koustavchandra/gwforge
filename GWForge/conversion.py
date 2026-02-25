from bilby.gw.conversion import *
import numpy
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
            output_samples[detector_frame_key] = output_samples[key] * (1 + output_samples["redshift"])
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

    mass_to_lambda = interp1d(mass_new, lambda_new, kind="cubic")

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

    Parameters:
    -----------
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

    Returns:
    --------
    float or array-like:
        The safe signal durations in seconds.

    Raises:
    -------
    ValueError:
        If any of the input arrays (mass_1, mass_2, spin_1z, spin_2z) is not a float array.
    RuntimeError:
        If the provided waveform approximant is not supported.
    """
    import numpy
    import lalsimulation
    import lal

    # Ensure that input parameters are float arrays
    for param_name, param in [
        ("mass_1", mass_1),
        ("mass_2", mass_2),
        ("spin_1z", spin_1z),
        ("spin_2z", spin_2z),
    ]:
        if not isinstance(param, (float, numpy.ndarray, list)) or not numpy.issubdtype(numpy.asarray(param).dtype, numpy.floating):
            raise ValueError(f"{param_name} must be a float or a float array")

    # Convert input parameters to appropriate units
    mass_1, mass_2, spin_1z, spin_2z = (
        numpy.atleast_1d(numpy.asarray(mass_1)) * lal.MSUN_SI,
        numpy.atleast_1d(numpy.asarray(mass_2)) * lal.MSUN_SI,
        numpy.atleast_1d(numpy.asarray(spin_1z)),
        numpy.atleast_1d(numpy.asarray(spin_2z)),
    )

    waveform_minimum_frequency = float(waveform_minimum_frequency)

    # Calculate the safe signal durations based on the waveform approximant
    if "IMRPhenom" in approximant:
        durations = [safety * lalsimulation.SimIMRPhenomXASDuration(m1, m2, s1z, s2z, waveform_minimum_frequency) for m1, m2, s1z, s2z in zip(mass_1, mass_2, spin_1z, spin_2z)]
    elif "SEOBNR" in approximant:
        durations = [safety * lalsimulation.SimIMRSEOBNRv5ROMTimeOfFrequency(m1, m2, s1z, s2z, waveform_minimum_frequency) for m1, m2, s1z, s2z in zip(mass_1, mass_2, spin_1z, spin_2z)]
    else:
        raise RuntimeError("Failed to compute durations for approximant {}".format(approximant))

    # Return the value directly if durations is a single-element array
    return durations[0] if len(durations) == 1 else durations
