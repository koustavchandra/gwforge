import numpy
import logging
import os
import fnmatch
import bilby


def remove_special_characters(
    input_string, characters_to_remove=["+", "-", "_", " ", "#"]
):
    """
    Remove specified special characters from a given input string.

    Parameters
    ----------
    input_string : str
        The input string from which to remove special characters.
    characters_to_remove : list of str, optional
        A list of special characters to remove from the input string. Defaults to ["+", "-", "_", " "].

    Returns
    -------
    str
        The input string with specified special characters removed.
    """
    result_string = "".join(
        char for char in input_string if char not in characters_to_remove
    )
    return result_string


def hdf_append(f, key, value):
    """
    Append a value to an HDF5 dataset or create a new dataset if the key does not exist.

    Parameters:
    ----------
    f : (h5py.File)
        An HDF5 file object.
    key : (str)
        The key to identify the dataset within the HDF5 file.
    value : (float or numpy.ndarray)
        The value to be appended to the dataset.

    If the dataset with the specified key already exists, the function appends the given value to the existing dataset.
    If the dataset does not exist, a new dataset is created with the specified key, and the value is stored.

    Note: The function ensures that the stored value is a 1-dimensional array.
    """
    if key in f:
        value = numpy.atleast_1d(value)
        tmp = numpy.concatenate([f[key][:], value])
        del f[key]
        f[key] = tmp
    else:
        # Convert scalar value to 1-dimensional array
        f[key] = numpy.atleast_1d(value)


# LaTeX axis labels for every parameter GWForge passes around. Module level so
# the Fisher plotting in :mod:`GWForge.fisher.plot` labels its axes the same way
# :func:`cornerplot` does rather than keeping a second copy.
GWLATEX_LABELS = {
    "luminosity_distance": r"$d_{L} [\mathrm{Mpc}]$",
    "geocent_time": r"$t_{c} [\mathrm{s}]$",
    "dec": r"$\delta [\mathrm{rad}]$",
    "ra": r"$\alpha [\mathrm{rad}]$",
    "a_1": r"$a_{1}$",
    "a_2": r"$a_{2}$",
    "phi_jl": r"$\phi_{JL} [\mathrm{rad}]$",
    "phase": r"$\phi [\mathrm{rad}]$",
    "psi": r"$\Psi [\mathrm{rad}]$",
    "iota": r"$\iota [\mathrm{rad}]$",
    "tilt_1": r"$\theta_{1} [\mathrm{rad}]$",
    "tilt_2": r"$\theta_{2} [\mathrm{rad}]$",
    "phi_12": r"$\phi_{12} [\mathrm{rad}]$",
    "mass_2": r"$m_{2} [M_{\odot}]$",
    "mass_1": r"$m_{1} [M_{\odot}]$",
    "total_mass": r"$M [M_{\odot}]$",
    "chirp_mass": r"$\mathcal{M} [M_{\odot}]$",
    "spin_1x": r"$S_{1x}$",
    "spin_1y": r"$S_{1y}$",
    "spin_1z": r"$S_{1z}$",
    "spin_2x": r"$S_{2x}$",
    "spin_2y": r"$S_{2y}$",
    "spin_2z": r"$S_{2z}$",
    "chi_1": r"$\chi_{1}$",
    "chi_2": r"$\chi_{2}$",
    "chi_p": r"$\chi_{\mathrm{p}}$",
    "chi_eff": r"$\chi_{\mathrm{eff}}$",
    "mass_ratio": r"$q$",
    "symmetric_mass_ratio": r"$\eta$",
    "inverted_mass_ratio": r"$1/q$",
    "cos_tilt_1": r"$\cos{\theta_{1}}$",
    "cos_tilt_2": r"$\cos{\theta_{2}}$",
    "redshift": r"$z$",
    "mass_1_source": r"$m_{1}^{\mathrm{source}} [M_{\odot}]$",
    "mass_2_source": r"$m_{2}^{\mathrm{source}} [M_{\odot}]$",
    "chirp_mass_source": r"$\mathcal{M}^{\mathrm{source}} [M_{\odot}]$",
    "total_mass_source": r"$M^{\mathrm{source}} [M_{\odot}]$",
    "cos_iota": r"$\cos{\iota}$",
    "theta_jn": r"$\theta_{JN} [\mathrm{rad}]$",
    "cos_theta_jn": r"$\cos{\theta_{JN}}$",
    "lambda_1": r"$\lambda_{1}$",
    "lambda_2": r"$\lambda_{2}$",
    "lambda_tilde": r"$\tilde{\lambda}$",
    "delta_lambda": r"$\delta\lambda$",
}


def split_duration(duration, size=4096.0):
    """
    Split a duration into chunks of a specified size.

    Parameters:
    -----------
    - duration: int
        The total duration to split
    - size: int
        The size of each chunk [Default:4096]

    Returns:
    list: A list containing chunks of the specified size.
    """
    result = []
    while duration > size:
        result.append(size)
        duration -= size
    result.append(duration)
    return result


def find_frame_files(directory, filePattern="*gwf", start_time=None, end_time=None):
    """
    Get a list of frame files within a given directory optionally filtered by time range
    Parameters:
    -----------
    directory : str
        Path to frames directory
    filepattern : str
        A Unix shell-style wildcard pattern to match files. [Default *gwf]
    start_time : float
        GPS start time for time range filter
    end_time : float
        GPS end time for time range filter
    """

    filenames, filepaths = [], []
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            file_path = os.path.join(path, filename)

            # Extract GPS start time and duration from filename
            parts = filename.split("-")
            gps_start_time = float(parts[1])
            duration = float(parts[-1].split(".")[0])

            # Check if the file overlaps with the specified time range.
            # Each bound is applied independently so that a one-sided range
            # (only start_time or only end_time) does not raise on a None
            # comparison.
            gps_end_time = gps_start_time + duration
            after_start = start_time is None or gps_end_time > start_time
            before_end = end_time is None or gps_start_time < end_time
            if after_start and before_end:
                filepaths.append(file_path)
                filenames.append(filename)

    return filenames, filepaths


def source_type(name):
    """Canonical source-type token: ``bhns`` for a neutron-star--black-hole binary.

    ``nsbh`` is accepted everywhere as an alias; everything else is lowercased.
    """
    name = name.lower()
    return "bhns" if name == "nsbh" else name


def split_odd_even(items):
    """
    Split a list into its odd-indexed and even-indexed sublists.

    Used by the workflow generator to schedule injection jobs for adjacent data
    segments in two non-overlapping passes: because neighbouring segments share
    frame files (window overlap), the even-indexed segments are made to depend on
    the odd-indexed ones so that no two adjacent segments are written concurrently.

    Parameters:
    -----------
    items: list
        The list to split.

    Returns:
    --------
    tuple(list, list): (odd-indexed items, even-indexed items)

    Example:
    --------
    >>> split_odd_even(['a', 'b', 'c', 'd', 'e'])
    (['b', 'd'], ['a', 'c', 'e'])
    """
    return list(items[1::2]), list(items[0::2])


def generate_frame_file_sublists(frame_files, window_size=3):
    """
    Generate sublists of frame files with a specified window size.

    Parameters:
    -----------
    frame_files: list
        A list of file paths representing frame files.
    window_size: int, optional
        The size of the sliding window to create sublists. [Defaults: 3]

    Returns:
    --------
    list: A list of sublists containing frame file paths.

    Example:
    --------
    >>> frame_files = ['file1.gwf', 'file2.gwf', 'file3.gwf', 'file4.gwf']
    >>> generate_frame_file_sublists(frame_files, window_size=2)
    [['file1.gwf', 'file2.gwf'], ['file2.gwf', 'file3.gwf'], ['file3.gwf', 'file4.gwf']]
    """
    sublists = []
    for i in range(len(frame_files) - window_size + 1):
        sublist = frame_files[i : i + window_size]
        sublists.append(sublist)
    return sublists


def update_ET_channels(channel_dict):
    updated_channels = {}
    for ifo, channel in channel_dict.items():
        if channel.startswith(ifo + ":") and ifo == "ET":
            # Extract the suffix after 'ifo:'
            suffix = channel[len(ifo) + 1 :]
            # Generate updated channels for ET1, ET2, ET3
            updated_channels.update(
                {"ET1": "ET1:" + suffix, "ET2": "ET2:" + suffix, "ET3": "ET3:" + suffix}
            )
        else:
            updated_channels[ifo] = channel
    return updated_channels


# Function to save frame files
def save_frame_files(ifo, start_time, duration, ifo_directory):
    from gwpy.timeseries import TimeSeries

    logging.info(f"Saving {ifo.name} frame files with injections")
    # Create a TimeSeries object for the interferometer data
    data = TimeSeries(
        data=ifo.time_domain_strain,
        times=ifo.time_array,
        name=f"{ifo.name}:INJ",
        channel=f"{ifo.name}:INJ",
    )
    # Iterate through the provided durations
    for dur in duration:
        end = start_time + dur
        # Crop the data based on the start and end times
        save_data = data.crop(start=start_time, end=end)
        # Write the cropped data to a frame file
        # The path must be positional and the format explicit: gwpy >=4 infers
        # the format from *args, so passing it as target= makes identification
        # fail with an IndexError.
        save_data.write(
            os.path.join(ifo_directory, f"{ifo.name}-{int(start_time)}-{int(dur)}.h5"),
            format="hdf5",
            overwrite=True,
        )
        # Delete the cropped data to free up memory
        del save_data
        # Update the start time for the next iteration
        start_time = end


pycbc_labels = {
    "mass_1": "mass1",
    "mass_2": "mass2",
    "spin_1x": "spin1x",
    "spin_1y": "spin1y",
    "spin_1z": "spin1z",
    "spin_2x": "spin2x",
    "spin_2y": "spin2y",
    "spin_2z": "spin2z",
    "lambda_1": "lambda1",
    "lambda_2": "lambda2",
    "geocent_time": "tc",
    "ra": "ra",
    "dec": "dec",
    "psi": "psi",
    "theta_jn": "inclination",
    "phase": "coa_phase",
    "luminosity_distance": "distance",
}

# TODO: Remove this dependency in future versions.
reference_prior_dict = {
    "ra": bilby.core.prior.analytical.Uniform(
        name="ra", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
    ),
    # dec is cosine-distributed on [-pi/2, pi/2] for an isotropic sky, not uniform on [0, pi]
    "dec": bilby.core.prior.analytical.Cosine(name="dec"),
    # theta_jn is sine-distributed on [0, pi] for isotropic orientations, not uniform
    "theta_jn": bilby.core.prior.analytical.Sine(name="theta_jn"),
    "psi": bilby.core.prior.analytical.Uniform(
        name="psi", minimum=0, maximum=numpy.pi, boundary="periodic"
    ),
    "luminosity_distance": bilby.gw.prior.UniformSourceFrame(
        name="luminosity_distance", minimum=10, maximum=1000
    ),
    "phase": bilby.core.prior.analytical.Uniform(
        name="phase", minimum=0, maximum=2 * numpy.pi, boundary="periodic"
    ),
}
