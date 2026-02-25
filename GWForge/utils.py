import numpy
import pylab
import h5py
import logging
import os
import fnmatch
import bilby
import seaborn as sns

sns.set_context("talk")
sns.set(font_scale=1.7)
sns.set_palette("colorblind")
sns.set_style("ticks")

pylab.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "stixgeneral",
        "mathtext.fontset": "stix",
    }
)

pylab.rcParams["axes.linewidth"] = 1

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)


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


def cornerplot(file, parameters=None, save=None):
    """
    Create a corner plot from samples stored in an HDF5 file.

    Parameters
    ----------
    file : str
        The path to the HDF5 file containing samples.
    parameters : list of str, optional
        A list of parameters for which to create the corner plot. If not provided,
        default parameters ['mass_1_source', 'mass_2_source', 'chi_eff', 'chi_p', 'theta_jn', 'redshift']
        will be used.
    save : str, optional
        The path to save the generated corner plot. If not provided, the plot will be displayed.
    """
    logging.info("Making Corner Plot")
    GWlatex_labels = {
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

    from corner import corner

    samples = {}
    with h5py.File(file, "r") as f:
        for key in f.keys():
            samples[key] = f[key][:]

    labels = []
    data = []
    if parameters:
        for parameter in parameters:
            data.append(samples[parameter])
            labels.append(GWlatex_labels[parameter])
    else:
        parameters = [
            "mass_1_source",
            "mass_2_source",
            "chi_eff",
            "chi_p",
            "theta_jn",
            "redshift",
        ]
        for parameter in parameters:
            data.append(samples[parameter])
            labels.append(GWlatex_labels[parameter])

    figure = pylab.figure(figsize=(3 * len(parameters), 3 * len(parameters)))
    defaults_kwargs = dict(
        bins=50,
        smooth=0.9,
        title_kwargs=dict(fontsize=16),
        color="#ca0020",
        fig=figure,
        quantiles=[0.16, 0.84],
        levels=(1 - numpy.exp(-0.5), 1 - numpy.exp(-2), 1 - numpy.exp(-9 / 2.0)),
        plot_density=False,
        plot_datapoints=True,
        fill_contours=True,
        max_n_ticks=3,
    )

    corner(numpy.asarray(data).T, labels=labels, **defaults_kwargs)
    if save:
        figure.savefig(save, bbox_inches="tight", dpi=100)
    logging.info("Done!")


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

            # Check if the file overlaps with the specified time range
            if (
                start_time is None
                or (
                    gps_start_time <= end_time
                    and gps_start_time + duration > start_time
                )
            ) and (
                end_time is None
                or (
                    gps_start_time < end_time
                    and gps_start_time + duration >= start_time
                )
            ):
                filepaths.append(file_path)
                filenames.append(filename)

    return filenames, filepaths


def filter_times_by_frame_files(times, frame_files):
    """
    Filter a list of injections based on frame file availability

    Parameters:
    -----------
    times: list
        List of tc
    frame_files:
        List of frame files
    """
    filtered_times = []

    for frame_file in frame_files:
        parts = frame_file.split("-")
        gps_start_time = int(parts[1])
        duration = int(parts[-1].split(".")[0])
        gps_end_time = gps_start_time + duration

        for time in times:
            if gps_start_time <= time <= gps_end_time:
                filtered_times.append(time)

    return filtered_times


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
        save_data.write(
            target=os.path.join(
                ifo_directory, f"{ifo.name}-{int(start_time)}-{int(dur)}.h5"
            ),
            overwrite=True,
        )
        # Delete the cropped data to free up memory
        del save_data
        # Update the start time for the next iteration
        start_time = end


def custom_optionxform(option):
    # Replace hyphens with underscores
    return option.replace("_", "-")


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
    "dec": bilby.core.prior.analytical.Uniform(
        name="dec", minimum=0, maximum=numpy.pi, boundary="periodic"
    ),
    "theta_jn": bilby.core.prior.analytical.Uniform(
        name="theta_jn", minimum=0, maximum=numpy.pi, boundary="periodic"
    ),
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
