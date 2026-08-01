"""Configuration handling for ``gwforge_fisher``.

Kept out of the executable so the parsing is importable and testable, and so the
CLI stays a thin driver like the other ``bin/`` scripts.
"""

import importlib
import json
import logging

import bilby
import h5py
import numpy

from ..conversion import get_safe_signal_durations
from ..ifo import antenna
from .parameters import DEFAULT_PARAMETERS
from .stepsize import TARGET_FRACTIONAL_CHANGE


def _optional_source_model(module_name, attribute):
    """Import a source model from an optional package, or None if absent."""
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        return None
    return getattr(module, attribute, None)


SOURCE_MODELS = {
    "lal_binary_black_hole": bilby.gw.source.lal_binary_black_hole,
    "lal_binary_neutron_star": bilby.gw.source.lal_binary_neutron_star,
    "gwsignal_binary_black_hole": getattr(
        bilby.gw.source, "gwsignal_binary_black_hole", None
    ),
    "lal_eccentric_binary_black_hole_no_spins": bilby.gw.source.lal_eccentric_binary_black_hole_no_spins,
    "EFPE_binary_black_hole": _optional_source_model(
        "pyEFPEHM", "EFPE_binary_black_hole"
    ),
}

PARAMETER_CONVERSIONS = {
    "lal_binary_black_hole": bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    "gwsignal_binary_black_hole": bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    "lal_eccentric_binary_black_hole_no_spins": bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    "lal_binary_neutron_star": bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    "EFPE_binary_black_hole": bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
}


def parse_collection(config, section, option, fallback=None):
    """Read a list or dict written with single quotes from an ini file.

    Parameters
    ----------
    config : configparser.ConfigParser
    section, option : str
    fallback : object or None
        Returned when the option is absent or empty.

    Returns
    -------
    object
    """
    if not config.has_option(section, option):
        return fallback
    raw = config.get(section, option).strip()
    if not raw:
        return fallback
    try:
        return json.loads(raw.replace("'", '"'))
    except json.JSONDecodeError as error:
        raise ValueError(
            "Could not parse [{}] {} = {!r}; expected a list or dict such as "
            "['CE20', 'CE40']".format(section, option, raw)
        ) from error


def read_injections(config, section="Injection_Parameters"):
    """Collect the sources to forecast.

    Accepts either parameters written inline in the configuration, or an
    ``injection-file`` pointing at a ``gwforge_population`` HDF5, optionally
    sliced with ``start-index`` and ``end-index``.

    Parameters
    ----------
    config : configparser.ConfigParser
    section : str

    Returns
    -------
    list of dict
        One dict of parameters per source.

    Raises
    ------
    ValueError
        If the section is missing or defines no source.
    """
    if not config.has_section(section):
        raise ValueError("The configuration has no [{}] section".format(section))

    if config.has_option(section, "injection-file"):
        path = config.get(section, "injection-file")
        start = config.getint(section, "start-index", fallback=0)
        end = config.getint(section, "end-index", fallback=None)
        return _read_injection_file(path, start, end)

    skip = {"injection-file", "start-index", "end-index"}
    parameters = {
        key: config.getfloat(section, key)
        for key in config.options(section)
        if key not in skip
    }
    if not parameters:
        raise ValueError(
            "[{}] defines neither an injection-file nor any inline "
            "parameters".format(section)
        )
    return [parameters]


def _read_injection_file(path, start=0, end=None):
    """Read a ``gwforge_population`` HDF5 into a list of parameter dicts."""
    logging.info("Reading injections from %s", path)
    with h5py.File(path, "r") as handle:
        # gwforge_population stores scalar-valued columns as length-1 datasets;
        # that is how the rest of the codebase detects an aligned-spin
        # population, and those columns are not per-source values.
        length = len(handle["mass_1"])
        columns = {
            key: handle[key][:] for key in handle.keys() if len(handle[key]) == length
        }
    if end is None:
        end = length
    return [
        {key: float(values[index]) for key, values in columns.items()}
        for index in range(start, min(end, length))
    ]


def waveform_settings(config, section="Waveform_Generator"):
    """Waveform arguments, source model and parameter conversion.

    Parameters
    ----------
    config : configparser.ConfigParser
    section : str

    Returns
    -------
    dict
        With keys ``waveform_arguments``, ``frequency_domain_source_model``,
        ``parameter_conversion`` and ``minimum_frequency``.

    Raises
    ------
    ValueError
        If the named source model is unknown or unavailable in this bilby.
    """
    minimum_frequency = config.getfloat(
        section, "waveform-minimum-frequency", fallback=20.0
    )
    reference_frequency = config.getfloat(
        section, "waveform-reference-frequency", fallback=minimum_frequency
    )
    approximant = config.get(section, "waveform-approximant", fallback="IMRPhenomXPHM")
    name = config.get(
        section, "frequency-domain-source-model", fallback="lal_binary_black_hole"
    )

    if name not in SOURCE_MODELS:
        raise ValueError(
            "Unknown frequency-domain-source-model {!r}; choose one of {}".format(
                name, ", ".join(sorted(SOURCE_MODELS))
            )
        )
    model = SOURCE_MODELS[name]
    if model is None:
        raise ValueError(
            "frequency-domain-source-model {!r} is not available in this "
            "installation of bilby ({})".format(name, bilby.__version__)
        )

    return {
        "waveform_arguments": {
            "waveform_approximant": approximant,
            "reference_frequency": reference_frequency,
            "minimum_frequency": minimum_frequency,
        },
        "frequency_domain_source_model": model,
        "parameter_conversion": PARAMETER_CONVERSIONS[name],
        "minimum_frequency": minimum_frequency,
    }


def fisher_settings(config, section="Fisher"):
    """Read the ``[Fisher]`` section.

    Parameters
    ----------
    config : configparser.ConfigParser
    section : str

    Returns
    -------
    dict
        With keys ``order``, ``earth_rotation``, ``finite_size``,
        ``parameters``, ``step_sizes``, ``calibrate`` and ``target_change``.
    """
    if not config.has_section(section):
        config.add_section(section)
    parameters = parse_collection(config, section, "parameters")
    return {
        "order": config.getint(section, "order", fallback=1),
        "earth_rotation": config.getboolean(section, "earth-rotation", fallback=True),
        "finite_size": config.getboolean(section, "finite-size", fallback=True),
        "parameters": list(parameters) if parameters else list(DEFAULT_PARAMETERS),
        "step_sizes": parse_collection(config, section, "step-sizes"),
        "calibrate": config.getboolean(section, "calibrate-steps", fallback=True),
        "estimate_noise": config.getboolean(section, "estimate-noise", fallback=True),
        "target_change": config.getfloat(
            section, "target-change", fallback=TARGET_FRACTIONAL_CHANGE
        ),
    }


def prior_settings(config, order, section="Priors"):
    """Read the ``[Priors]`` section and enforce the order-dependent rule.

    Order 1 is a closed-form covariance and needs no prior; order 2 samples a
    truncated Taylor expansion, which is only meaningful inside bounds, so a
    prior file is required there.

    Parameters
    ----------
    config : configparser.ConfigParser
    order : int
        Expansion order, from :func:`fisher_settings`.
    section : str

    Returns
    -------
    dict
        With keys ``prior_file`` and ``enforce_physicality``.

    Raises
    ------
    ValueError
        If ``order >= 2`` and no ``prior-file`` is given.
    """
    prior_file = None
    enforce = False
    if config.has_section(section):
        prior_file = config.get(section, "prior-file", fallback=None)
        enforce = config.getboolean(section, "enforce-physicality", fallback=False)

    if order >= 2 and not prior_file:
        raise ValueError(
            "order = {} samples the DALI posterior and needs bounds, so "
            "[Priors] prior-file is required. Point it at a bilby .prior file "
            "covering {} for every source.".format(order, "the Fisher parameters")
        )
    return {"prior_file": prior_file, "enforce_physicality": enforce}


def _aligned_spin(parameters, component):
    """Aligned spin component, however the parameters happen to spell it.

    ``convert_to_lal_binary_black_hole_parameters`` produces ``a_1`` and
    ``cos_tilt_1`` rather than ``spin_1z``, so reading ``spin_1z`` alone
    silently treats every binary as non-spinning and shortens the estimated
    duration.
    """
    for key in ("spin_{}z".format(component), "chi_{}".format(component)):
        value = parameters.get(key)
        if value is not None:
            return float(value)
    magnitude = parameters.get("a_{}".format(component))
    if magnitude is None:
        return 0.0
    cosine = parameters.get("cos_tilt_{}".format(component))
    if cosine is None:
        tilt = parameters.get("tilt_{}".format(component))
        cosine = numpy.cos(float(tilt)) if tilt is not None else 1.0
    return float(magnitude) * float(cosine)


def segment_duration(parameters, minimum_frequency, approximant, maximum=None):
    """Segment length long enough to hold the signal.

    Getting this wrong is not a subtle error: a segment shorter than the inspiral wraps it, and the
    resulting derivatives are meaningless.

    Parameters
    ----------
    parameters : dict
        Source parameters; needs component masses and aligned spins, which are
        derived from whatever mass parametrisation is present.
    minimum_frequency : float
        Waveform starting frequency in Hz.
    approximant : str
    maximum : float or None
        Cap on the returned duration, in seconds.

    Returns
    -------
    float
        Duration in seconds, rounded up to a power of two.
    """
    converted, _ = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(
        dict(parameters)
    )
    try:
        duration = get_safe_signal_durations(
            mass_1=numpy.atleast_1d(float(converted["mass_1"])),
            mass_2=numpy.atleast_1d(float(converted["mass_2"])),
            spin_1z=numpy.atleast_1d(_aligned_spin(converted, 1)),
            spin_2z=numpy.atleast_1d(_aligned_spin(converted, 2)),
            waveform_minimum_frequency=minimum_frequency,
            approximant=approximant,
        )
        duration = float(numpy.atleast_1d(duration)[0])
    except Exception as error:
        masses = antenna.chirp_mass_and_symmetric_mass_ratio(converted)
        if masses is None:
            raise ValueError(
                "Cannot size a segment for {}: {} rejected it and the masses "
                "could not be determined for the fallback estimate.".format(
                    approximant, get_safe_signal_durations.__name__
                )
            ) from error
        duration = float(
            numpy.atleast_1d(antenna.time_to_coalescence(minimum_frequency, *masses))[0]
        )
        logging.debug(
            "%s could not size %s (%s); using the 3.5PN estimate of %.1f s",
            get_safe_signal_durations.__name__,
            approximant,
            error,
            duration,
        )
    duration = duration + 2.0
    duration = 2.0 ** numpy.ceil(numpy.log2(max(duration, 4.0)))
    if maximum is not None:
        duration = min(duration, maximum)
    return float(duration)
