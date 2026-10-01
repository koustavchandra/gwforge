"""Configuration handling for ``gwforge_population_fisher``.

Kept out of the executable so the parsing is importable and testable, and so the
CLI stays a thin driver -- the same split :mod:`GWForge.fisher.config` uses.

The five analyses differ only in ``[Model] blocks`` and ``[Fisher]
free-parameters``, so one config reader covers redshift-only through spectral
sirens without branching in the driver.
"""

import logging

from ..cosmology import COSMOLOGY_PARAMETERS, FlatwCDM, astropy_cosmology
from ..fisher.config import parse_collection
from ..population.mass import BGP_PARAMETERS
from .base import JointPopulationModel
from .model import (
    BrokenPowerLawTwoPeakMass,
    DefaultSpin,
    MadauDickinsonRedshift,
    SpectralSirenModel,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Population blocks a config may ask for. ``cosmology`` is not a block of its
# own -- it turns the model inside out, from a source-frame product into the
# detector-frame :class:`~GWForge.population_fisher.model.SpectralSirenModel`.
BLOCKS = ("mass", "redshift", "spin", "cosmology")

# Which model class owns which block, for the fit stage.
BLOCK_CLASSES = {
    "mass": BrokenPowerLawTwoPeakMass,
    "redshift": MadauDickinsonRedshift,
    "spin": DefaultSpin,
}


def catalogue_settings(config, section="Catalogue"):
    """Read the ``[Catalogue]`` section.

    Returns
    -------
    dict
        ``population_file``, ``snr_file``, ``detectors``, ``snr_thresholds``.
    """
    for option in ("population-file", "snr-file"):
        if not config.has_option(section, option):
            raise ValueError(
                "[{}] {} is required: the population Fisher needs both the "
                "injections and their SNRs.".format(section, option)
            )
    return {
        "population_file": config.get(section, "population-file"),
        "snr_file": config.get(section, "snr-file"),
        "detectors": parse_collection(config, section, "detectors", None),
        "snr_thresholds": [
            float(value)
            for value in parse_collection(
                config, section, "snr-thresholds", [10.0]
            )
        ],
    }


def cosmology_settings(config, section="Model"):
    """Build the fiducial :class:`GWForge.cosmology.FlatwCDM`.

    ``cosmology = <name>`` seeds it from an ``astropy`` realization, folding
    massive neutrinos into the matter density -- this background has only matter
    and dark energy, and at the redshifts a gravitational-wave catalogue reaches
    that is a far smaller error than pretending the neutrinos are not there.
    Explicit ``H0``/``Om0``/``w0`` override whatever the name supplied.

    Returns
    -------
    GWForge.cosmology.FlatwCDM
    """
    name = config.get(section, "cosmology", fallback="Planck18")
    try:
        cosmology = FlatwCDM.from_astropy(astropy_cosmology(name))
    except ValueError:
        cosmology = FlatwCDM()
    settings = cosmology.parameters()
    for parameter in COSMOLOGY_PARAMETERS:
        if config.has_option(section, parameter):
            settings[parameter] = config.getfloat(section, parameter)
    return FlatwCDM(**settings)


def build_model(config, section="Model"):
    """Assemble the population model the config describes.

    Parameters
    ----------
    config : configparser.ConfigParser
    section : str

    Returns
    -------
    tuple
        ``(model, blocks, sub_models, cosmology)``. ``sub_models`` maps block
        name to the source-frame model, which the fit stage needs even when the
        Fisher runs on the detector-frame wrapper.

    Raises
    ------
    ValueError
        For an unknown block, or ``cosmology`` without ``mass`` and ``redshift``.
    """
    blocks = parse_collection(config, section, "blocks", ["mass", "redshift"])
    unknown = [name for name in blocks if name not in BLOCKS]
    if unknown:
        raise ValueError(
            "Unknown [{}] block(s) {}; choose from {}.".format(
                section, unknown, list(BLOCKS)
            )
        )
    cosmology = cosmology_settings(config, section)
    maximum_redshift = config.getfloat(section, "maximum-redshift", fallback=10.0)

    sub_models = {}
    if "mass" in blocks:
        # Fall back to the model's own signature defaults -- the GWTC-5.0
        # medians -- rather than to a second set of literals here. mmin-2 and
        # delta-m-2 must be passed through or the secondary's independent taper
        # silently collapses onto the primary's.
        sub_models["mass"] = BrokenPowerLawTwoPeakMass(
            mmin=config.getfloat(section, "mmin", fallback=BGP_PARAMETERS["mmin"]),
            m_high=config.getfloat(
                section, "m-high", fallback=BGP_PARAMETERS["m_high"]
            ),
            maximum_mass=config.getfloat(
                section, "maximum-mass", fallback=BGP_PARAMETERS["maximum_mass"]
            ),
            mmin_2=config.getfloat(
                section, "mmin-2", fallback=BGP_PARAMETERS["mmin_2"]
            ),
            delta_m_2=config.getfloat(
                section, "delta-m-2", fallback=BGP_PARAMETERS["delta_m_2"]
            ),
            fiducial=_fiducial(config, section, "mass", BrokenPowerLawTwoPeakMass),
        )
    if "redshift" in blocks:
        sub_models["redshift"] = MadauDickinsonRedshift(
            maximum_redshift=maximum_redshift,
            cosmology=cosmology,
            fiducial=_fiducial(config, section, "redshift", MadauDickinsonRedshift),
        )
    if "spin" in blocks:
        # amax and t_min, both pinned in the GWTC-5.0 analysis. There is no
        # minimum-spin: the truncated Gaussian starts at 0 by construction.
        sub_models["spin"] = DefaultSpin(
            fiducial=_fiducial(config, section, "spin", DefaultSpin),
            maximum_spin=config.getfloat(section, "maximum-spin", fallback=1.0),
            minimum_cosine_tilt=config.getfloat(
                section, "minimum-cosine-tilt", fallback=-1.0
            ),
        )

    if "cosmology" in blocks:
        missing = [name for name in ("mass", "redshift") if name not in sub_models]
        if missing:
            raise ValueError(
                "A spectral-siren forecast needs the {} block(s) too: the "
                "cosmology is measured through source-frame mass features "
                "redshifting, so there has to be a mass model for them to be "
                "features of.".format(missing)
            )
        model = SpectralSirenModel(
            mass_model=sub_models["mass"],
            redshift_model=sub_models["redshift"],
            spin_model=sub_models.get("spin"),
            cosmology=cosmology,
            cosmology_parameters=parse_collection(
                config, section, "cosmology-parameters", list(COSMOLOGY_PARAMETERS)
            ),
        )
    else:
        model = JointPopulationModel(
            [sub_models[name] for name in blocks if name in sub_models]
        )
    return model, list(blocks), sub_models, cosmology


def _fiducial(config, section, block, model_class):
    """Fiducial hyper-parameters for one block, defaults filled in.

    A config may set only the values it cares about; anything omitted falls back
    to the class default, so a config need not restate all eleven BGP numbers to
    change one.
    """
    supplied = parse_collection(
        config, section, "{}-parameters".format(block), None
    )
    if supplied is None:
        return None
    defaults = dict(model_class().fiducial)
    unknown = [name for name in supplied if name not in defaults]
    if unknown:
        raise ValueError(
            "[{}] {}-parameters has unknown name(s) {}; {} takes {}.".format(
                section, block, unknown, model_class.__name__, sorted(defaults)
            )
        )
    defaults.update({name: float(value) for name, value in supplied.items()})
    return defaults


def fit_settings(config, section="Fit"):
    """Which blocks to refit to the injections, and which to only check.

    ``redshift`` defaults to ``True`` and should stay that way:
    :class:`GWForge.population.redshift.Redshift` convolves the star-formation
    rate with a time-delay distribution, so the generated redshifts are not
    distributed as the direct Madau-Dickinson form and the config's
    ``(gamma, kappa, z_peak)`` are not the fiducials of the model being
    forecast.

    ``mass`` and ``spin`` default to ``False``, because ``gwforge_population``
    draws those straight from the same densities this package evaluates. Turning
    them on is a consistency check, not a correction -- the fit must return the
    generation values.

    Returns
    -------
    dict
        Block name to bool.
    """
    return {
        "redshift": config.getboolean(section, "fit-redshift", fallback=True),
        "mass": config.getboolean(section, "fit-mass", fallback=False),
        "spin": config.getboolean(section, "fit-spin", fallback=False),
    }


def fisher_settings(config, section="Fisher"):
    """Read the ``[Fisher]`` section.

    ``configurations`` maps a name to *extra* parameters pinned on top of
    ``fixed-parameters``, so one config file can produce a degeneracy ladder --
    the joint fit, then LCDM with ``w0`` pinned, then ``H_0`` alone. That ladder
    is usually what separates two published forecasts, far more than the physics
    does, so it is worth reporting rather than choosing one rung silently.

    Returns
    -------
    dict
        ``free_parameters``, ``fixed_parameters``, ``configurations``,
        ``score_method``.
    """
    fixed = parse_collection(config, section, "fixed-parameters", {}) or {}
    configurations = parse_collection(config, section, "configurations", None)
    if configurations is None:
        configurations = {"": {}}
    return {
        "free_parameters": parse_collection(config, section, "free-parameters", None),
        "fixed_parameters": {name: float(value) for name, value in fixed.items()},
        "configurations": {
            name: {key: float(value) for key, value in (extra or {}).items()}
            for name, extra in configurations.items()
        },
        "score_method": config.get(section, "score-method", fallback="analytic"),
    }


def output_settings(config, section="Output", label="population_fisher"):
    """Read the ``[Output]`` section.

    Returns
    -------
    dict
        ``output_directory``, ``label``, ``plot_corner``.
    """
    return {
        "output_directory": config.get(
            section, "output-directory", fallback="population_fisher_output"
        ),
        "label": config.get(section, "label", fallback=label),
        "plot_corner": config.getboolean(section, "plot-corner", fallback=True),
        # A twenty-one-parameter corner is 441 panels and legible in none of
        # them; name the handful worth looking at instead.
        "plot_parameters": parse_collection(config, section, "plot-parameters", None),
    }
