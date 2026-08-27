"""The GWTC-5.0 medians are the defaults, and there is only one copy of them.

Every place that needs a fiducial BGP mass or Default BBH spin population used
to carry its own literal, and they had drifted: the Fisher model held a pre-O4b
*shape* on top of an O4b *support*, which is not a population any fit produced,
and the generation side had no BGP defaults at all. These tests pin all of it to
the model modules themselves, so the next drift fails here rather than in a
forecast six steps downstream.
"""

import configparser
import json
import os

import bilby
import numpy
import pytest

from GWForge.population import mass as mass_module
from GWForge.population import spin as spin_module
from GWForge.population.mass import BGP_PARAMETERS
from GWForge.population.redshift import GWTC5_LOCAL_MERGER_RATE_DENSITY
from GWForge.population.spin import DEFAULT_BBH_SPIN_PARAMETERS
from GWForge.population.mass import Mass
from GWForge.population.spin import Spin
from GWForge.population_fisher import BrokenPowerLawTwoPeakMass, DefaultSpin
from GWForge.population_fisher.config import build_model

PACKAGE_DIRECTORY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "GWForge"
)
CONFIGURATION_DIRECTORY = os.path.join(
    PACKAGE_DIRECTORY, "population_fisher", "configuration_files"
)
# The ``gwforge_population`` examples that ship with the package. These are
# what the docs point users at, so a broken one is a broken first experience.
POPULATION_CONFIGURATION_DIRECTORY = os.path.join(
    PACKAGE_DIRECTORY, "population", "population_configuration_files"
)
BGP_INI = os.path.join(POPULATION_CONFIGURATION_DIRECTORY, "bgp-gwtc5.ini")

# Every shipped population example, by filename.
POPULATION_CONFIGURATIONS = sorted(
    name
    for name in os.listdir(POPULATION_CONFIGURATION_DIRECTORY)
    if name.endswith(".ini")
)


@pytest.fixture(autouse=True)
def seeded():
    """The samplers draw from bilby priors as well as numpy's global RNG."""
    numpy.random.seed(20260824)
    bilby.core.utils.random.seed(20260824)


def parse(path):
    config = configparser.ConfigParser()
    config.read(path)
    return config


# ---------------------------------------------------------------------------
# The defaults are reachable at all
# ---------------------------------------------------------------------------


def test_bgp_samples_without_being_given_parameters():
    """``Mass("BGP", n)`` used to raise ``KeyError: alpha_1``.

    There was no BGP entry anywhere, so it fell through to the PowerLaw+Peak
    dictionary and died inside ``p_m1``.
    """
    samples = Mass("BGP", 2000).sample()
    assert len(samples["mass_1_source"]) == 2000
    assert samples["mass_1_source"].min() >= BGP_PARAMETERS["mmin"]
    secondaries = samples["mass_1_source"] * samples["mass_ratio"]
    assert secondaries.min() >= BGP_PARAMETERS["mmin_2"]


def test_default_spin_samples_without_being_given_parameters():
    samples = Spin("Default", 2000).sample()
    for index in (1, 2):
        magnitudes = samples["a_{}".format(index)]
        assert numpy.all((magnitudes >= 0.0) & (magnitudes <= 1.0))
        assert numpy.all(numpy.isfinite(samples["tilt_{}".format(index)]))


def test_default_spin_accepts_the_empty_dictionary_the_cli_passes():
    """``gwforge_population`` passes ``{}`` when an ini has no ``spin-parameters``.

    It logs "Assuming default values" as it does so, which was untrue until
    there were any.
    """
    assert Spin("Default", 5, {}).parameters["mu_t"] == DEFAULT_BBH_SPIN_PARAMETERS["mu_t"]


# ---------------------------------------------------------------------------
# One copy of the numbers
# ---------------------------------------------------------------------------


def test_generation_mass_defaults_are_the_medians():
    assert Mass("BGP", 1).parameters == BGP_PARAMETERS


def test_generation_spin_defaults_are_the_medians():
    parameters = Spin("Default", 1).parameters
    for name, value in DEFAULT_BBH_SPIN_PARAMETERS.items():
        assert parameters[name] == value


def test_fisher_mass_defaults_are_the_medians():
    model = BrokenPowerLawTwoPeakMass()
    assert model.fiducial == {
        name: BGP_PARAMETERS[name] for name in model.parameter_names
    }
    # Support parameters, which are fixed by construction rather than fitted.
    assert model.mmin == BGP_PARAMETERS["mmin"]
    assert model.mmin_2 == BGP_PARAMETERS["mmin_2"]
    assert model.delta_m_2 == BGP_PARAMETERS["delta_m_2"]
    assert model.m_high == BGP_PARAMETERS["m_high"]
    assert model.maximum_mass == BGP_PARAMETERS["maximum_mass"]


def test_fisher_spin_defaults_are_the_medians():
    model = DefaultSpin()
    assert model.fiducial == {
        name: DEFAULT_BBH_SPIN_PARAMETERS[name] for name in model.parameter_names
    }


def test_the_shipped_ini_matches_the_module():
    """The ini and the code are the two things that drifted apart last time."""
    config = parse(BGP_INI)
    mass = json.loads(config.get("Mass", "mass-parameters").replace("'", '"'))
    spin = json.loads(config.get("Spin", "spin-parameters").replace("'", '"'))
    assert mass == BGP_PARAMETERS
    assert spin == DEFAULT_BBH_SPIN_PARAMETERS
    assert config.getfloat("Redshift", "local-merger-rate-density") == pytest.approx(
        GWTC5_LOCAL_MERGER_RATE_DENSITY
    )


# Every shipped Fisher analysis config, by filename. Read from the directory
# rather than listed, so a config added later is covered without anyone
# remembering to add it here.
ANALYSIS_CONFIGURATIONS = sorted(
    name for name in os.listdir(CONFIGURATION_DIRECTORY) if name.endswith(".ini")
)


@pytest.mark.parametrize("name", ANALYSIS_CONFIGURATIONS)
def test_analysis_configs_use_the_current_spin_names(name):
    """``sigma_squared_chi`` is a *variance* and belongs to the Beta models.

    ``DefaultSpin`` takes ``sigma_chi``, a standard deviation, and a free
    ``mu_t``; ``_fiducial`` rejects anything else, so these configs used to
    fail to build.
    """
    path = os.path.join(CONFIGURATION_DIRECTORY, name)
    assert "sigma_squared_chi" not in open(path).read()
    config = parse(path)
    spin = json.loads(config.get("Model", "spin-parameters").replace("'", '"'))
    assert spin == DEFAULT_BBH_SPIN_PARAMETERS


# Spin models whose bilby prior is built by numerical integration --
# ``AlignedSpin`` with a non-uniform magnitude prior evaluates 10,000 scalar
# ``quad`` calls per component, about 90 s each, no matter how few samples are
# asked for. Their mass models are still exercised below.
SLOW_SPIN_MODELS = ("beta-aligned",)


@pytest.mark.parametrize("name", POPULATION_CONFIGURATIONS)
def test_the_shipped_population_examples_sample(name):
    """Every example in ``population_configuration_files`` still runs.

    ``precessing_bbh_powerlawpeak.ini`` shipped for a long time passing
    ``sigma_squared_chi`` to ``spin-model = Default``, which raises -- the
    example aborted at the spin step. Nothing covered this directory, so only a
    user running the file would ever have found out.

    ``Redshift`` is left out on purpose: it derives its own sample count from a
    full year of observing, which costs about a minute per configuration. The
    failure this guards against is a model/parameter-name mismatch, and that
    lives entirely in ``Mass`` and ``Spin``.
    """
    config = parse(os.path.join(POPULATION_CONFIGURATION_DIRECTORY, name))

    def parameters(section, option):
        """``None`` when the option is absent, the way the executable reads it."""
        if not config.has_option(section, option):
            return None
        return json.loads(config.get(section, option).replace("'", '"'))

    masses = Mass(
        mass_model=config.get("Mass", "mass-model"),
        number_of_samples=200,
        parameters=parameters("Mass", "mass-parameters"),
    ).sample()
    assert len(masses["mass_1_source"]) == 200

    spin_model = config.get("Spin", "spin-model")
    if spin_model.lower() in SLOW_SPIN_MODELS:
        pytest.skip(
            "{} builds its prior by numerical integration; see "
            "SLOW_SPIN_MODELS".format(spin_model)
        )
    spins = Spin(
        spin_model=spin_model,
        number_of_samples=200,
        parameters=parameters("Spin", "spin-parameters"),
    ).sample()
    assert spins
    assert all(len(value) == 200 for value in spins.values())


# ---------------------------------------------------------------------------
# Regressions in the Fisher config path
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ANALYSIS_CONFIGURATIONS)
def test_every_shipped_analysis_config_builds(name):
    """Each shipped example assembles a model whose free parameters exist.

    These are what the docs point users at, so a config naming a parameter the
    model does not have is a broken first experience rather than a subtle bug.
    """
    config = parse(os.path.join(CONFIGURATION_DIRECTORY, name))
    model, blocks, _, _ = build_model(config)
    assert blocks
    free = json.loads(
        config.get("Fisher", "free-parameters", fallback="[]").replace("'", '"')
    )
    unknown = [name for name in free if name not in model.parameter_names]
    assert not unknown, "{} frees {}, which the model does not have".format(
        name, unknown
    )


def test_a_config_with_a_spin_block_builds():
    """``build_model`` passed ``minimum_spin=``, which ``DefaultSpin`` never took.

    Every config naming the ``spin`` block raised ``TypeError`` at construction.
    """
    config = parse(os.path.join(CONFIGURATION_DIRECTORY, "mass_redshift_spin.ini"))
    _, blocks, sub_models, _ = build_model(config)
    assert "spin" in blocks
    assert sub_models["spin"].maximum_spin == 1.0
    assert sub_models["spin"].minimum_cosine_tilt == -1.0
    assert sub_models["spin"].fiducial == DEFAULT_BBH_SPIN_PARAMETERS


def test_the_built_mass_model_keeps_the_secondary_taper_independent():
    """The config used to override the class defaults and drop mmin_2/delta_m_2.

    That silently collapsed the secondary's taper onto the primary's, undoing
    the model fix on every Fisher run.
    """
    config = configparser.ConfigParser()
    config.add_section("Model")
    config.set("Model", "blocks", "['mass']")
    _, _, sub_models, _ = build_model(config)
    mass = sub_models["mass"]
    assert mass.mmin == BGP_PARAMETERS["mmin"]
    assert mass.mmin_2 == BGP_PARAMETERS["mmin_2"] != mass.mmin
    assert mass.delta_m_2 == BGP_PARAMETERS["delta_m_2"] != mass.fiducial["delta_m"]
    assert mass.m_high == BGP_PARAMETERS["m_high"]


# ---------------------------------------------------------------------------
# Nothing else moved
# ---------------------------------------------------------------------------


def test_power_law_peak_keeps_its_historical_defaults():
    assert Mass("PowerLaw+Peak", 1).parameters == {
        "alpha": 3.37,
        "beta": 0.76,
        "delta_m": 5.23,
        "mmin": 4.89,
        "mmax": 88.81,
        "lam": 0.04,
        "mpp": 33.60,
        "sigpp": 4.59,
    }


def test_models_without_an_entry_fall_back_as_before():
    """The PowerLaw+Peak set stays the catch-all, so no other model changed."""
    assert "multipeak" not in mass_module.DEFAULT_PARAMETERS
    assert (
        Mass("MultiPeak", 1).parameters == mass_module.POWER_LAW_PEAK_PARAMETERS
    )
    assert "nonspinning" not in spin_module.DEFAULT_PARAMETERS
    assert "mu_chi" not in Spin("nonspinning", 1).parameters


def test_the_defaults_cannot_be_mutated_through_an_instance():
    """``Mass`` assigned the shared default dictionary without copying it."""
    instance = Mass("BGP", 1)
    instance.parameters["alpha_1"] = -99.0
    assert BGP_PARAMETERS["alpha_1"] == 1.456442737
    assert Mass("BGP", 1).parameters["alpha_1"] == 1.456442737

    given = dict(DEFAULT_BBH_SPIN_PARAMETERS)
    Spin("Default", 1, given)
    assert given == DEFAULT_BBH_SPIN_PARAMETERS


def test_a_partial_dictionary_is_not_filled_in():
    """``p_m1`` takes ``**kwargs``, so a per-key fallback would hide a typo."""
    with pytest.raises((KeyError, TypeError)):
        Mass("BGP", 10, {"alpha_1": 2.0}).sample()
    with pytest.raises(ValueError, match="mu_t"):
        Spin("Default", 10, {"mu_chi": 0.1}).sample()
