"""``gwforge_fisher`` configuration parsing, without evaluating a waveform.

These guard the pyEFPEHM plumbing: ``f22_start`` must reach the waveform only
when set (pyEFPEHM takes an explicit ``None`` over its own default) and must
reach the segment sizing whenever it is set.
"""

import configparser

import pytest

from GWForge.fisher.config import SOURCE_MODELS, segment_duration, waveform_settings

SOURCE = dict(mass_1=30.0, mass_2=25.0, chi_1=0.0, chi_2=0.0)


def _config(**options):
    config = configparser.ConfigParser()
    config["Waveform_Generator"] = {"waveform-approximant": "IMRPhenomXPHM", **options}
    return config


def test_f22_start_is_passed_only_when_set():
    assert "f22_start" not in waveform_settings(_config())["waveform_arguments"]
    arguments = waveform_settings(_config(**{"waveform-f22-start": "10"}))["waveform_arguments"]
    assert arguments["f22_start"] == 10.0


def test_extra_waveform_arguments_are_merged():
    arguments = waveform_settings(
        _config(**{"waveform-arguments-extra": "{'f22_ref': 10.0}"})
    )["waveform_arguments"]
    assert arguments["f22_ref"] == 10.0


def test_approximant_is_required():
    config = configparser.ConfigParser()
    config["Waveform_Generator"] = {}
    with pytest.raises(ValueError, match="waveform-approximant"):
        waveform_settings(config)


def test_efpehm_source_model_resolves():
    pyEFPEHM = pytest.importorskip("pyEFPEHM")
    assert SOURCE_MODELS["EFPEHM_binary_black_hole"] is pyEFPEHM.EFPE_binary_black_hole


def test_segment_duration_honours_f22_start_and_bounds():
    from_20 = segment_duration(SOURCE, 20.0, "IMRPhenomXPHM")
    from_10 = segment_duration(SOURCE, 20.0, "IMRPhenomXPHM", f22_start=10.0)
    assert from_10 > from_20
    assert segment_duration(SOURCE, 20.0, "IMRPhenomXPHM", minimum=64.0) == 64.0
    assert segment_duration(SOURCE, 20.0, "IMRPhenomXPHM", f22_start=10.0, maximum=4.0) == 4.0
