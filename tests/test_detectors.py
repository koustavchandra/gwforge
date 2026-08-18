import numpy
import pytest

import bilby
from GWForge.ifo.detectors import IFO, Network


def test_initialise_ifo_returns_interferometer():
    ifo = IFO(name="CE40").initialise_ifo()
    assert isinstance(ifo, bilby.gw.detector.Interferometer)


def test_initialise_ifo_seed_is_reproducible():
    # With the seed honoured, subsequent draws from the global RNG must be
    # identical across two seeded initialisations. (The old code called
    # numpy.random.seed(None), ignoring the argument and destroying reproducibility.)
    IFO(name="CE40").initialise_ifo(seed=42)
    first = numpy.random.rand(5)
    IFO(name="CE40").initialise_ifo(seed=42)
    second = numpy.random.rand(5)
    assert numpy.allclose(first, second)


def test_network_returns_interferometer_list():
    network = Network(ifos=["CE40", "CE20"]).initialise_ifos()
    assert isinstance(network, bilby.gw.detector.InterferometerList)
    assert len(network) == 2


@pytest.mark.parametrize("name", ["H1", "L1", "V1"])
def test_current_generation_detectors_fall_back_to_bilby(name):
    """Non-XG detectors must come from bilby, with a usable PSD.

    These have no ``.ifo`` file and no bundled noise curve, so they take the
    ``FileNotFoundError`` fallback. That fallback was unreachable for a while:
    the XG branch loads its noise curve *before* the ``.ifo`` file (because
    ``load_xg_interferometer`` needs a valid PSD up front), and bilby's PSD
    loader raises ``ValueError`` rather than ``FileNotFoundError`` for a missing
    file -- so every one of these was reported as "not implemented", breaking
    HLV for all four CLIs.
    """
    ifo = IFO(name=name).initialise_ifo()
    assert isinstance(ifo, bilby.gw.detector.Interferometer)
    assert ifo.name == name
    assert ifo.power_spectral_density.power_spectral_density_interpolated(100.0) > 0


def test_hlv_network():
    network = Network(ifos=["H1", "L1", "V1"]).initialise_ifos()
    assert [ifo.name for ifo in network] == ["H1", "L1", "V1"]


def test_asharp_psd_is_attached():
    """The A# branch lives inside the fallback, so it broke along with it."""
    ifo = IFO(name="H1", asharp=True).initialise_ifo()
    assert ifo.power_spectral_density.asd_file.endswith("Asharp-asd.txt")
    plain = IFO(name="H1").initialise_ifo()
    assert ifo.power_spectral_density.power_spectral_density_interpolated(
        100.0
    ) != plain.power_spectral_density.power_spectral_density_interpolated(100.0)


def test_triangular_detector_expands_to_three_arms():
    """ET is one .ifo file but three interferometers, each carrying the PSD."""
    network = Network(ifos=["ET"]).initialise_ifos()
    assert [ifo.name for ifo in network] == ["ET1", "ET2", "ET3"]
    for ifo in network:
        assert ifo.power_spectral_density.power_spectral_density_interpolated(100.0) > 0


def test_unknown_detector_is_reported_clearly():
    with pytest.raises(ValueError, match="not implemented"):
        IFO(name="NOT_A_DETECTOR").initialise_ifo()
