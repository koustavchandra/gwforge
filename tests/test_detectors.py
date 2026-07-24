import numpy
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
