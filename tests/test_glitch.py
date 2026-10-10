"""Sine-Gaussian glitches: waveform, SNR scaling, reproducible draws, and the CLI on a file."""

import subprocess
import sys
from pathlib import Path

import h5py
import numpy
import pytest
from gwpy.timeseries import TimeSeries

from GWForge.glitch.sine_gaussian import EDGE, add_glitches, draw_glitches, envelope_width, optimal_snr, sine_gaussian
from GWForge.ifo.detectors import Network

RATE = 4096


@pytest.fixture(scope="module")
def ce40():
    return Network(ifos=["CE40"]).initialise_ifos()[0]


def test_sine_gaussian_peaks_at_time_and_dies_within_edge():
    t = numpy.arange(-2, 2, 1 / RATE)
    tau = envelope_width(50.0, 10.0)
    h = sine_gaussian(t, 0.0, 50.0, 10.0)
    assert h[numpy.argmin(numpy.abs(t))] == pytest.approx(1.0)
    assert numpy.abs(h[numpy.abs(t) > EDGE * tau]).max() < 1e-6
    # Morlet-Gabor: the other quadrature is just a phase shift
    assert sine_gaussian(t, 0.0, 50.0, 10.0, phase=-numpy.pi / 2) == pytest.approx(numpy.exp(-((t / tau) ** 2)) * numpy.sin(2 * numpy.pi * 50.0 * t))


def test_add_glitches_hits_requested_snr(ce40):
    samples = numpy.zeros(64 * RATE)
    glitches = {
        "time": numpy.array([20.0, 40.0]),
        "frequency": numpy.array([60.0, 300.0]),
        "q": numpy.array([12.0, 4.0]),
        "snr": numpy.array([50.0, 8.0]),
        "phase": numpy.array([0.0, 1.0]),
    }
    amplitudes = add_glitches(samples, 0.0, RATE, glitches, ce40)
    assert amplitudes.shape == (2,) and numpy.all(amplitudes > 0)
    for lo, hi, snr in ((0, 32 * RATE, 50.0), (32 * RATE, 64 * RATE, 8.0)):
        assert optimal_snr(samples[lo:hi], RATE, ce40) == pytest.approx(snr, rel=1e-2)
    assert samples[numpy.argmax(numpy.abs(samples[: 32 * RATE]))] == pytest.approx(amplitudes[0])


def test_draws_are_reproducible_and_clear_of_the_edges():
    span = dict(start=1e9, end=1e9 + 4096, rate=100, snr_range=[8, 100], frequency_range=[20, 1000], q_range=[2, 40])
    a = draw_glitches(numpy.random.default_rng([1, 0, 1_000_000_000]), **span)
    b = draw_glitches(numpy.random.default_rng([1, 0, 1_000_000_000]), **span)
    c = draw_glitches(numpy.random.default_rng([2, 0, 1_000_000_000]), **span)
    assert all(numpy.array_equal(a[key], b[key]) for key in a)
    assert len(a["time"]) > 50 and not numpy.array_equal(a["time"], c["time"])
    tau = envelope_width(a["frequency"], a["q"])
    assert numpy.all(a["time"] - EDGE * tau >= span["start"]) and numpy.all(a["time"] + EDGE * tau <= span["end"])
    assert 20 <= a["frequency"].min() and a["frequency"].max() <= 1000
    assert 8 <= a["snr"].min() and a["snr"].max() <= 100


def test_snr_power_law_and_discrete_phases():
    span = dict(start=0.0, end=3600.0, rate=20000, snr_range=[5, 200], frequency_range=[50, 500], q_range=[2, 6])
    g = draw_glitches(numpy.random.default_rng(0), snr_index=2.3, phases=[0, numpy.pi], **span)
    assert len(g["snr"]) > 10000 and 5 <= g["snr"].min() and g["snr"].max() <= 200
    assert 1 + len(g["snr"]) / numpy.log(g["snr"] / 5).sum() == pytest.approx(2.3, abs=0.2)
    assert set(numpy.unique(g["phase"])) == {0.0, numpy.pi}
    # index 1 is the log-uniform draw
    g = draw_glitches(numpy.random.default_rng(0), snr_index=1.0, **span)
    assert numpy.histogram(numpy.log(g["snr"]), bins=4)[0].std() / len(g["snr"]) < 0.02


def test_cli_glitches_a_zero_noise_file(tmp_path):
    start, duration = 1_000_000_000, 64
    (tmp_path / "CE40").mkdir()
    TimeSeries(numpy.zeros(duration * RATE), t0=start, dt=1 / RATE, name="CE40:INJ", channel="CE40:INJ").write(
        str(tmp_path / "CE40" / f"CE40-{start}-{duration}.h5"), format="hdf5"
    )
    ini = tmp_path / "glitches.ini"
    ini.write_text("[IFOS]\ndetectors = ['CE40']\n[Blip]\nrate = {'CE40': 1000}\nseed = 3\n[Tomte]\nrate = 1000\nfrequency-range = [30, 65]\nphases = [3.141592653589793]\n")
    script = Path(__file__).resolve().parent.parent / "bin" / "gwforge_glitch"
    result = subprocess.run(
        [sys.executable, str(script), "--config-file", str(ini), "--data-directory", str(tmp_path), "--gps-start-time", str(start), "--gps-end-time", str(start + duration)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    data = TimeSeries.read(str(tmp_path / "CE40" / f"CE40-{start}-{duration}.h5"), "CE40:INJ")
    with h5py.File(tmp_path / f"glitches-{start}.h5") as f:
        table = {key: f["CE40"][key][:] for key in f["CE40"]}
    assert set(table) == {"time", "frequency", "q", "snr", "phase", "amplitude", "class"}
    assert len(table["time"]) > 0 and all(len(v) == len(table["time"]) for v in table.values())
    tomte = table["class"] == b"Tomte"
    assert tomte.any() and (~tomte).any() and set(table["class"][~tomte]) == {b"Blip"}
    assert numpy.all((30 <= table["frequency"][tomte]) & (table["frequency"][tomte] <= 65)) and numpy.all(table["phase"][tomte] == numpy.pi)
    assert numpy.any(data.value != 0)
    # the loudest sample sits inside a recorded glitch
    loudest = data.times.value[numpy.argmax(numpy.abs(data.value))]
    assert numpy.min(numpy.abs(table["time"] - loudest)) < EDGE * envelope_width(table["frequency"], table["q"]).max()
