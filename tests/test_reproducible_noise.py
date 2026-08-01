"""GWForge's in-house noise generator vs PyCBC's, and its own invariants.

The headline question this file answers: **is GWForge's reimplementation
identical to PyCBC's, or merely similar?** It is identical to FFT round-off --
``max|dx| / rms`` lands around 1e-13, roughly nine orders of magnitude below the
statistical scatter you would see from a genuinely different realisation.

The PyCBC comparisons skip cleanly if PyCBC is not installed; the invariant tests
(reproducibility, chunk stitching, recovered PSD, zero noise) do not depend on it
and always run.
"""

import os

import numpy
import pytest

from GWForge.ifo import reproducible_noise

pycbc_noise = pytest.importorskip("pycbc.noise.reproduceable")
pycbc_psd = pytest.importorskip("pycbc.psd")
from pycbc.types import FrequencySeries  # noqa: E402

NOISE_CURVES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "GWForge",
    "ifo",
    "noise_curves",
)
CE40_ASD = os.path.join(NOISE_CURVES, "CE40-asd.txt")
ET_PSD = os.path.join(NOISE_CURVES, "ET-psd.txt")

GPS_START = 1893024018


def _load_psd(path, is_asd, delta_f, sample_rate, low_frequency_cutoff):
    """Tabulated curve -> uniform PSD grid, via GWForge's own interpolator."""
    table = numpy.loadtxt(path)
    frequencies, values = table[:, 0], table[:, -1]
    if is_asd:
        values = values**2
    length = int(sample_rate / 2 / delta_f) + 1
    return reproducible_noise.interpolate_psd(
        frequencies,
        values,
        length=length,
        delta_f=delta_f,
        low_frequency_cutoff=low_frequency_cutoff,
    )


def _relative_difference(reference, test):
    reference = numpy.asarray(reference)
    return numpy.abs(reference - numpy.asarray(test)).max() / reference.std()


# -- identity with PyCBC ------------------------------------------------------


@pytest.mark.parametrize(
    "path, is_asd, low_frequency_cutoff",
    [(CE40_ASD, True, 6.0), (ET_PSD, False, 6.0)],
)
@pytest.mark.parametrize(
    "sample_rate, filter_duration",
    [
        (1024, 8),  # delta_f (1/16) >= 1/(2*T): the interpolate-first branch
        (2048, 4),  # same branch, different resolution
        (1024, 32),  # delta_f (1/16) < 1/(2*T): the resize-to-segment branch
    ],
)
def test_colored_noise_matches_pycbc(
    path, is_asd, low_frequency_cutoff, sample_rate, filter_duration
):
    """Both ``filter_duration`` branches, for an ASD-file and a PSD-file curve."""
    delta_f = 1.0 / 16
    psd = _load_psd(path, is_asd, delta_f, sample_rate, low_frequency_cutoff)

    reference = pycbc_noise.colored_noise(
        FrequencySeries(psd.copy(), delta_f=delta_f),
        start_time=GPS_START,
        end_time=GPS_START + 8,
        seed=7,
        sample_rate=sample_rate,
        low_frequency_cutoff=low_frequency_cutoff,
        filter_duration=filter_duration,
    )
    samples, epoch = reproducible_noise.colored_noise(
        psd.copy(),
        delta_f=delta_f,
        start_time=GPS_START,
        end_time=GPS_START + 8,
        seed=7,
        sample_rate=sample_rate,
        low_frequency_cutoff=low_frequency_cutoff,
        filter_duration=filter_duration,
    )

    assert epoch == GPS_START
    assert len(samples) == len(reference)
    assert _relative_difference(reference.numpy(), samples) < 1e-10


def test_white_noise_matches_pycbc():
    """The white-noise draw itself, before any colouring."""
    reference = pycbc_noise.normal(GPS_START, GPS_START + 4, sample_rate=1024, seed=3)
    samples, epoch = reproducible_noise.normal(
        GPS_START, GPS_START + 4, sample_rate=1024, seed=3
    )
    assert epoch == GPS_START
    numpy.testing.assert_array_equal(samples, reference.numpy())


def test_interpolate_psd_matches_pycbc():
    reference = pycbc_psd.from_txt(
        filename=CE40_ASD,
        length=513,
        delta_f=1.0 / 16,
        low_freq_cutoff=6.0,
        is_asd_file=True,
    )
    table = numpy.loadtxt(CE40_ASD)
    mine = reproducible_noise.interpolate_psd(
        table[:, 0],
        table[:, -1] ** 2,
        length=513,
        delta_f=1.0 / 16,
        low_frequency_cutoff=6.0,
    )
    numpy.testing.assert_allclose(mine, reference.numpy(), rtol=1e-14, atol=0)


def test_resample_matches_pycbc():
    psd = _load_psd(CE40_ASD, True, 1.0 / 16, 1024, 6.0)
    reference = pycbc_psd.interpolate(
        FrequencySeries(psd.copy(), delta_f=1.0 / 16), 1.0 / 256
    )
    mine = reproducible_noise._resample(psd, 1.0 / 16, 1.0 / 256)
    numpy.testing.assert_allclose(mine, reference.numpy(), rtol=1e-14, atol=0)


def test_inverse_spectrum_truncation_matches_pycbc():
    """Tight in band; round-off-limited in the discarded roll-off.

    In band the typical bin agrees to an ulp. A handful below the cutoff --
    where ``1/|F|^2`` sits on a near-cancellation and numpy's and FFTW's
    summation orders diverge -- differ by up to ~1e-6 relative. Those bins are
    zeroed by :func:`colored_noise` before use, so the check is split: tight in
    the analysis band, loose below it.

    Nothing here asserts bit-identity. How many bins land on exactly the same
    double depends on the FFTW build and its SIMD codelets, so that count is a
    property of the machine rather than of the port.
    """
    delta_f = 1.0 / 256
    low_frequency_cutoff = 6.0
    psd = _load_psd(CE40_ASD, True, delta_f, 1024, low_frequency_cutoff)
    psd[psd == 0] = psd.max()
    reference = pycbc_psd.inverse_spectrum_truncation(
        FrequencySeries(psd.copy(), delta_f=delta_f),
        int(8 * 1024),
        low_frequency_cutoff=low_frequency_cutoff,
        trunc_method="hann",
    ).numpy()
    mine = reproducible_noise._inverse_spectrum_truncation(
        psd.copy(), delta_f, int(8 * 1024), low_frequency_cutoff
    )

    frequencies = numpy.arange(len(mine)) * delta_f
    band = frequencies >= low_frequency_cutoff
    numpy.testing.assert_allclose(mine[band], reference[band], rtol=1e-11, atol=0)

    nonzero = reference != 0
    relative = numpy.zeros_like(reference)
    relative[nonzero] = numpy.abs(mine[nonzero] - reference[nonzero]) / numpy.abs(
        reference[nonzero]
    )
    # The line above bounds every in-band bin; this says the typical one is far
    # tighter still, sitting at the one-ulp level rather than near that bound.
    assert numpy.median(relative[band]) < 1e-13
    assert relative[~band].max() < 1e-5


# -- reproducibility invariants ----------------------------------------------


def test_same_seed_and_gps_reproduce():
    kwargs = dict(sample_rate=1024, seed=11)
    first, _ = reproducible_noise.normal(GPS_START, GPS_START + 4, **kwargs)
    second, _ = reproducible_noise.normal(GPS_START, GPS_START + 4, **kwargs)
    numpy.testing.assert_array_equal(first, second)


def test_different_seed_gives_different_noise():
    first, _ = reproducible_noise.normal(
        GPS_START, GPS_START + 4, sample_rate=1024, seed=11
    )
    second, _ = reproducible_noise.normal(
        GPS_START, GPS_START + 4, sample_rate=1024, seed=12
    )
    assert not numpy.allclose(first, second)


def test_adjacent_chunks_stitch_seamlessly():
    """The property the Condor workflow depends on.

    Two chunks written by independent jobs must concatenate into exactly what a
    single long call would have produced -- otherwise every 4096 s boundary in a
    long observation carries a discontinuity.
    """
    long_span, _ = reproducible_noise.normal(
        GPS_START, GPS_START + 8, sample_rate=1024, seed=5
    )
    first, _ = reproducible_noise.normal(
        GPS_START, GPS_START + 4, sample_rate=1024, seed=5
    )
    second, _ = reproducible_noise.normal(
        GPS_START + 4, GPS_START + 8, sample_rate=1024, seed=5
    )
    numpy.testing.assert_array_equal(numpy.concatenate([first, second]), long_span)


def test_white_noise_stitches_across_a_block_boundary():
    """The white draw must stitch exactly even across a block boundary.

    ``BLOCK_SAMPLES / sample_rate`` = 1600 s at 1024 Hz, so a span straddling a
    multiple of 1600 s exercises the block-concatenation path rather than
    sitting inside a single block.
    """
    sample_rate = 1024
    block_duration = reproducible_noise.BLOCK_SAMPLES // sample_rate
    boundary = int(numpy.ceil(GPS_START / block_duration)) * block_duration

    whole, _ = reproducible_noise.normal(
        boundary - 4, boundary + 4, sample_rate=sample_rate, seed=5
    )
    first, _ = reproducible_noise.normal(
        boundary - 4, boundary, sample_rate=sample_rate, seed=5
    )
    second, _ = reproducible_noise.normal(
        boundary, boundary + 4, sample_rate=sample_rate, seed=5
    )
    numpy.testing.assert_array_equal(numpy.concatenate([first, second]), whole)


def test_colored_chunk_converges_to_the_long_span():
    """Coloured chunks approach a longer span as the filter lengthens.

    Unlike the white draw, coloured chunks are *not* slices of one another: the
    colouring is a circular convolution over the whole padded span, so a chunk of
    a different length is a different (statistically equivalent) realisation. The
    discrepancy is set by ``filter_duration`` and shrinks as it grows -- roughly
    2e-2 of the RMS at 8 s down to 1e-3 at 64 s. If it did *not* shrink, the
    difference would be a bug rather than a truncation effect.

    PyCBC behaves identically here, so this is not a GWForge regression; it is
    why ``GWForge.utils.split_duration`` pins a fixed 4096 s chunk length, which
    is what actually makes a generated dataset reproducible.
    """
    sample_rate = 1024
    psd = _load_psd(CE40_ASD, True, 1.0 / 16, sample_rate, 6.0)

    def discrepancy(filter_duration):
        kwargs = dict(
            delta_f=1.0 / 16,
            seed=5,
            sample_rate=sample_rate,
            low_frequency_cutoff=6.0,
            filter_duration=filter_duration,
        )
        whole, _ = reproducible_noise.colored_noise(
            psd.copy(), start_time=GPS_START, end_time=GPS_START + 16, **kwargs
        )
        first, _ = reproducible_noise.colored_noise(
            psd.copy(), start_time=GPS_START, end_time=GPS_START + 8, **kwargs
        )
        return _relative_difference(whole[: len(first)], first)

    short = discrepancy(8)
    long = discrepancy(64)
    assert short < 5e-2
    assert long < 5e-3
    assert long < short


# -- physical sanity, independent of PyCBC ------------------------------------


def test_recovered_psd_matches_input():
    """Welch-estimate the PSD back out: a units/normalisation check.

    Nothing here refers to PyCBC, so this catches a normalisation error that a
    faithful-but-wrong port would share with its reference.
    """
    sample_rate = 1024
    delta_f = 1.0 / 16
    low_frequency_cutoff = 20.0
    psd = _load_psd(CE40_ASD, True, delta_f, sample_rate, low_frequency_cutoff)

    samples, _ = reproducible_noise.colored_noise(
        psd.copy(),
        delta_f=delta_f,
        start_time=GPS_START,
        end_time=GPS_START + 512,
        seed=17,
        sample_rate=sample_rate,
        low_frequency_cutoff=low_frequency_cutoff,
        filter_duration=8,
    )

    # Welch with 8 s Hann segments, 50% overlap.
    segment = 8 * sample_rate
    window = numpy.hanning(segment)
    step = segment // 2
    periodograms = []
    for start in range(0, len(samples) - segment + 1, step):
        chunk = samples[start : start + segment] * window
        spectrum = numpy.abs(numpy.fft.rfft(chunk)) ** 2
        periodograms.append(2.0 * spectrum / (sample_rate * (window**2).sum()))
    estimated = numpy.median(periodograms, axis=0) / numpy.log(2)

    frequencies = numpy.arange(len(estimated)) * (sample_rate / segment)
    expected = numpy.interp(
        frequencies, numpy.arange(len(psd)) * delta_f, psd, left=0, right=0
    )
    band = (frequencies > 30) & (frequencies < 400)
    ratio = estimated[band] / expected[band]
    assert 0.9 < numpy.median(ratio) < 1.1


def test_zero_noise():
    samples, epoch = reproducible_noise.zero_noise(
        GPS_START, GPS_START + 4, sample_rate=1024
    )
    assert epoch == GPS_START
    assert len(samples) == 4 * 1024
    assert numpy.count_nonzero(samples) == 0
