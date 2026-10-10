"""Sine-Gaussian glitches, added incoherently: each detector gets its own.

A Morlet-Gabor wavelet is a Gaussian envelope times a complex exponential. Its
real part is the sine-Gaussian below with a phase offset, so one function covers
both choices of glitch model.
"""
import numpy

# Envelope widths kept clear of a file edge; exp(-EDGE**2) ~ 2e-16 is below double precision.
EDGE = 6


def envelope_width(frequency, q):
    """Gaussian envelope width ``tau = Q / (sqrt(2) pi f0)`` in seconds."""
    return q / (numpy.sqrt(2) * numpy.pi * frequency)


def sine_gaussian(times, time, frequency, q, amplitude=1.0, phase=0.0):
    """``h(t) = A exp(-((t - t0) / tau)^2) cos(2 pi f0 (t - t0) + phi)``.

    Parameters
    ----------
    times : numpy.ndarray
        Sample times in seconds.
    time, frequency, q : float
        Centre time (s), central frequency (Hz) and quality factor of the burst.
    amplitude, phase : float, optional
        Peak strain and phase at the centre; ``phase=-pi/2`` gives the sine quadrature.
    """
    tau = envelope_width(frequency, q)
    return amplitude * numpy.exp(-(((times - time) / tau) ** 2)) * numpy.cos(2 * numpy.pi * frequency * (times - time) + phase)


def optimal_snr(samples, sample_rate, interferometer):
    """Optimal SNR ``sqrt(4 df sum |h(f)|^2 / S(f))`` of ``samples`` against the detector's PSD.

    The transform is zero-padded to four times the input so a short burst is
    resolved in frequency; the band below the detector's minimum frequency is dropped.
    """
    n = 4 * len(samples)
    f = numpy.fft.rfftfreq(n, 1.0 / sample_rate)
    hf = numpy.fft.rfft(samples, n=n) / sample_rate
    psd = interferometer.power_spectral_density.power_spectral_density_interpolated(f)
    band = (f >= interferometer.minimum_frequency) & numpy.isfinite(psd)
    return numpy.sqrt(4 * (f[1] - f[0]) * numpy.sum(numpy.abs(hf[band]) ** 2 / psd[band]))


def draw_glitches(rng, start, end, rate, snr_range, frequency_range, q_range, snr_index=1.0, phases=None):
    """Draw a glitch population for one detector over ``[start, end]``.

    The number is Poisson with mean ``rate * (end - start) / 3600`` (``rate`` is per
    hour), the frequency is log-uniform in ``frequency_range``, the SNR follows the
    bounded power law ``p(snr) ~ snr**-snr_index`` on ``snr_range`` (``snr_index=1``
    is log-uniform), Q is uniform in ``q_range``, the phase is drawn from ``phases``
    when given (else uniform), and the time is uniform but kept ``EDGE`` envelope
    widths inside the span so no burst is cut by a file edge.

    The SNR is the optimal SNR against the detector's PSD. Omicron SNRs quoted in
    glitch catalogues are lower than this for the same burst.

    Returns
    -------
    dict of numpy.ndarray
        Keys ``time, frequency, q, snr, phase``, one entry per glitch.
    """
    n = rng.poisson(rate * (end - start) / 3600.0)
    frequency = numpy.exp(rng.uniform(*numpy.log(frequency_range), size=n))
    q = rng.uniform(*q_range, size=n)
    tau = envelope_width(frequency, q)
    low, high = snr_range
    u = rng.uniform(size=n)
    if numpy.isclose(snr_index, 1.0):
        snr = low * (high / low) ** u
    else:
        k = 1.0 - snr_index
        snr = (low**k + u * (high**k - low**k)) ** (1.0 / k)
    return {
        "time": rng.uniform(start + EDGE * tau, end - EDGE * tau),
        "frequency": frequency,
        "q": q,
        "snr": snr,
        "phase": rng.choice(numpy.asarray(phases, dtype=float), size=n) if phases is not None else rng.uniform(0, 2 * numpy.pi, size=n),
    }


def add_glitches(samples, epoch, sample_rate, glitches, interferometer):
    """Add every glitch in ``glitches`` to ``samples`` in place and return their peak strains.

    Each unit burst is built on ``time +- EDGE tau`` only and scaled so its
    optimal SNR against ``interferometer``'s PSD equals the requested ``snr``.
    """
    amplitudes = numpy.empty(len(glitches["time"]))
    for k, (time, frequency, q, snr, phase) in enumerate(zip(*(glitches[key] for key in ("time", "frequency", "q", "snr", "phase")))):
        half = int(numpy.ceil(EDGE * envelope_width(frequency, q) * sample_rate))
        centre = int(round((time - epoch) * sample_rate))
        lo, hi = max(centre - half, 0), min(centre + half + 1, len(samples))
        unit = sine_gaussian(epoch + numpy.arange(lo, hi) / sample_rate, time, frequency, q, phase=phase)
        amplitudes[k] = snr / optimal_snr(unit, sample_rate, interferometer)
        samples[lo:hi] += amplitudes[k] * unit
    return amplitudes
