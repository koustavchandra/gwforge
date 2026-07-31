"""
Reproducible coloured Gaussian noise and zero noise, in pure numpy.
Probably equivalent to ``pycbc.noise.reproduceable`` (agrees to FFT round-off), but without the PyCBC dependency.
The noise is reproducible in the sense that a given
``(seed, GPS time)`` always yields the same samples, so adjacent requests stitch together seamlessly.
This makes it suitable for generating long observations in chunks.
"""

import numpy
BLOCK_SAMPLES = 1638400
DEFAULT_FILTER_DURATION = 128


def block(seed, sample_rate):
    """Return one fixed-length block of white Gaussian samples.

    Parameters
    ----------
    seed : int
        Block seed. Reduced modulo ``2**32`` for ``RandomState``.
    sample_rate : float
        Sampling rate in Hz. Sets the variance to ``sample_rate / 2`` so that the
        one-sided PSD of the result is unity.

    Returns
    -------
    numpy.ndarray
        ``BLOCK_SAMPLES`` normally distributed samples.
    """
    generator = numpy.random.RandomState(seed % 2**32)
    variance = sample_rate / 2.0
    return generator.normal(size=BLOCK_SAMPLES, scale=variance**0.5)


def normal(start, end, sample_rate=16384, seed=0):
    """Return white Gaussian noise covering ``[start, end)`` in GPS seconds.

    The samples depend only on ``(seed, absolute GPS time)``: adjacent requests
    stitch together seamlessly, which is what lets the workflow generate a long
    observation as independent per-chunk Condor jobs.

    Returns
    -------
    tuple of (numpy.ndarray, float)
        The samples and the GPS epoch of the first one (which is ``start``).
    """
    block_duration = BLOCK_SAMPLES / sample_rate
    first = int(numpy.floor(start / block_duration))
    last = int(numpy.floor(end / block_duration))

    # The data evenly divides, so the last block would be superfluous.
    if end % block_duration == 0:
        last -= 1

    offset = numpy.random.RandomState(seed).randint(-(2**50), 2**50)
    samples = numpy.concatenate(
        [
            block(index + offset, sample_rate)
            for index in numpy.arange(first, last + 1, 1)
        ]
    )
    return time_slice(samples, first * block_duration, sample_rate, start, end)


def time_slice(samples, epoch, sample_rate, start, end):
    """Crop ``samples`` (starting at GPS ``epoch``) to ``[start, end)``.

    Returns
    -------
    tuple of (numpy.ndarray, float)
        The cropped samples and their GPS epoch, which is ``start``.
    """
    first = int(numpy.rint((start - epoch) * sample_rate))
    last = int(numpy.rint((end - epoch) * sample_rate))
    if first < 0 or last > len(samples):
        raise ValueError(
            "requested slice [{}, {}) is outside the generated span "
            "[{}, {})".format(start, end, epoch, epoch + len(samples) / sample_rate)
        )
    return samples[first:last], start


def zero_noise(start_time, end_time, sample_rate=16384):
    """Return an all-zero strain series covering ``[start_time, end_time)``.

    Same return contract as :func:`colored_noise`, so callers can switch between
    the two without special-casing anything.
    """
    length = int(numpy.rint((end_time - start_time) * sample_rate))
    return numpy.zeros(length), float(start_time)


def interpolate_psd(frequencies, psd, length, delta_f, low_frequency_cutoff):
    """Interpolate a tabulated PSD onto a uniform grid, log-log.

    Port of ``pycbc.psd.read.from_numpy_arrays``. Interpolation is linear in
    ``(log f, log S)``, which is the right space for noise curves spanning many
    decades. Bins below ``low_frequency_cutoff`` are set to zero.

    Parameters
    ----------
    frequencies, psd : array_like
        Tabulated one-sided PSD. ``frequencies`` must be ascending and positive.
    length : int
        Number of samples in the output (DC to Nyquist inclusive).
    delta_f : float
        Frequency resolution of the output in Hz.
    low_frequency_cutoff : float
        Bins below this are zeroed.

    Returns
    -------
    numpy.ndarray
        The PSD on ``numpy.arange(length) * delta_f``.
    """
    frequencies = numpy.asarray(frequencies, dtype=float)
    psd = numpy.asarray(psd, dtype=float)

    if frequencies[0] > low_frequency_cutoff:
        raise ValueError(
            "lowest frequency in the input PSD ({} Hz) is above the requested "
            "low-frequency cutoff ({} Hz)".format(frequencies[0], low_frequency_cutoff)
        )

    kmin = int(low_frequency_cutoff / delta_f)
    flow = kmin * delta_f

    # Start one tabulated point below the cutoff so the first interpolated bin is
    # bracketed rather than extrapolated.
    if frequencies[0] == low_frequency_cutoff:
        start = 0
    else:
        start = max(0, numpy.searchsorted(frequencies, flow) - 1)
    if frequencies[start + 1] == low_frequency_cutoff:
        start += 1

    frequencies = frequencies[start:]
    psd = psd[start:]

    if (length - 1) * delta_f > frequencies[-1]:
        length = int(frequencies[-1] / delta_f + 1)

    interpolated = numpy.zeros(length, dtype=numpy.float64)
    # numpy.interp clamps outside the input range, matching PyCBC's
    # scipy.interpolate.interp1d(fill_value=(slog[0], slog[-1])).
    interpolated[kmin:] = numpy.exp(
        numpy.interp(
            numpy.log(numpy.arange(kmin, length) * delta_f),
            numpy.log(frequencies),
            numpy.log(psd),
        )
    )
    return interpolated


def _resample(series, delta_f, new_delta_f, length=None):
    """Linearly resample a uniformly sampled series to a new ``delta_f``.

    Port of ``pycbc.psd.interpolate``.
    """
    if length is None:
        length = int(numpy.rint((len(series) - 1) * delta_f / new_delta_f + 1))
    return numpy.interp(
        numpy.arange(length) * new_delta_f,
        numpy.arange(len(series)) * delta_f,
        series,
    )


def _inverse_spectrum_truncation(psd, delta_f, max_filter_length, low_frequency_cutoff):
    """Coarse-grain a PSD by truncating its inverse ASD's impulse response.

    Port of ``pycbc.psd.estimate.inverse_spectrum_truncation`` restricted to the
    two options GWForge uses: ``which_spectrum='invasd'`` and
    ``trunc_method='hann'``. See arXiv:gr-qc/0509116 for the method.

    The goal is to remove the sharp spectral features that
    would otherwise give the colouring filter an impulse response longer than the
    data segment.
    """
    n_time = (len(psd) - 1) * 2

    inverse = numpy.zeros(len(psd), dtype=complex)
    kmin = 1
    if low_frequency_cutoff:
        kmin = int(low_frequency_cutoff / delta_f)
    inverse[kmin : n_time // 2] = 1.0 / psd[kmin : n_time // 2]
    inverse[: n_time // 2] = inverse[: n_time // 2] ** 0.5

    impulse = numpy.fft.irfft(inverse, n=n_time) * n_time * delta_f

    start = max_filter_length // 2
    end = n_time - max_filter_length // 2
    if end < start:
        raise ValueError(
            "max_filter_length ({}) exceeds the PSD's time-domain length "
            "({})".format(max_filter_length, n_time)
        )

    window = numpy.hanning(max_filter_length)
    impulse[0:start] *= window[-start:]
    impulse[end:n_time] *= window[0 : max_filter_length // 2]
    if start < end:
        impulse[start:end] = 0

    delta_t = 1.0 / (n_time * delta_f)
    truncated = numpy.fft.rfft(impulse) * delta_t
    return 1.0 / numpy.abs(truncated * numpy.conj(truncated))


def colored_noise(
    psd,
    delta_f,
    start_time,
    end_time,
    seed=0,
    sample_rate=16384,
    low_frequency_cutoff=1.0,
    filter_duration=DEFAULT_FILTER_DURATION,
    scale=1.0,
):
    """Return Gaussian noise coloured by ``psd``, covering ``[start, end)``.

    Numerically equivalent to ``pycbc.noise.reproduceable.colored_noise`` (agrees
    to FFT round-off), including the seed/GPS reproducibility contract described
    in the module docstring.

    Parameters
    ----------
    psd : array_like
        One-sided PSD sampled uniformly from DC at spacing ``delta_f``.
    delta_f : float
        Frequency resolution of ``psd`` in Hz.
    start_time, end_time : float
        GPS bounds of the requested data. The returned series covers
        ``[start_time, end_time)``.
    seed : int
        Noise seed. A given ``(seed, GPS time)`` always yields the same samples.
    sample_rate : float
        Output sampling rate in Hz. Must be held constant to keep continuity
        between disjoint spans.
    low_frequency_cutoff : float
        The PSD is zeroed below this, so the output carries no power there.
    filter_duration : float
        Length in seconds of the colouring filter.
    scale : float
        Extra multiplicative factor on the amplitude spectrum.

    Returns
    -------
    tuple of (numpy.ndarray, float)
        The coloured samples and their GPS epoch (which is ``start_time``).
    """
    psd = numpy.array(psd, dtype=float, copy=True)
    target_length = int(sample_rate / delta_f) // 2 + 1
    old_length = len(psd)
    if target_length > old_length:
        psd = numpy.concatenate([psd, numpy.zeros(target_length - old_length)])
    else:
        psd = psd[:target_length]
    maximum = psd.max()
    psd[old_length - 1 :] = psd[old_length - 2]
    psd[psd == 0] = maximum

    filter_length = int(filter_duration * sample_rate)
    white_duration = int(end_time - start_time) + 2 * filter_duration

    if delta_f >= 1.0 / (2.0 * filter_duration):
        filter_delta_f = 1.0 / (2.0 * filter_duration)
        psd = _resample(psd, delta_f, filter_delta_f)
        delta_f = filter_delta_f
        psd = 1.0 / _inverse_spectrum_truncation(
            1.0 / psd, delta_f, filter_length, low_frequency_cutoff
        )
        n_time = (len(psd) - 1) * 2
        impulse = numpy.fft.irfft(psd.astype(complex), n=n_time) * n_time * delta_f
        impulse = numpy.roll(impulse, filter_length)
        padded_length = int(white_duration * sample_rate)
        if padded_length > len(impulse):
            impulse = numpy.concatenate(
                [impulse, numpy.zeros(padded_length - len(impulse))]
            )
        else:
            impulse = impulse[:padded_length]
        impulse = numpy.roll(impulse, -filter_length)
        delta_t = 1.0 / sample_rate
        psd = numpy.fft.rfft(impulse) * delta_t
        delta_f = 1.0 / (len(impulse) * delta_t)
    else:
        psd = _resample(psd, delta_f, 1.0 / white_duration)
        delta_f = 1.0 / white_duration
        psd = 1.0 / _inverse_spectrum_truncation(
            1.0 / psd, delta_f, filter_length, low_frequency_cutoff
        )

    kmin = int(low_frequency_cutoff / delta_f)
    psd[:kmin] = 0
    amplitude_spectrum = numpy.abs(psd) ** 0.5

    white, white_epoch = normal(
        start_time - filter_duration,
        end_time + filter_duration,
        sample_rate=sample_rate,
        seed=seed,
    )
    delta_t = 1.0 / sample_rate
    colored = numpy.fft.rfft(white) * delta_t
    colored = colored * amplitude_spectrum * scale
    colored = numpy.fft.irfft(colored, n=len(white)) / delta_t
    return time_slice(colored, white_epoch, sample_rate, start_time, end_time)


def noise_from_interferometer(
    interferometer,
    start_time,
    end_time,
    seed=0,
    sample_rate=16384,
    noise_type="gaussian",
    low_frequency_cutoff=None,
    filter_duration=DEFAULT_FILTER_DURATION,
    psd_delta_f=1.0 / 16,
):
    """Generate noise for a bilby ``Interferometer``.

    The PSD is taken from ``interferometer.power_spectral_density`` as *arrays*
    (``frequency_array`` / ``psd_array``) rather than by re-reading the file. Those
    arrays are populated no matter how the PSD was supplied -- PSD file, ASD file,
    or raw arrays -- so this works for ET (which ships a PSD file) and for every
    user-supplied ``psd-dict`` entry, both of which the old file-sniffing path
    could not handle.

    Parameters
    ----------
    interferometer : bilby.gw.detector.Interferometer
        Detector whose PSD and ``minimum_frequency`` are used.
    start_time, end_time : float
        GPS bounds of the requested data.
    seed : int
        Noise seed; ignored for zero noise.
    sample_rate : float
        Output sampling rate in Hz.
    noise_type : str
        ``'gaussian'`` or ``'zero'``.
    low_frequency_cutoff : float or None
        Defaults to the interferometer's ``minimum_frequency``.
    filter_duration : float
        Length in seconds of the colouring filter.
    psd_delta_f : float
        Resolution the tabulated PSD is interpolated onto before colouring.

    Returns
    -------
    tuple of (numpy.ndarray, float)
        The strain samples and their GPS epoch.
    """
    noise_type = str(noise_type).lower()
    if "zero" in noise_type:
        return zero_noise(start_time, end_time, sample_rate=sample_rate)
    if "gaussian" not in noise_type:
        raise ValueError(
            "unsupported noise type {!r}: expected 'gaussian' or 'zero'".format(
                noise_type
            )
        )

    if low_frequency_cutoff is None:
        low_frequency_cutoff = interferometer.minimum_frequency

    power_spectral_density = interferometer.power_spectral_density
    length = int(sample_rate / 2 / psd_delta_f) + 1
    psd = interpolate_psd(
        power_spectral_density.frequency_array,
        power_spectral_density.psd_array,
        length=length,
        delta_f=psd_delta_f,
        low_frequency_cutoff=low_frequency_cutoff,
    )
    return colored_noise(
        psd,
        delta_f=psd_delta_f,
        start_time=start_time,
        end_time=end_time,
        seed=seed,
        sample_rate=sample_rate,
        low_frequency_cutoff=low_frequency_cutoff,
        filter_duration=filter_duration,
    )
