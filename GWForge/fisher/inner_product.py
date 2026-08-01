
import numpy
from bilby.gw.utils import noise_weighted_inner_product


def inner_product(aa, bb, interferometer):
    """Noise-weighted inner product of two frequency-domain series.

    Wraps :func:`bilby.gw.utils.noise_weighted_inner_product` -- the same
    routine bilby uses for its own SNRs, so a Fisher SNR and an injection SNR
    are guaranteed to agree -- and takes the real part, which is what the
    Fisher matrix is built from.

    Only bins inside the interferometer's ``frequency_mask`` contribute. Outside
    it bilby leaves the PSD at an arbitrary sentinel value, so including those
    bins would silently corrupt the integral.

    Parameters
    ----------
    aa, bb : array_like
        Frequency-domain series on the interferometer's frequency grid. ``aa``
        is the one that gets conjugated. Either may be complex.
    interferometer : bilby.gw.detector.Interferometer
        Supplies the PSD, the frequency mask and the segment duration.

    Returns
    -------
    float
        The real inner product.
    """
    mask = interferometer.frequency_mask
    psd = interferometer.power_spectral_density_array[mask]
    return float(
        numpy.real(
            noise_weighted_inner_product(
                numpy.asarray(aa)[mask],
                numpy.asarray(bb)[mask],
                psd,
                interferometer.strain_data.duration,
            )
        )
    )


def optimal_snr_squared(signal, interferometer):
    """:math:`\\langle h \\mid h \\rangle` for a signal in one detector.

    Parameters
    ----------
    signal : array_like
        Frequency-domain detector response.
    interferometer : bilby.gw.detector.Interferometer
        Detector whose PSD is used.

    Returns
    -------
    float
    """
    return inner_product(signal, signal, interferometer)


def inner_product_matrix(left, right, interferometer):
    """Gram matrix of inner products between two stacks of series.

    Building the Fisher matrix means taking ``len(left) * len(right)`` inner
    products that all share the same PSD weighting, so the weighting is applied
    once and the rest is a single matrix product. For an 11-parameter Fisher
    over a 65537-bin grid this is roughly two orders of magnitude faster than
    looping over :func:`inner_product`.

    Parameters
    ----------
    left, right : array_like
        Arrays of shape ``(n_left, n_frequencies)`` and
        ``(n_right, n_frequencies)``. ``left`` is conjugated.
    interferometer : bilby.gw.detector.Interferometer
        Supplies the PSD, the frequency mask and the segment duration.

    Returns
    -------
    numpy.ndarray
        Real array of shape ``(n_left, n_right)``.
    """
    mask = interferometer.frequency_mask
    psd = interferometer.power_spectral_density_array[mask]
    duration = interferometer.strain_data.duration

    left = numpy.atleast_2d(numpy.asarray(left))[:, mask]
    right = numpy.atleast_2d(numpy.asarray(right))[:, mask]
    weighted = right / psd
    return (4.0 / duration) * numpy.real(numpy.conj(left) @ weighted.T)
