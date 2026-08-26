r"""Fisher matrix and doublet-DALI tensors for a detector network.

Expanding the Gaussian log-likelihood
:math:`\log L = -\tfrac{1}{2}\langle d - h(\theta)\mid d - h(\theta)\rangle`
about the fiducial parameters with :math:`d = h_0`, and writing
:math:`\Delta h = h_{,a}\Delta^a + \tfrac{1}{2}h_{,ab}\Delta^a\Delta^b + \ldots`,

.. math::

   \log L = -\tfrac{1}{2} F_{ab}\Delta^a\Delta^b
            - \tfrac{1}{2} G_{abc}\Delta^a\Delta^b\Delta^c
            - \tfrac{1}{8} H_{abcd}\Delta^a\Delta^b\Delta^c\Delta^d
            + \ldots

with

.. math::

   F_{ab} = \langle h_{,a}\mid h_{,b}\rangle, \quad
   G_{abc} = \langle h_{,ab}\mid h_{,c}\rangle, \quad
   H_{abcd} = \langle h_{,ab}\mid h_{,cd}\rangle .

The :math:`\tfrac{1}{2}` and :math:`\tfrac{1}{8}` come from the :math:`2\times
\tfrac{1}{2}` and :math:`(\tfrac{1}{2})^2` cross terms, reproducing eq. (15) of
Sellentin, Quartin & Amendola. Truncating at :math:`F` is the Fisher matrix of
Vallisneri; keeping all three is doublet-DALI.

The quartic term is :math:`-\tfrac{1}{8}\lVert h_{,ab}\Delta^a\Delta^b\rVert^2`,
which is negative-definite, so the doublet posterior is normalisable by
construction.

Third-order DALI ("triplet") is deliberately not implemented: there is no group
of terms beyond the three above that can be written with first and second
derivatives alone. See :func:`check_order`.
"""

import logging

import numpy

from .derivatives import WaveformDerivatives
from .parameters import DEFAULT_PARAMETERS, strip_shadowed_parameters
from .stepsize import TARGET_FRACTIONAL_CHANGE
from .inner_product import inner_product_matrix

#: Orders this module knows how to build.
SUPPORTED_ORDERS = (1, 2)

#: Above this the Fisher matrix is too ill-conditioned for the inverse to mean
#: anything in double precision, and the covariance is reported with a warning.
CONDITION_NUMBER_WARNING = 1e15


def check_order(order):
    """Validate the requested expansion order.

    Parameters
    ----------
    order : int
        ``1`` for the Fisher matrix, ``2`` for doublet-DALI.

    Returns
    -------
    int

    Raises
    ------
    NotImplementedError
        For ``order = 3``. Triplet-DALI needs genuine *third* derivatives of the
        waveform -- the terms :math:`\\langle h_{,a}\\mid h_{,bcd}\\rangle`,
        :math:`\\langle h_{,ab}\\mid h_{,cde}\\rangle` and
        :math:`\\langle h_{,abc}\\mid h_{,def}\\rangle` of Sellentin+ eq. (16).
        Grouping the expansion by highest derivative, the Fisher term plus
        :math:`G` and :math:`H` already exhausts everything expressible with
        first and second derivatives, so there is no consistent third order that
        avoids them.
    ValueError
        For any other unsupported value.
    """
    order = int(order)
    if order == 3:
        raise NotImplementedError(
            "order = 3 (triplet-DALI) requires third derivatives of the waveform "
            "-- the <h_,a|h_,bcd>, <h_,ab|h_,cde> and <h_,abc|h_,def> terms of "
            "Sellentin, Quartin & Amendola (arXiv:1401.6892) eq. (16). GWForge "
            "implements order 1 (Fisher) and order 2 (doublet-DALI), which "
            "together cover everything expressible with first and second "
            "derivatives."
        )
    if order not in SUPPORTED_ORDERS:
        raise ValueError(
            "order must be one of {}, got {}".format(SUPPORTED_ORDERS, order)
        )
    return order


class FisherMatrix:
    """Network Fisher matrix, and DALI tensors when ``order = 2``.

    Every tensor is a sum over detectors of noise-weighted inner products, so
    the network result is the sum of the per-detector results (Vallisneri
    eq. 13). ET contributes as three separate interferometers, which is how
    :class:`GWForge.ifo.detectors.Network` already presents it.

    Usage
    -----
    >>> fisher = FisherMatrix(ifos, waveform_generator, injection_parameters)
    >>> sigmas = numpy.sqrt(numpy.diag(fisher.covariance))

    Attributes
    ----------
    names : list of str
        Parameter names in matrix order.
    matrix : numpy.ndarray
        ``(n, n)`` network Fisher matrix.
    per_detector : dict
        Detector name to its own ``(n, n)`` Fisher matrix.
    doublet, quartic : numpy.ndarray or None
        :math:`G_{abc}` and :math:`H_{abcd}`; ``None`` when ``order < 2``.
    optimal_snrs : dict
        Detector name to optimal SNR, plus ``"network"``.
    """

    def __init__(
        self,
        interferometers,
        waveform_generator,
        parameters,
        fisher_parameters=None,
        earth_rotation=True,
        finite_size=True,
        step_sizes=None,
        order=1,
        calibrate=True,
        target_change=TARGET_FRACTIONAL_CHANGE,
        estimate_noise=True,
    ):
        """
        Parameters
        ----------
        interferometers : bilby.gw.detector.InterferometerList
            Network with strain data already set.
        waveform_generator : bilby.gw.WaveformGenerator
        parameters : dict
            Fiducial source parameters.
        fisher_parameters : sequence of str or None
            Defaults to :data:`GWForge.fisher.derivatives.DEFAULT_PARAMETERS`.
        earth_rotation, finite_size : bool
            Response options.
        step_sizes : dict or None
            Finite-difference step overrides, or starting guesses when
            ``calibrate`` is True.
        order : int
            ``1`` or ``2``; see :func:`check_order`.
        calibrate, target_change, estimate_noise :
            Automatic step selection and waveform-noise measurement; see
            :class:`GWForge.fisher.derivatives.WaveformDerivatives`.
        """
        self.order = check_order(order)
        # An injection read from a population file carries every mass and spin
        # parametrisation at once, and bilby's conversion prefers the component
        # masses -- so displacing chirp_mass would displace a key the waveform
        # ignores and hand back a zero derivative without raising.
        parameters, shadowed = strip_shadowed_parameters(
            parameters, fisher_parameters or DEFAULT_PARAMETERS
        )
        if shadowed:
            logging.info(
                "Dropped %s from the source parameters: they duplicate a Fisher "
                "parameter and the waveform would have used them instead.",
                ", ".join(shadowed),
            )
        self.parameters = dict(parameters)
        self.derivatives = WaveformDerivatives(
            interferometers,
            waveform_generator,
            parameters,
            fisher_parameters=fisher_parameters,
            earth_rotation=earth_rotation,
            finite_size=finite_size,
            step_sizes=step_sizes,
            order=self.order,
            calibrate=calibrate,
            target_change=target_change,
            estimate_noise=estimate_noise,
        )
        self.names = self.derivatives.names
        self._accumulate(interferometers)

    def _accumulate(self, interferometers):
        """Sum the per-detector tensors over the network."""
        count = len(self.names)
        self.matrix = numpy.zeros((count, count))
        self.per_detector = {}
        self.optimal_snrs = {}
        self.doublet = numpy.zeros((count,) * 3) if self.order >= 2 else None
        self.quartic = numpy.zeros((count,) * 4) if self.order >= 2 else None

        for interferometer in interferometers:
            response, first, second = self.derivatives.derivatives(interferometer)

            fisher = inner_product_matrix(first, first, interferometer)
            # Symmetrise: <h_,a|h_,b> is real and symmetric analytically, but
            # the two triangles are accumulated from different roundings.
            fisher = 0.5 * (fisher + fisher.T)
            self.per_detector[interferometer.name] = fisher
            self.matrix += fisher

            snr_squared = inner_product_matrix(
                response[numpy.newaxis, :], response[numpy.newaxis, :], interferometer
            )[0, 0]
            self.optimal_snrs[interferometer.name] = numpy.sqrt(max(snr_squared, 0.0))

            if self.order >= 2:
                flat = second.reshape(count * count, -1)
                self.doublet += inner_product_matrix(
                    flat, first, interferometer
                ).reshape(count, count, count)
                self.quartic += inner_product_matrix(
                    flat, flat, interferometer
                ).reshape(count, count, count, count)
            # Free the largest array before moving to the next detector.
            del first, second

        self.optimal_snrs["network"] = numpy.sqrt(
            sum(value**2 for value in self.optimal_snrs.values())
        )

    @property
    def covariance(self):
        """Inverse Fisher matrix, the linearised-signal covariance."""
        matrix, _ = covariance(self.matrix, names=self.names)
        return matrix

    @property
    def sigmas(self):
        """One-sigma marginalised uncertainties, as a name to value dict."""
        diagonal = numpy.diag(self.covariance)
        return dict(zip(self.names, numpy.sqrt(numpy.abs(diagonal))))

    def with_gaussian_prior(self, prior_sigmas):
        """Add Gaussian prior information to the Fisher matrix.

        A Gaussian prior of width :math:`\\sigma_a` contributes
        :math:`1/\\sigma_a^2` to the diagonal -- the only way a prior can enter a
        Fisher matrix, since a uniform prior carries no information in its
        interior.

        Parameters
        ----------
        prior_sigmas : dict
            Parameter name to prior standard deviation. Names not in
            :attr:`names` are ignored.

        Returns
        -------
        numpy.ndarray
            A new ``(n, n)`` matrix; :attr:`matrix` is left untouched.
        """
        combined = self.matrix.copy()
        for name, sigma in prior_sigmas.items():
            if name in self.names:
                combined[self.names.index(name), self.names.index(name)] += (
                    1.0 / float(sigma) ** 2
                )
        return combined


def covariance(fisher, names=None, condition_number_warning=CONDITION_NUMBER_WARNING):
    """Invert a Fisher matrix and report how much to trust the result.

    Fisher matrices for a gravitational-wave network are badly scaled -- a
    chirp-mass entry and a time entry differ by many orders of magnitude -- so
    the matrix is normalised by its diagonal before inversion and un-normalised
    afterwards. GWFast's ``CovMatr`` does the same thing for the same reason.

    Parameters
    ----------
    fisher : numpy.ndarray
        ``(n, n)`` symmetric positive-definite matrix.
    names : sequence of str or None
        Parameter names, used only to make the warnings readable.
    condition_number_warning : float
        Warn above this condition number.

    Returns
    -------
    tuple
        ``(covariance, diagnostics)``. ``diagnostics`` holds
        ``condition_number`` and ``inversion_error``, the latter being
        ``max|Sigma F - I|``.

    Raises
    ------
    numpy.linalg.LinAlgError
        If the matrix is singular even after normalisation.
    """
    fisher = numpy.asarray(fisher, dtype=float)
    diagonal = numpy.diag(fisher)
    if numpy.any(diagonal <= 0):
        bad = (
            [name for name, value in zip(names, diagonal) if value <= 0]
            if names is not None
            else numpy.flatnonzero(diagonal <= 0).tolist()
        )
        raise numpy.linalg.LinAlgError(
            "Fisher matrix has a non-positive diagonal entry for {}; the signal "
            "does not depend on those parameters at this point.".format(bad)
        )

    scale = 1.0 / numpy.sqrt(diagonal)
    normalised = fisher * numpy.outer(scale, scale)
    condition_number = numpy.linalg.cond(normalised)

    inverse = numpy.linalg.inv(normalised)
    result = inverse * numpy.outer(scale, scale)
    result = 0.5 * (result + result.T)

    inversion_error = numpy.max(
        numpy.abs(result @ fisher - numpy.eye(len(fisher)))
    )
    if condition_number > condition_number_warning:
        logging.warning(
            "Fisher matrix is ill-conditioned (condition number %.3e, inversion "
            "error %.3e). The covariance is unreliable; consider dropping a "
            "parameter or adding a prior.",
            condition_number,
            inversion_error,
        )
    return result, {
        "condition_number": float(condition_number),
        "inversion_error": float(inversion_error),
    }
