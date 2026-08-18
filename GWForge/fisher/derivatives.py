r"""First and second derivatives of the detector response.

Most of the work is deciding which derivatives have to be numerical. A bilby
source model returns :math:`\tilde h_+` and :math:`\tilde h_\times` from the
intrinsic parameters plus ``theta_jn``, ``phase`` and ``luminosity_distance``,
and has **no** dependence on ``ra``, ``dec``, ``psi`` or ``geocent_time`` -- so
those four, plus ``luminosity_distance``, are exact. See
:mod:`GWForge.fisher.antenna_derivatives` for the geometry and
:mod:`GWForge.fisher.stepsize` for how the numerical steps are chosen.

Everything else is differenced with fourth-order five-point stencils, which give
the first and second derivative from the same five evaluations. Note ``phase``
is *not* analytic for a higher-mode waveform: each :math:`(\ell, m)` carries
:math:`e^{-i m \phi_c}` and the mode decomposition is not exposed.

Mixed second derivatives between a numerical and an analytic parameter are taken
by differencing the *exact analytic* derivative along the numerical stencil.
That matters when ``earth_rotation`` is on, because the masses then reach the
response through :math:`\tau(f)` as well as through the waveform.
"""

import logging

import numpy

from ..ifo import antenna
from .antenna_derivatives import AntennaDerivatives
from .inner_product import inner_product
from .parameters import (
    ANALYTIC_PARAMETERS,
    DEFAULT_PARAMETERS,
    bounded_step,
    step_size,
)
from .stepsize import (
    NOISE_MARGIN,
    TARGET_FRACTIONAL_CHANGE,
    _NOISE_WARNING_ERROR,
    calibrate_step,
    estimate_noise,
    sample_around,
)

def _masked_frequencies(interferometer):
    """In-band frequencies and the mask that selects them."""
    mask = interferometer.frequency_mask
    return interferometer.frequency_array[mask], mask


def _times_to_coalescence(parameters, frequencies, earth_rotation, interferometer):
    """:math:`\\tau(f)` for the rotating response, or None if not needed."""
    if not earth_rotation:
        return None
    masses = antenna.chirp_mass_and_symmetric_mass_ratio(parameters)
    if masses is None:
        logging.warning(
            "%s: cannot determine the masses from the source parameters, so the "
            "time-frequency relation is unavailable; falling back to a static "
            "antenna pattern for the derivatives.",
            interferometer.name,
        )
        return None
    return antenna.time_to_coalescence(frequencies, *masses)


def analytic_derivatives(
    interferometer,
    polarizations,
    parameters,
    earth_rotation=True,
    finite_size=True,
    order=1,
):
    r"""Exact response and derivatives for the five closed-form parameters.

    Writing :math:`h = P\,S` with :math:`P = F_+\tilde h_+ + F_\times\tilde h_\times`
    and :math:`S = e^{-2\pi i f\tau}`,

    .. math::

       \partial_a h = \big(\partial_a P - 2\pi i f\,\partial_a\tau\,P\big) S,

    and the second derivative follows from one more product rule. ``ra``,
    ``dec``, ``psi`` and ``geocent_time`` enter through
    :class:`AntennaDerivatives`; ``luminosity_distance`` is handled separately
    because :math:`h \propto 1/d_L` makes every one of its derivatives a rescaling
    of a derivative we already have.

    Parameters
    ----------
    interferometer : bilby.gw.detector.Interferometer
    polarizations : dict
        ``{'plus': ..., 'cross': ...}`` on the interferometer's frequency grid.
    parameters : dict
        Source parameters.
    earth_rotation, finite_size : bool
        Response options; must match those used elsewhere.
    order : int
        ``2`` also returns second derivatives.

    Returns
    -------
    tuple
        ``(response, first, second)``. ``first`` maps a parameter name to a
        full-length complex array; ``second`` maps a sorted name pair to one.
        Both are zero outside the interferometer's frequency mask.
    """
    frequencies, mask = _masked_frequencies(interferometer)
    times_to_coalescence = _times_to_coalescence(
        parameters, frequencies, earth_rotation, interferometer
    )
    rotating = earth_rotation and times_to_coalescence is not None

    antenna_derivatives = AntennaDerivatives(
        interferometer,
        ra=parameters["ra"],
        dec=parameters["dec"],
        geocent_time=parameters["geocent_time"],
        psi=parameters["psi"],
        frequencies=frequencies,
        start_time=interferometer.strain_data.start_time,
        times_to_coalescence=times_to_coalescence,
        earth_rotation=rotating,
        finite_size=finite_size,
    )

    plus, cross = polarizations["plus"][mask], polarizations["cross"][mask]
    phase_factor = -2j * numpy.pi * frequencies

    def projected(key):
        return antenna_derivatives.plus[key] * plus + antenna_derivatives.cross[key] * cross

    shift = numpy.exp(phase_factor * antenna_derivatives.delay["value"])
    band_response = projected("value") * shift

    def embed(values):
        full = numpy.zeros(len(mask), dtype=complex)
        full[mask] = values
        return full

    geometric = ("ra", "dec", "psi", "geocent_time")
    band_first = {
        name: (
            projected(name)
            + phase_factor * antenna_derivatives.delay[name] * projected("value")
        )
        * shift
        for name in geometric
    }

    response = embed(band_response)
    first = {name: embed(values) for name, values in band_first.items()}
    distance = parameters["luminosity_distance"]
    first["luminosity_distance"] = -response / distance

    second = {}
    if order >= 2:
        for i, left in enumerate(geometric):
            for right in geometric[i:]:
                pair = tuple(sorted((left, right)))
                delay_left = antenna_derivatives.delay[left]
                delay_right = antenna_derivatives.delay[right]
                band = (
                    projected(pair)
                    + phase_factor * delay_left * projected(right)
                    + phase_factor * delay_right * projected(left)
                    + phase_factor * antenna_derivatives.delay[pair] * projected("value")
                    + phase_factor**2 * delay_left * delay_right * projected("value")
                ) * shift
                second[pair] = embed(band)
        # h = g(other parameters) / d_L, so every distance second derivative is
        # a rescaling of something already computed.
        for name in geometric:
            second[tuple(sorted((name, "luminosity_distance")))] = -first[name] / distance
        second[("luminosity_distance", "luminosity_distance")] = (
            2.0 * response / distance**2
        )

    return response, first, second


#: Offsets and weights of the fourth-order five-point first-derivative stencil,
#: in units of the step size. The centre point has zero weight here but is
#: reused by the second-derivative stencil below.
_FIRST_STENCIL = ((-2, 1.0), (-1, -8.0), (1, 8.0), (2, -1.0))
_FIRST_DENOMINATOR = 12.0

#: Fourth-order five-point second-derivative stencil, sharing the same points.
_SECOND_STENCIL = ((-2, -1.0), (-1, 16.0), (0, -30.0), (1, 16.0), (2, -1.0))
_SECOND_DENOMINATOR = 12.0


class WaveformDerivatives:
    r"""First and second derivatives of the detector response for one source.

    Combines the exact closed-form derivatives of :func:`analytic_derivatives`
    with fourth-order finite differences for the parameters the waveform
    actually depends on.

    The waveform is evaluated once per stencil point and cached, because the
    polarizations do not depend on which detector is looking; detectors are then
    projected one at a time by :meth:`derivatives`. That ordering matters for
    memory: a second-derivative array is ``n_parameters**2`` frequency series, so
    holding one for every detector at once reaches gigabytes at a high sampling
    frequency, while the cached polarizations are only ``2`` series per stencil
    point.

    For eleven parameters this is ``4 * n_numerical + 1`` waveform evaluations at
    ``order = 1``, plus four per numerical pair at ``order = 2``.

    Attributes
    ----------
    names : list of str
        Parameter names, in the order used by the Fisher matrix.
    analytic, numerical : list of str
        The split of ``names`` by differentiation method.
    steps : dict
        Finite-difference step actually used for each numerical parameter.
    waveform_evaluations : int
        How many times the waveform generator was called.
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
            Any generator; the waveform is treated as a black box.
        parameters : dict
            Fiducial source parameters, containing every requested name.
        fisher_parameters : sequence of str or None
            Parameters to differentiate with respect to. Defaults to
            :data:`DEFAULT_PARAMETERS`.
        earth_rotation, finite_size : bool
            Passed through to :func:`GWForge.ifo.antenna.detector_response`.
        step_sizes : dict or None
            Overrides for :data:`DEFAULT_STEP_SIZES`. These are starting guesses
            unless ``calibrate`` is False, in which case they are used as given.
        order : int
            ``1`` for first derivatives only, ``2`` to also build second
            derivatives.
        calibrate : bool
            Rescale each numerical step per source with :func:`calibrate_step`.
            On by default because a fixed step cannot serve a whole population:
            the same fractional step that is comfortably on the stable plateau
            for a 28 solar-mass binary is far below the noise floor for a
            neutron star. Costs up to two extra waveform evaluations per
            numerical parameter.
        target_change : float
            Target fractional waveform change for the calibration.
        estimate_noise : bool
            Measure the waveform model's own numerical noise with
            :func:`estimate_noise` and raise the step if it is too small to
            clear it. Costs about :data:`NOISE_SAMPLE_POINTS` extra waveform
            evaluations per numerical parameter, and never lowers a step.

        Raises
        ------
        KeyError
            If a requested parameter is missing from ``parameters``.
        """
        self.interferometers = interferometers
        self.waveform_generator = waveform_generator
        self.parameters = dict(parameters)
        self.names = list(fisher_parameters or DEFAULT_PARAMETERS)
        self.earth_rotation = earth_rotation
        self.finite_size = finite_size
        self.order = order
        self.estimate_noise = estimate_noise

        missing = [name for name in self.names if name not in self.parameters]
        if missing:
            raise KeyError(
                "Fisher parameters missing from the injection: {}".format(
                    ", ".join(sorted(missing))
                )
            )

        self.index = {name: position for position, name in enumerate(self.names)}
        self.analytic = [name for name in self.names if name in ANALYTIC_PARAMETERS]
        self.numerical = [name for name in self.names if name not in ANALYTIC_PARAMETERS]
        self.steps = {
            name: step_size(name, self.parameters[name], step_sizes)
            for name in self.numerical
        }

        #: Approximate relative error of each numerical derivative, caused by
        #: the waveform model's own numerical noise. Zero means the model is
        #: smooth at the scale of its step; a value approaching one means that
        #: derivative -- and so that row of the Fisher matrix -- is not
        #: determined. See :func:`estimate_noise`.
        self.noise_ratios = dict.fromkeys(self.numerical, 0.0)
        self.waveform_evaluations = 0
        self._cache = {}
        if calibrate:
            self._calibrate(target_change)
        self._fill_cache()

    def _calibrate(self, target):
        """Choose every numerical step, then check the model can support it."""
        _, base = self._point({})
        reference = numpy.concatenate([base["plus"], base["cross"]])
        norm = numpy.linalg.norm(reference)
        if norm == 0.0:
            raise ValueError(
                "The waveform is identically zero at the fiducial parameters; "
                "check the injection and the minimum frequency."
            )

        for name in self.numerical:
            def change(step, name=name):
                _, shifted = self._point({name: step})
                difference = (
                    numpy.concatenate([shifted["plus"], shifted["cross"]]) - reference
                )
                return numpy.linalg.norm(difference) / norm

            self.steps[name], self.noise_ratios[name] = calibrate_step(
                change, name, self.parameters[name], self.steps[name], target=target
            )
            if self.estimate_noise:
                self._apply_noise_estimate(name, base)
        logging.debug(
            "Calibrated finite-difference steps: %s",
            {name: "{:.3e}".format(step) for name, step in self.steps.items()},
        )

    def _overlap(self, name, offset, reference, interferometer):
        r"""Noise-weighted overlap :math:`\langle h(x + \delta) \mid h_0\rangle`.

        The scalar functional the noise estimator runs on. It has to be *smooth*
        in the parameter, which rules out the ``||dh|| / ||h||`` used for the
        fractional-change target -- that has a kink at zero. An inner product
        with the fiducial waveform is smooth, and is weighted by the PSD exactly
        the way the Fisher matrix is, so the noise it reports is the noise that
        actually reaches the answer.

        Kept complex on purpose. The real part is stationary at the fiducial
        point for a phase-like parameter -- for ``phase`` the overlap goes as
        ``cos(2 dphi)`` -- so a real functional alone would report a vanishing
        slope there and a meaningless noise ratio. The imaginary part is not
        stationary, and using both sidesteps the issue for any parameter.
        """
        _, shifted = self._point({name: offset})
        return complex(
            inner_product(shifted["plus"], reference["plus"], interferometer)
            + inner_product(shifted["cross"], reference["cross"], interferometer),
            inner_product(1j * shifted["plus"], reference["plus"], interferometer)
            + inner_product(1j * shifted["cross"], reference["cross"], interferometer),
        )

    def _apply_noise_estimate(self, name, base):
        """Raise the step if the model's own noise says it is too small.

        Only ever *raises* it. For a smooth waveform the noise is at the level
        of machine precision, the required step comes out far below the one the
        fractional-change target already picked, and nothing changes -- which is
        what keeps the validated IMRPhenomXHM results intact. It binds only for
        a model whose output genuinely scatters.
        """
        interferometer = self.interferometers[0]
        step = self.steps[name]

        def real_part(offset):
            return self._overlap(name, offset, base, interferometer).real

        def imaginary_part(offset):
            return self._overlap(name, offset, base, interferometer).imag

        # Probe well inside the calibrated step: the noise has to be measured
        # where the signal is small enough not to swamp it.
        half_width = step / 100.0
        real_noise, real_status = estimate_noise(sample_around(real_part, half_width))
        imaginary_noise, imaginary_status = estimate_noise(
            sample_around(imaginary_part, half_width)
        )
        if real_status != 1 and imaginary_status != 1:
            return
        noise = numpy.hypot(real_noise, imaginary_noise)
        if noise <= 0.0:
            return

        centre = self._overlap(name, 0.0, base, interferometer)
        slope = abs(self._overlap(name, step, base, interferometer) - centre) / step
        if slope <= 0.0:
            return

        # Require the change across one step to clear the noise by the margin,
        # which is the same as asking for a derivative error of 1 / NOISE_MARGIN.
        required = bounded_step(
            name, self.parameters[name], NOISE_MARGIN * noise / slope
        )
        self.noise_ratios[name] = max(
            self.noise_ratios[name], float(noise / (slope * step))
        )
        if required > step:
            self.steps[name] = required
            self.noise_ratios[name] = float(noise / (slope * required))
            logging.debug(
                "%s: raised the step from %.3e to %.3e; measured noise %.3e "
                "leaves a derivative error of %.1f%%",
                name,
                step,
                required,
                noise,
                100 * self.noise_ratios[name],
            )
        if self.noise_ratios[name] > _NOISE_WARNING_ERROR:
            logging.warning(
                "%s: the waveform model's own numerical noise leaves this "
                "derivative with a relative error of about %.0f%%, and no step "
                "size fixes it -- clearing the noise would need a step large "
                "enough to leave the linear regime. Treat %s with caution.",
                name,
                100 * self.noise_ratios[name],
                name,
            )

    def _point(self, offsets):
        """Cache key, shifted parameters and polarizations at a stencil point."""
        key = tuple(sorted(offsets.items()))
        if key not in self._cache:
            shifted = dict(self.parameters)
            for name, delta in offsets.items():
                shifted[name] = shifted[name] + delta
            self.waveform_evaluations += 1
            self._cache[key] = (
                shifted,
                self.waveform_generator.frequency_domain_strain(dict(shifted)),
            )
        return self._cache[key]

    def _fill_cache(self):
        """Evaluate the waveform at every stencil point this order needs."""
        self._point({})
        for name in self.numerical:
            for multiple, _ in _FIRST_STENCIL:
                self._point({name: multiple * self.steps[name]})
        if self.order >= 2:
            for position, left in enumerate(self.numerical):
                for right in self.numerical[position + 1 :]:
                    for sign_left in (1, -1):
                        for sign_right in (1, -1):
                            self._point(
                                {
                                    left: sign_left * self.steps[left],
                                    right: sign_right * self.steps[right],
                                }
                            )

    def _analytic_at(self, interferometer, offsets, order):
        """Closed-form response and derivatives at one cached stencil point."""
        shifted, polarizations = self._point(offsets)
        return analytic_derivatives(
            interferometer,
            polarizations,
            shifted,
            earth_rotation=self.earth_rotation,
            finite_size=self.finite_size,
            order=order,
        )

    def derivatives(self, interferometer):
        """Response and derivative arrays for one detector.

        Parameters
        ----------
        interferometer : bilby.gw.detector.Interferometer
            Must belong to the network this object was built with.

        Returns
        -------
        tuple
            ``(response, first, second)`` where ``first`` has shape
            ``(n_parameters, n_frequencies)`` and ``second`` has shape
            ``(n_parameters, n_parameters, n_frequencies)`` -- or is ``None``
            when ``order < 2``.
        """
        index, count = self.index, len(self.names)
        response, analytic_first, analytic_second = self._analytic_at(
            interferometer, {}, self.order
        )
        length = len(response)

        # Along each numerical axis the same five points give the first
        # derivative, the diagonal second derivative and -- by differencing the
        # closed-form derivatives evaluated there -- the mixed
        # numerical-analytic second derivatives.
        axis_order = 1
        axis = {
            name: {
                multiple: self._analytic_at(
                    interferometer, {name: multiple * self.steps[name]}, axis_order
                )
                for multiple, _ in _FIRST_STENCIL
            }
            for name in self.numerical
        }

        first = numpy.zeros((count, length), dtype=complex)
        for name in self.analytic:
            first[index[name]] = analytic_first[name]
        for name in self.numerical:
            first[index[name]] = self._stencil_first(
                {multiple: axis[name][multiple][0] for multiple, _ in _FIRST_STENCIL},
                self.steps[name],
            )

        if self.order < 2:
            return response, first, None

        second = numpy.zeros((count, count, length), dtype=complex)
        for pair, values in analytic_second.items():
            if pair[0] in index and pair[1] in index:
                second[index[pair[0]], index[pair[1]]] = values
                second[index[pair[1]], index[pair[0]]] = values

        for name in self.numerical:
            step, position = self.steps[name], index[name]
            points = {multiple: axis[name][multiple][0] for multiple, _ in _FIRST_STENCIL}
            points[0] = response
            second[position, position] = self._stencil_second(points, step)

            # Mixed numerical x analytic: difference the *exact* analytic
            # derivative along the numerical axis. This keeps the term by which
            # the masses reach the response through tau(f) when earth_rotation
            # is on, which applying a frozen projection to differenced
            # polarizations would drop.
            for other in self.analytic:
                if other == "luminosity_distance":
                    continue
                mixed = self._stencil_first(
                    {
                        multiple: axis[name][multiple][1][other]
                        for multiple, _ in _FIRST_STENCIL
                    },
                    step,
                )
                second[position, index[other]] = mixed
                second[index[other], position] = mixed

        for position, left in enumerate(self.numerical):
            for right in self.numerical[position + 1 :]:
                corners = {
                    (sign_left, sign_right): self._analytic_at(
                        interferometer,
                        {
                            left: sign_left * self.steps[left],
                            right: sign_right * self.steps[right],
                        },
                        1,
                    )[0]
                    for sign_left in (1, -1)
                    for sign_right in (1, -1)
                }
                mixed = (
                    corners[(1, 1)]
                    - corners[(1, -1)]
                    - corners[(-1, 1)]
                    + corners[(-1, -1)]
                ) / (4.0 * self.steps[left] * self.steps[right])
                second[index[left], index[right]] = mixed
                second[index[right], index[left]] = mixed

        # h = g(everything else) / d_L, so the whole distance row is a rescaling
        # of the first derivatives.
        if "luminosity_distance" in index:
            position = index["luminosity_distance"]
            distance = self.parameters["luminosity_distance"]
            for name in self.numerical:
                mixed = -first[index[name]] / distance
                second[position, index[name]] = mixed
                second[index[name], position] = mixed

        return response, first, second

    @staticmethod
    def _stencil_first(points, step):
        """Fourth-order first derivative from the five-point stencil."""
        total = sum(weight * points[multiple] for multiple, weight in _FIRST_STENCIL)
        return total / (_FIRST_DENOMINATOR * step)

    @staticmethod
    def _stencil_second(points, step):
        """Fourth-order second derivative from the same five points."""
        total = sum(weight * points[multiple] for multiple, weight in _SECOND_STENCIL)
        return total / (_SECOND_DENOMINATOR * step**2)
