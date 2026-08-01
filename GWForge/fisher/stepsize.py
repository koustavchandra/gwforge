r"""Choosing finite-difference steps, and measuring what limits them.

Two things bound a usable step from opposite sides. Too large and the difference
stops being local; too small and it drowns in the waveform model's own numerical
noise. The constants here are measured rather than guessed, and the measurements
are recorded beside them because they are the whole justification.

:func:`calibrate_step` handles the upper bound by targeting a fixed fractional
change in the waveform. :func:`estimate_noise` handles the lower one, using the
ECNoise procedure of Moré and Wild, *Estimating Computational Noise*, SIAM J.
Sci. Comput. 33(3):1292--1314 (2011): a smooth function's k-th finite differences
shrink like :math:`h^k` while noise does not shrink at all, so the level at which
the difference table stops changing is the noise.
"""

import logging

import numpy

from .parameters import bounded_step

#: Target for :func:`calibrate_step`: the fractional change
#: :math:`\\lVert\\Delta h\\rVert/\\lVert h\\rVert` a single step should produce.
#:
#: This is the single most important accuracy knob in the module, and 1e-2 is
#: measured rather than guessed. Below roughly 1e-3 the Fisher matrix for
#: ``IMRPhenomXHM`` starts to drift and then diverges as the step shrinks --
#: sigma(eta) moves by 150% over a decade -- because LAL places each ``(l, m)``
#: mode's onset at a mass-dependent frequency, so differencing across that edge
#: contributes a spurious term growing like ``1/step``. Above 1e-3 the Fisher is
#: flat to a fraction of a percent for every parameter, and truncation is
#: harmless: a fractional change of 1e-2 means a phase change of ~1e-2 rad per
#: step, so the fourth-order stencil truncates at ``(1e-2)**4 / 30 ~ 3e-10``.
TARGET_FRACTIONAL_CHANGE = 1e-2


#: How far the change across one step must exceed the model's measured noise.
#: The derivative's relative error is roughly the reciprocal, so 100 buys ~1%.
NOISE_MARGIN = 100.0

#: Derivative relative error above which a parameter is worth warning about.
_NOISE_WARNING_ERROR = 0.05

#: Acceptance window and iteration budget for the fractional-change search.
#:
#: The window used to be a factor of two either side, which let the search stop
#: almost anywhere within a factor of four and made the chosen step depend on
#: where it started. That mattered: for SEOBNRv5HM the usable plateau in
#: ``chi_1`` begins around a step of 3e-4, and a loose window let some starting
#: guesses settle just below it, moving ``sigma(chi_1)`` by a factor of four.
#: The fractional change is very nearly linear in the step, so tightening the
#: window costs at most an iteration or two.
_CALIBRATION_TOLERANCE_LOW = 0.8
_CALIBRATION_TOLERANCE_HIGH = 1.25
_CALIBRATION_ITERATIONS = 5

#: Points used by :func:`estimate_noise`. Moré and Wild find ``m = 6`` (seven
#: points) enough for a stochastic function but recommend ``m = 8`` for a
#: deterministic one whose noise comes from adaptive internal tolerances --
#: which is exactly what a LAL or EOB waveform is.
NOISE_SAMPLE_POINTS = 9

#: Successive noise estimates must agree within this factor before one is
#: accepted. From the reference implementation; it also sets the accuracy one
#: can expect, so two independent estimates may legitimately differ by 16.
_NOISE_RATIO_TOLERANCE = 4.0


def estimate_noise(values):
    r"""Estimate the numerical noise of a computed function from its values.

    Implements the ECNoise procedure of Moré and Wild, *Estimating Computational
    Noise*, SIAM J. Sci. Comput. 33(3):1292--1314 (2011), which is the standard
    way to measure how much of a simulated function's output is numerical
    scatter rather than signal.

    The idea is that for a smooth :math:`f` the :math:`k`-th finite difference
    over a spacing :math:`h` shrinks like :math:`h^k`, while iid noise does not
    shrink at all -- it merely gets amplified by a known binomial factor. So
    build the difference table and look for the level at which the estimates
    stop changing:

    .. math::

       \sigma_k^2 = \frac{\gamma_k}{m + 1 - k} \sum_i \big(\Delta^k f_i\big)^2,
       \qquad \gamma_k = \frac{(k!)^2}{(2k)!}.

    Parameters
    ----------
    values : array_like
        Function values at equally spaced points. At least four, preferably
        :data:`NOISE_SAMPLE_POINTS`. The spacing itself never enters.

    Returns
    -------
    tuple
        ``(noise, status)``. ``status`` is ``1`` when the noise was detected,
        ``2`` when the spacing was too *small* (half the first differences came
        out exactly zero -- widen it by ~100), and ``3`` when it was too *large*
        (no level settled, or the values already differ in their first digit --
        narrow it by ~100). ``noise`` is ``0.0`` unless ``status`` is ``1``.
    """
    values = numpy.asarray(values, dtype=float)
    count = len(values)
    if count < 4:
        raise ValueError("estimate_noise needs at least four values")

    smallest, largest = values.min(), values.max()
    scale = max(abs(largest), abs(smallest))
    if scale > 0 and (largest - smallest) / scale > 0.1:
        # The values differ in the first digit, so the signal swamps any noise.
        return 0.0, 3

    table = values.copy()
    levels = numpy.zeros(count - 1)
    changes_sign = numpy.zeros(count - 1, dtype=bool)
    gamma = 1.0
    for order in range(1, count):
        table = numpy.diff(table)
        if order == 1 and numpy.count_nonzero(table == 0) >= count / 2:
            return 0.0, 2
        gamma *= 0.5 * order / (2 * order - 1)
        levels[order - 1] = numpy.sqrt(gamma * numpy.mean(table**2))
        changes_sign[order - 1] = table.min() * table.max() < 0.0

    for index in range(count - 3):
        window = levels[index : index + 3]
        if window.max() <= _NOISE_RATIO_TOLERANCE * window.min() and changes_sign[index]:
            return float(levels[index]), 1
    return 0.0, 3


def sample_around(evaluate, half_width, points=NOISE_SAMPLE_POINTS):
    """Evaluate a function at equally spaced points centred on zero.

    The spacing convention follows Moré and Wild's driver: ``half_width`` is the
    half-width of the span, not the gap between points.

    Parameters
    ----------
    evaluate : callable
        ``evaluate(offset) -> float``.
    half_width : float
        Half the total span.
    points : int
        Number of evaluations.

    Returns
    -------
    numpy.ndarray
    """
    middle = (points + 1) // 2 - 1
    offsets = [2.0 * (index - middle) / (points - 1) * half_width for index in range(points)]
    return numpy.array([evaluate(offset) for offset in offsets])


def calibrate_step(
    evaluate,
    name,
    value,
    initial,
    target=TARGET_FRACTIONAL_CHANGE,
):
    """Rescale a step so one displacement changes the waveform by ``target``.

    The fractional change is linear in the step over many decades, so a single
    probe fixes the scale; a second probe confirms it and catches the
    non-linear case. See :data:`TARGET_FRACTIONAL_CHANGE` for why the target is
    what it is.

    Parameters
    ----------
    evaluate : callable
        ``evaluate(step) -> float`` returning
        :math:`\\lVert\\Delta h\\rVert / \\lVert h\\rVert` for that displacement.
    name : str
        Parameter name, used for the bound check and for logging.
    value : float
        Fiducial value.
    initial : float
        Starting step.
    target : float
        Desired fractional change.

    Returns
    -------
    tuple
        ``(step, noise_ratio)``. The step is bounded. The noise ratio is always
        ``0.0`` here -- measuring it is
        :meth:`WaveformDerivatives._apply_noise_estimate`'s job, which uses
        :func:`estimate_noise`. Falls back to ``initial`` if the waveform does
        not depend on the parameter at all.
    """
    step = initial
    change = None
    for _ in range(_CALIBRATION_ITERATIONS):
        change = evaluate(step)
        if not numpy.isfinite(change) or change <= 0.0:
            logging.warning(
                "The waveform does not respond to %s at a step of %.3e; keeping "
                "the default step. The Fisher matrix will be singular in this "
                "parameter.",
                name,
                step,
            )
            return step, 0.0
        if _CALIBRATION_TOLERANCE_LOW * target <= change <= _CALIBRATION_TOLERANCE_HIGH * target:
            break
        step = bounded_step(name, value, step * target / change)

    return step, 0.0
