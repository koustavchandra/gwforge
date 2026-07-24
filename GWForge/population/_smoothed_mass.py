"""Pure-numpy reimplementations of the ``gwpopulation`` smoothed mass models.

GWForge only ever used a thin slice of ``gwpopulation`` for its mass models:
the ``SinglePeak``/``MultiPeak``/``BrokenPowerLaw`` ``SmoothedMassDistribution``
classes, and only their ``m1s``/``qs`` grids plus ``p_m1``/``p_q`` evaluations.
These are the standard Talbot & Thrane (2018) analytic forms with a Planck-taper
low-mass smoothing window (see the GWTC-3 population paper, arXiv:2111.03634).

Reimplementing them here removes the runtime dependency on ``gwpopulation`` (and
its heavy transitive stack). The classes intentionally mirror the ``gwpopulation``
API (same constructor signature, ``m1s``/``qs`` attributes, ``p_m1``/``p_q``
methods, same keyword names) so they are drop-in replacements and can be validated
directly against ``gwpopulation`` -- see ``tests/test_population_smoothed_mass.py``.

The functional forms match ``gwpopulation`` exactly:

* ``powerlaw`` / ``truncnorm`` -- normalised power-law and truncated-normal pdfs.
* ``smoothing`` -- one-sided Planck-taper window on ``(mmin, mmin + delta_m]``
  (Talbot & Thrane 2018, Eqs. 7-8; note the sign fix documented there).
* ``two_component_single`` / ``three_component_single`` /
  ``double_power_law_primary_mass`` -- the primary-mass models.
* ``p_m1`` normalises by the trapezoidal integral over the ``m1s`` grid;
  ``p_q`` normalises via linear interpolation of the per-``m1`` mass-ratio norm.
"""

import numpy
from scipy.special import expit, ndtr


def powerlaw(xx, alpha, high, low):
    """Normalised power-law pdf ``p(x) ∝ x**alpha`` on ``[low, high]``.

    Matches ``gwpopulation.utils.powerlaw``. ``alpha`` is the spectral index
    (``alpha == -1`` uses the logarithmic normalisation).
    """
    xx = numpy.asarray(xx, dtype=float)
    with numpy.errstate(divide="ignore", invalid="ignore"):
        norm = numpy.where(
            numpy.asarray(alpha) == -1,
            1.0 / numpy.log(high / low),
            (1.0 + alpha)
            / (
                numpy.asarray(high) ** (1.0 + alpha)
                - numpy.asarray(low) ** (1.0 + alpha)
            ),
        )
        prob = numpy.power(xx, alpha)
    prob = prob * norm
    prob = prob * ((xx <= high) & (xx >= low))
    return prob


def truncnorm(xx, mu, sigma, high, low):
    """Truncated-normal pdf on ``[low, high]``.

    Matches ``gwpopulation.utils.truncnorm`` (which reproduces scipy's truncated
    normal normalisation via the error function).
    """
    if numpy.any(numpy.asarray(sigma) <= 0):
        raise ValueError("sigma must be greater than 0")
    xx = numpy.asarray(xx, dtype=float)
    zz = (xx - mu) / sigma
    aa = (numpy.asarray(low) - mu) / sigma
    bb = (numpy.asarray(high) - mu) / sigma
    # Normalisation is the Gaussian mass between the truncation bounds.
    norm = ndtr(bb) - ndtr(aa)
    prob = numpy.exp(-(zz**2) / 2.0) / numpy.sqrt(2.0 * numpy.pi) / sigma
    prob = prob / norm
    prob = numpy.nan_to_num(prob) * (xx >= low) * (xx <= high)
    return prob


def smoothing(masses, mmin, mmax, delta_m):
    """One-sided Planck-taper window applied to a mass pdf.

    ``S`` rises smoothly from 0 to 1 over ``(mmin, mmin + delta_m]`` and is a hard
    step at ``mmax``. Reproduces ``gwpopulation``'s ``BaseSmoothedMassDistribution
    .smoothing`` (Talbot & Thrane 2018, Eqs. 7-8). ``delta_m == 0`` -> no taper.
    """
    masses = numpy.asarray(masses, dtype=float)
    if numpy.all(delta_m == 0.0):
        return numpy.ones(numpy.shape(masses))
    shifted_mass = numpy.nan_to_num((masses - mmin) / delta_m, nan=0)
    shifted_mass = numpy.clip(shifted_mass, 1e-6, 1 - 1e-6)
    exponent = 1.0 / shifted_mass - 1.0 / (1.0 - shifted_mass)
    window = expit(-exponent)
    window *= (masses >= mmin) * (masses <= mmax)
    return window


def two_component_single(
    mass, alpha, mmin, mmax, lam, mpp, sigpp, gaussian_mass_maximum=100
):
    """Power law + single Gaussian peak primary-mass model (unsmoothed)."""
    p_pow = powerlaw(mass, alpha=-alpha, high=mmax, low=mmin)
    p_norm = truncnorm(mass, mu=mpp, sigma=sigpp, high=gaussian_mass_maximum, low=mmin)
    return (1 - lam) * p_pow + lam * p_norm


def three_component_single(
    mass,
    alpha,
    mmin,
    mmax,
    lam,
    lam_1,
    mpp_1,
    sigpp_1,
    mpp_2,
    sigpp_2,
    gaussian_mass_maximum=100,
):
    """Power law + two Gaussian peaks primary-mass model (unsmoothed)."""
    p_pow = powerlaw(mass, alpha=-alpha, high=mmax, low=mmin)
    p_norm1 = truncnorm(
        mass, mu=mpp_1, sigma=sigpp_1, high=gaussian_mass_maximum, low=mmin
    )
    p_norm2 = truncnorm(
        mass, mu=mpp_2, sigma=sigpp_2, high=gaussian_mass_maximum, low=mmin
    )
    return (1 - lam) * p_pow + lam * lam_1 * p_norm1 + lam * (1 - lam_1) * p_norm2


def double_power_law_primary_mass(mass, alpha_1, alpha_2, mmin, mmax, break_fraction):
    """Broken power-law primary-mass model (unsmoothed)."""
    m_break = mmin + break_fraction * (mmax - mmin)
    correction = powerlaw(m_break, alpha=-alpha_2, low=m_break, high=mmax) / powerlaw(
        m_break, alpha=-alpha_1, low=mmin, high=m_break
    )
    low_part = powerlaw(mass, alpha=-alpha_1, low=mmin, high=m_break)
    high_part = powerlaw(mass, alpha=-alpha_2, low=m_break, high=mmax)
    prob = low_part * (mass < m_break) * correction + high_part * (mass >= m_break)
    return prob / (1 + correction)


def left_truncated_normal(mass, mu, sigma, low):
    """Normal pdf truncated below at ``low`` and unbounded above (N_lt).

    Used by the BGP mass model (GWTC-5.0 Eq. B12). Implemented as the truncated
    normal with an infinite upper bound.
    """
    return truncnorm(mass, mu=mu, sigma=sigma, high=numpy.inf, low=low)


def broken_power_law(mass, alpha_1, alpha_2, m_break, mmin, m_high):
    r"""Broken power law of GWTC-5.0 Eq. B10, normalised via Eq. B11.

    .. math::
        p_{BP}(m) = \frac{1}{N} \begin{cases}
            (m/m_{break})^{-\alpha_1} & m_{min} \le m < m_{break} \\
            (m/m_{break})^{-\alpha_2} & m_{break} \le m < m_{high}
        \end{cases}

    with normalisation :math:`N = m_{break}\left[
    \frac{1 - (m_{min}/m_{break})^{1-\alpha_1}}{1-\alpha_1}
    + \frac{(m_{high}/m_{break})^{1-\alpha_2} - 1}{1-\alpha_2}\right]`.
    """
    mass = numpy.asarray(mass, dtype=float)
    norm = m_break * (
        (1 - (mmin / m_break) ** (1 - alpha_1)) / (1 - alpha_1)
        + ((m_high / m_break) ** (1 - alpha_2) - 1) / (1 - alpha_2)
    )
    low_part = (mass / m_break) ** (-alpha_1) * (mass >= mmin) * (mass < m_break)
    high_part = (mass / m_break) ** (-alpha_2) * (mass >= m_break) * (mass < m_high)
    return (low_part + high_part) / norm


def broken_power_law_two_peak(
    mass,
    alpha_1,
    alpha_2,
    m_break,
    mmin,
    m_high,
    lam_0,
    lam_1,
    mpp_1,
    sigpp_1,
    mpp_2,
    sigpp_2,
):
    r"""Broken power law + two left-truncated Gaussian peaks (GWTC-5.0 Eq. B12).

    The mixture (before the low-mass Planck taper, which the base class applies):

    .. math::
        \pi(m) \propto \lambda_0\, p_{BP}(m)
                     + \lambda_1\, N_{lt}(m | \mu_1, \sigma_1)
                     + (1 - \lambda_0 - \lambda_1)\, N_{lt}(m | \mu_2, \sigma_2).

    This is the "BGP" (Broken power law + Gaussian Peaks) fiducial BBH mass model
    used in the GWTC-4.0/5.0 population analyses.
    """
    pbp = broken_power_law(mass, alpha_1, alpha_2, m_break, mmin, m_high)
    peak_1 = left_truncated_normal(mass, mpp_1, sigpp_1, low=mmin)
    peak_2 = left_truncated_normal(mass, mpp_2, sigpp_2, low=mmin)
    return lam_0 * pbp + lam_1 * peak_1 + (1 - lam_0 - lam_1) * peak_2


class BaseSmoothedMassDistribution:
    """Shared machinery for the smoothed primary-mass + power-law mass-ratio models.

    Mirrors ``gwpopulation.models.mass.BaseSmoothedMassDistribution`` closely enough
    to be a drop-in replacement for GWForge's usage: it exposes ``m1s``/``qs`` grids
    and ``p_m1``/``p_q`` with identical keyword names and normalisation conventions.
    """

    #: set by subclasses to one of the module-level primary-mass functions
    primary_model = None

    def __init__(self, mmin=2, mmax=100, normalization_shape=(1000, 500)):
        self.mmin = mmin
        self.mmax = mmax
        self.m1s = numpy.linspace(mmin, mmax, normalization_shape[0])
        self.qs = numpy.linspace(0.001, 1, normalization_shape[1])
        self.m1s_grid, self.qs_grid = numpy.meshgrid(self.m1s, self.qs)

    def p_m1(self, dataset, **kwargs):
        mmin = kwargs.get("mmin", self.mmin)
        delta_m = kwargs.pop("delta_m", 0)
        p_m = self.primary_model(dataset["mass_1"], **kwargs)
        p_m = p_m * smoothing(
            dataset["mass_1"], mmin=mmin, mmax=self.mmax, delta_m=delta_m
        )
        return p_m / self.norm_p_m1(delta_m=delta_m, **kwargs)

    def norm_p_m1(self, delta_m, **kwargs):
        mmin = kwargs.get("mmin", self.mmin)
        if delta_m == 0:
            return 1
        p_m = self.primary_model(self.m1s, **kwargs)
        p_m = p_m * smoothing(self.m1s, mmin=mmin, mmax=self.mmax, delta_m=delta_m)
        return numpy.nan_to_num(numpy.trapezoid(p_m, self.m1s))

    def p_q(self, dataset, beta, mmin, delta_m):
        p_q = powerlaw(dataset["mass_ratio"], beta, 1, mmin / dataset["mass_1"])
        p_q = p_q * smoothing(
            dataset["mass_1"] * dataset["mass_ratio"],
            mmin=mmin,
            mmax=dataset["mass_1"],
            delta_m=delta_m,
        )
        # The per-m1 norm vanishes as m1 -> mmin (the q-support collapses); the
        # resulting 0/0 is expected and cleaned up by nan_to_num below.
        with numpy.errstate(divide="ignore", invalid="ignore"):
            p_q = p_q / self._norm_p_q(
                beta=beta, mmin=mmin, delta_m=delta_m, masses=dataset["mass_1"]
            )
        return numpy.nan_to_num(p_q)

    def _norm_p_q(self, beta, mmin, delta_m, masses):
        """Per-``m1`` mass-ratio norm, linearly interpolated onto ``masses``."""
        p_q = powerlaw(self.qs_grid, beta, 1, mmin / self.m1s_grid)
        p_q = p_q * smoothing(
            self.m1s_grid * self.qs_grid, mmin=mmin, mmax=self.m1s_grid, delta_m=delta_m
        )
        norms = numpy.nan_to_num(numpy.trapezoid(p_q, self.qs, axis=0))
        return numpy.interp(masses, self.m1s, norms)


class SinglePeakSmoothedMassDistribution(BaseSmoothedMassDistribution):
    primary_model = staticmethod(two_component_single)


class MultiPeakSmoothedMassDistribution(BaseSmoothedMassDistribution):
    primary_model = staticmethod(three_component_single)


class BrokenPowerLawSmoothedMassDistribution(BaseSmoothedMassDistribution):
    primary_model = staticmethod(double_power_law_primary_mass)


class BrokenPowerLawTwoPeakSmoothedMassDistribution(BaseSmoothedMassDistribution):
    """BGP model: broken power law + two peaks, with low-mass Planck taper and a
    power-law mass ratio (GWTC-5.0 App. B.4, Eqs. B10-B14)."""

    primary_model = staticmethod(broken_power_law_two_peak)
