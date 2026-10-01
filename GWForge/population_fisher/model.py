r"""The population models GWForge forecasts on, with their analytic scores.

Four sectors live here -- the BGP mass function, the Default BBH spin
distribution, Madau-Dickinson redshift evolution, and the detector-frame
spectral-siren composite that folds a cosmology into the first three. They share
one module because they share one contract
(:class:`GWForge.population_fisher.base.PopulationModel`) and because splitting
them by sector bought nothing: every consumer wants the set, not one of them.

The densities are **not** reimplemented here. The mass model calls
:mod:`GWForge.population._smoothed_mass`, the redshift model calls
:func:`GWForge.population.redshift.madau_dickinson_psi_of_z`, and the spin model
calls :func:`GWForge.population.spin.truncated_normal` -- the kernel both halves
of the ``Default`` model are built from, and the one the density checks in
``validation/`` use. So a forecast and the catalogue it runs on cannot describe
two different populations. What is new in every case is the *score*,
:math:`\partial_\Lambda \ln p`, in closed form.

The sections below are the four sectors' own documentation, unchanged from when
they were separate modules.

Mass: the BGP model
===================

The BGP mass model -- broken power law plus two peaks -- with analytic scores.

This is the fiducial BBH mass model of the GWTC-4.0/5.0 population analyses
(`arXiv:2605.27226 <https://arxiv.org/abs/2605.27226>`_ Eqs. B10-B14), and the
one :class:`GWForge.population.mass.Mass` samples under ``mass-model = BGP``. The
density is *not* reimplemented here: it is
:mod:`GWForge.population._smoothed_mass`, called directly, so a Fisher forecast
and the catalogue it is run on cannot describe two different populations. What
is new is the score.

.. math::

   \pi(m_1) = \frac{1}{Z}\Big[\lambda_0\, p_{\rm BP}(m_1)
        + \lambda_1\, N_{\rm lt}(m_1|\mu_1,\sigma_1)
        + (1 - \lambda_0 - \lambda_1)\, N_{\rm lt}(m_1|\mu_2,\sigma_2)\Big]
        S(m_1|m_{\min},\delta_m),

.. math::

   p(m_2 \mid m_1)
     = \frac{m_2^{\beta}\, S(m_2 \mid m_{2,\rm low}, \delta_{m,2})}{G(m_1)},
   \qquad
   G(m_1) = \int_{m_{2,\rm low}}^{m_1} m_2^{\beta}\, S(m_2)\, dm_2.

The primary and secondary tapers are **independent** (GWTC-5.0 Tab. 5): the
secondary's prior is :math:`m_{2,\rm low} \sim U(3, m_{1,\rm low})`, so it turns
on below the primary, and the O4b medians separate them clearly -- 4.49/3.12
against 3.49/5.60. Both secondary parameters are constructor arguments here, for
the same reason :math:`m_{\min}` is: they move the edge of the support.

``log_prob`` returns the joint density in :math:`(m_1, m_2)`, not in
:math:`(m_1, q)`. The two differ by the Jacobian :math:`m_1`, which carries no
hyper-parameter and so cancels from every score -- but it does not cancel from
``log_prob``, which the maximum-likelihood fits compare, and the detector-frame
model in :mod:`GWForge.population_fisher.model` wants component masses
anyway because that is where the :math:`(1+z)` factors live. Events may be given
as ``mass_2_source`` or as ``mass_ratio``; both describe the same density.

Which parameters are free, and why the rest are not
---------------------------------------------------

:attr:`BrokenPowerLawTwoPeakMass.parameter_names` holds eleven names. The three
that are missing -- ``mmin``, ``m_high`` and ``maximum_mass`` -- are
**constructor arguments**, not hyper-parameters, for two independent reasons.

*Statistically*, they are hard cutoffs. Inside the power-law interior their
scores are the same number for every event, so the centred score
:math:`s_k - \langle s \rangle` of the Term-I Fisher vanishes identically and
the matrix acquires a zero row. The previous generation of this code discovered
this the hard way and pinned exactly these parameters in every science run.

*Numerically*, they are the endpoints of the normalisation integrals, so their
derivatives carry Leibniz boundary terms and cannot be checked against a finite
difference on a fixed grid -- as ``m_high`` crosses a node, the trapezoid jumps
by the integrand times the node spacing, so the finite difference is a staircase
rather than a derivative. Every score this module *does* expose is smooth in its
parameter and is validated against a finite difference in the test suite.

To study sensitivity to them, build a second model with different values. That
is honest about what is happening -- a different support is a different model,
not a different point in one.

``delta_m`` is free and is smooth, because the Planck taper drives the integrand
continuously to zero at :math:`m_{\min}`: there is no edge for a boundary term
to sit on. It is also the parameter that most deserves watching. The taper is a
sharp feature at a *fixed source-frame* mass, i.e. a standard scale, so it feeds
:math:`H_0` directly in a spectral-siren forecast. The previous generation of
this code measured :math:`\sigma(H_0)` swinging by a factor of 77 across
:math:`\delta_m \in [0, 4]\,M_\odot`; that sensitivity should be scanned rather
than assumed.

Grids
-----

One fixed grid over :math:`[m_{\min}, m_{\rm max, grid}]` serves both
normalisations: the primary-mass :math:`Z` and the cumulative mass-ratio norm
:math:`G`.

*Fixed*, because a grid that moved with the hyper-parameters would make
``d/dLambda`` of the discretised integral differ from the discretised
``d/dLambda``, and every score here is the latter.

*Clustered towards* :math:`m_{\min}` rather than uniform -- see
:func:`mass_grid` -- because :math:`G` is interpolated at each event's
:math:`m_1`, and just above :math:`m_{\min}` it is a very small number rising
very fast.

Spin: the Default BBH model
===========================

The GWTC-5.0 **Default BBH** spin population, with analytic scores.

This is Eqs. B15-B16 of `arXiv:2605.27226 <https://arxiv.org/abs/2605.27226>`_,
the model :class:`GWForge.population.spin.Spin` samples under
``spin-model = Default``. The sampler has no density attached to it, so unlike
the mass and redshift blocks the density here is written out rather than reused;
the tests check it against ``scipy.stats`` on the same parameters the sampler
uses.

Magnitudes are independent and identically distributed truncated Gaussians
(Eq. B15),

.. math::

   \pi(\chi_1, \chi_2 \mid \mu_\chi, \sigma_\chi)
     = N_{[0, a_{\max}]}(\chi_1 \mid \mu_\chi, \sigma_\chi)\,
       N_{[0, a_{\max}]}(\chi_2 \mid \mu_\chi, \sigma_\chi),

and the cosine tilts are identically but **not independently** distributed
(Eq. B16),

.. math::

   \pi(\cos\theta_1, \cos\theta_2 \mid \mu_t, \sigma_t, \xi)
     = \xi\, N_{[t_{\min}, 1]}(\cos\theta_1 \mid \mu_t, \sigma_t)\,
             N_{[t_{\min}, 1]}(\cos\theta_2 \mid \mu_t, \sigma_t)
     + \frac{1 - \xi}{(1 - t_{\min})^2}.

Three things that are easy to get wrong
---------------------------------------

**The magnitude is a truncated Gaussian, not a Beta**, and ``sigma_chi`` is a
*standard deviation*. The GWTC-3-era model this package used to implement takes
a Beta parameterised by a mean and a **variance**; the two are not
interconvertible, and the Beta is forced to vanish or diverge at both endpoints
where the truncated Gaussian is finite and generally nonzero.

**The tilt Gaussian's mean is free.** :math:`\mu_t \sim U(-1, 1)`, with an O4b
median of 0.279 -- the paper says outright that it "allow[s] for the location of
the Gaussian subpopulation to vary". A model with :math:`\mu_t` pinned at 1
cannot produce the peak away from alignment that the data prefer.

**The tilt mixture is over the binary, not the component.** One Bernoulli draw
decides whether *both* tilts come from the Gaussian, so the two are correlated.
Factorising it into two per-component mixtures gives identical marginals and the
wrong joint -- the difference is a genuine correlation between
:math:`\cos\theta_1` and :math:`\cos\theta_2` that vanishes only at
:math:`\xi \in \{0, 1\}`. Because this leaves the marginals untouched, only a
joint statistic catches it; the test suite uses the correlation coefficient.

Frames and Jacobians
--------------------

The sampler returns ``tilt_1``/``tilt_2``, not their cosines. The Jacobian
:math:`|\sin t|` relating the two carries no hyper-parameter, so it cancels from
every score; ``log_prob`` is quoted in :math:`\cos t` and the model accepts
either ``cos_tilt_i`` or ``tilt_i``.

Redshift: Madau-Dickinson
=========================

Madau-Dickinson merger-rate density, with analytic scores.

The observed rate per unit redshift is

.. math::

   \frac{dN}{dz} \propto \psi(z)\,\frac{dV_c}{dz}\,\frac{1}{1+z},

so, normalised over :math:`[0, z_{\max}]`,

.. math::

   \ln p(z \mid \gamma, \kappa, z_p) = \ln\psi(z) + \ln\frac{dV_c}{dz}
       - \ln(1 + z) - \ln Z(\gamma, \kappa, z_p),

with :math:`\psi` the Madau-Dickinson shape that
:func:`GWForge.population.redshift.madau_dickinson_psi_of_z` already implements
and this module calls directly:

.. math::

   \psi(z) = \frac{(1+z)^\gamma}{1 + \left(\frac{1+z}{1+z_p}\right)^{\kappa}}
             \left[1 + (1 + z_p)^{-\kappa}\right].

Two simplifications fall out and are used below.

* The :math:`dV_c/dz\,/(1+z)` measure carries no hyper-parameter, so it drops
  out of every score exactly. It is *not* dropped from ``log_prob``, which the
  maximum-likelihood fit in :mod:`GWForge.population_fisher.fit` needs whole.
* The :math:`[1 + (1+z_p)^{-\kappa}]` factor -- there to make
  :math:`\psi(0) = 1` -- is independent of :math:`z`, so it multiplies the
  numerator and :math:`Z` alike and cancels from ``log_prob`` identically. The
  scores below are therefore of the *shape* only.

What is left is the same structure everywhere:

.. math::

   \frac{\partial \ln p}{\partial\Lambda} = g(z) - \mathbb{E}_p[g],

where :math:`\mathbb{E}_p` is the expectation under :math:`p(z)` itself, coming
from :math:`\partial_\Lambda \ln Z`. It is evaluated on the same trapezoidal
grid that normalises ``log_prob``, so the density and its derivative cannot
drift apart.

Note on fiducials
-----------------

``GWForge.population.redshift.Redshift`` always convolves :math:`\psi` with a
formation-to-merger time-delay distribution, so a population it generates is
*not* distributed as :math:`\psi(z)\,dV_c/dz/(1+z)`. The fiducial
:math:`(\gamma, \kappa, z_p)` for a Fisher forecast on such a catalogue must
therefore be **fitted** to the injected redshifts rather than read off the
generation config -- see :func:`GWForge.population_fisher.fit.fit_model`.

Spectral sirens
===============

Spectral sirens: the population in the detector frame, with cosmology free.

The detector measures redshifted masses, so a feature at a fixed *source-frame*
mass -- the low-mass peak, the taper at :math:`m_{\min}` -- appears at a
detector-frame mass that grows with redshift. Combined with the measured
:math:`d_L`, that pins the distance-redshift relation and hence
:math:`(H_0, \Omega_{m,0}, w_0)`.

Why the detector frame
----------------------

This is the one structural choice that has to be right. Promoting the cosmology
to hyper-parameters and building the density over the **observables**
:math:`(m_1^{\rm det}, m_2^{\rm det}, d_L)` means :math:`P_{\rm det}` depends on
no hyper-parameter -- it is a threshold on an SNR computed from detector-frame
quantities alone. The centred-score identity
:math:`\partial_\Lambda \ln\alpha = \mathbb{E}_{\rm det}[s]` that supplies the
selection term in :mod:`GWForge.population_fisher.population_fisher_term_I` then holds *exactly*.

Push the events to the source frame at a trial cosmology instead and both break:
the detection probability acquires a cosmology dependence, and the per-event
measurement uncertainties become functions of the hyper-parameters.

The density
-----------

The source-frame model gives :math:`p(m_1^{\rm src}, m_2^{\rm src})` and
:math:`p(z)`, both already normalised. Changing variables to the observables
costs one Jacobian factor per mass and one for the distance,

.. math::

   \ln p(m_1^{\rm det}, m_2^{\rm det}, d_L)
     = \ln p_{\rm mass}(m_1^{\rm src}, m_2^{\rm src})
     + \ln p(z)
     - 2\ln(1 + z)
     - \ln\frac{dd_L}{dz}
     \;[+ \ln p_{\rm spin}],

with :math:`z = z(d_L; H_0, \Omega_{m,0}, w_0)` and
:math:`m^{\rm src} = m^{\rm det}/(1+z)`. Spins are frame invariant and simply
add.

The redshift measure is rebuilt at the trial cosmology rather than taken from
:class:`GWForge.population_fisher.model.MadauDickinsonRedshift`, whose
:math:`dV_c/dz` is fixed: it is precisely the cosmology dependence of that
measure that carries the information.

The cosmology columns
---------------------

Every one is closed form. Writing :math:`\Lambda` for a cosmology parameter,
everything reaches it through :math:`z(d_L)` at fixed distance, plus an explicit
dependence of the measure at fixed :math:`z`:

.. math::

   \frac{\partial \ln p}{\partial \Lambda}
     = \left[
        \frac{\partial \ln p_{\rm mass}}{\partial m_1^{\rm src}}
        \left(\frac{-m_1^{\rm src}}{1+z}\right)
      + \frac{\partial \ln p_{\rm mass}}{\partial m_2^{\rm src}}
        \left(\frac{-m_2^{\rm src}}{1+z}\right)
      + \frac{\partial \ln \psi}{\partial z}
      + \frac{\partial}{\partial z}\ln\frac{dV_c}{dz}
      - \frac{3}{1+z}
      - \frac{\partial}{\partial z}\ln\frac{dd_L}{dz}
     \right] \frac{\partial z}{\partial \Lambda}
     + \left.\frac{\partial}{\partial\Lambda}\ln\frac{dV_c}{dz}\right|_z
     - \left.\frac{\partial}{\partial\Lambda}\ln\frac{dd_L}{dz}\right|_z
     - \frac{\partial \ln Z_z}{\partial \Lambda},

with :math:`\partial z/\partial\Lambda` from implicit differentiation of
:math:`d_L(z;\Lambda) = \text{const}` (see
:meth:`GWForge.cosmology.FlatwCDM.dredshift_dparameter`) and every remaining
piece from :meth:`GWForge.cosmology.FlatwCDM.derivatives`.

Nothing here is finite-differenced. The previous generation of this code
finite-differenced the :math:`\Omega_{m,0}` and :math:`w_0` columns and hit an
error floor set by rebuilding a spline per trial cosmology, such that shrinking
the step from 1e-4 to 1e-5 made them *two orders of magnitude worse*.

A caveat that must be quoted with any result
---------------------------------------------

The taper at :math:`m_{\min}` is a sharp feature at a fixed source-frame mass,
i.e. a standard scale, and it feeds :math:`H_0` directly. How sharp it is --
``delta_m`` -- therefore sets how much cosmological information the low-mass
edge carries, and :math:`\sigma(H_0)` moves by nearly two orders of magnitude
across :math:`\delta_m \in [0, 4]\,M_\odot`. Because GWForge *generates* the
catalogue at a known ``delta_m``, that number is under control here rather than
being set by whatever a maximum-likelihood fit does against a hard injection
cutoff -- but the sensitivity is real and should be scanned, not assumed.

"""

import numpy
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicHermiteSpline
from scipy.special import ndtr

from ..cosmology import (
    COSMOLOGY_PARAMETERS,
    FlatwCDM,
    differential_comoving_volume,
)
from ..population._smoothed_mass import (
    broken_power_law,
    left_truncated_normal,
    smoothing,
)
from ..population.mass import BGP_PARAMETERS
from ..population.redshift import madau_dickinson_psi_of_z
from ..population.spin import DEFAULT_BBH_SPIN_PARAMETERS, truncated_normal
from .base import PopulationModel



# =========================================================================
# Mass: the BGP model
# =========================================================================

# Nodes for the shared mass grid, which is quadratically clustered towards
# ``mmin`` (see :func:`mass_grid`). Over the default [5, 200] support this puts
# the spacing at 1e-6 Msun at the foot of the taper and 0.02 Msun at the top,
# so both the delta_m = 4.8 ramp and a sigma = 1.5 peak are resolved by several
# hundred nodes. Measured drift in the scores between 20000 and 80000 nodes is
# below 1e-8 relative.
MASS_GRID_NODES = 20000

# Below this the ``(1 - alpha)`` denominators in the broken-power-law
# normalisation are evaluated by series rather than directly.
POWER_LAW_SERIES_THRESHOLD = 1e-8

# 1 / sqrt(2 pi).
_INVERSE_ROOT_TWO_PI = 1.0 / numpy.sqrt(2.0 * numpy.pi)


def mass_grid(mmin, maximum_mass, nodes=MASS_GRID_NODES):
    r"""The shared normalisation grid: quadratically clustered towards ``mmin``.

    :math:`m = m_{\min} + (m_{\max} - m_{\min})\,s^2` for
    :math:`s \in [0, 1]`.

    Uniform spacing is the obvious choice and the wrong one. The mass-ratio
    normalisation :math:`G(m_1) = \int_{m_{\min}}^{m_1} m_2^{\beta} S(m_2)dm_2`
    is interpolated at each event's :math:`m_1`, and just above :math:`m_{\min}`
    it is a very small number rising very fast, so an interpolant on a uniform
    grid has a large *relative* error exactly there. Measured against a finite
    difference at the same node count, ``d ln p / d m_1`` at
    :math:`m_1 = m_{\min} + 0.73` was wrong by 3.6e-1 on a uniform grid and by
    4.7e-2 on this one; the cubic-Hermite interpolation in
    :meth:`BrokenPowerLawTwoPeakMass._cumulative_secondary` then takes it to
    6.7e-5.

    The grid is fixed -- it does not depend on any hyper-parameter -- which is
    what lets ``d/dLambda`` of the discretised integral equal the discretised
    ``d/dLambda`` exactly, with no Leibniz boundary terms.

    Parameters
    ----------
    mmin, maximum_mass : float
    nodes : int

    Returns
    -------
    numpy.ndarray
    """
    fraction = numpy.linspace(0.0, 1.0, int(nodes))
    return mmin + (maximum_mass - mmin) * fraction**2


def _power_law_integral(ratio, alpha):
    r""":math:`A(x, \alpha) = \int_x^1 t^{-\alpha} dt = (1 - x^{1-\alpha})/(1-\alpha)`.

    Written with ``expm1`` so the :math:`\alpha \to 1` cancellation is exact
    rather than catastrophic, with a series branch at the removable singularity
    itself (where :math:`A \to -\ln x`).
    """
    exponent = 1.0 - alpha
    logarithm = numpy.log(ratio)
    if abs(exponent) < POWER_LAW_SERIES_THRESHOLD:
        return -logarithm * (
            1.0
            + 0.5 * exponent * logarithm
            + (exponent * logarithm) ** 2 / 6.0
        )
    return -numpy.expm1(exponent * logarithm) / exponent


def _power_law_integral_dalpha(ratio, alpha):
    r"""Derivative of :func:`_power_law_integral` with respect to the slope.

    .. math::
        \frac{\partial A}{\partial\alpha}
          = \frac{u x^{u}\ln x + (1 - x^{u})}{u^2},
        \qquad u = 1 - \alpha.

    Series branch as :math:`u \to 0`, where the limit is
    :math:`\tfrac{1}{2}\ln^2 x`.
    """
    exponent = 1.0 - alpha
    logarithm = numpy.log(ratio)
    if abs(exponent) < POWER_LAW_SERIES_THRESHOLD:
        return 0.5 * logarithm**2 + exponent * logarithm**3 / 3.0
    power = ratio**exponent
    return (exponent * power * logarithm - numpy.expm1(exponent * logarithm)) / (
        exponent**2
    )


def broken_power_law_norm(alpha_1, alpha_2, m_break, mmin, m_high):
    r"""Normalisation of :func:`GWForge.population._smoothed_mass.broken_power_law`.

    .. math::

       N = m_b\Big[A(m_{\min}/m_b, \alpha_1) - A(m_{\rm high}/m_b, \alpha_2)\Big].

    Returns
    -------
    float
    """
    return m_break * (
        _power_law_integral(mmin / m_break, alpha_1)
        - _power_law_integral(m_high / m_break, alpha_2)
    )


def broken_power_law_norm_derivatives(alpha_1, alpha_2, m_break, mmin, m_high):
    r"""``dN/dLambda`` for the broken power law, as a dict.

    Only the three interior parameters are returned; ``mmin`` and ``m_high``
    move the support and are constructor arguments here (see the module
    docstring).

    .. math::

       \frac{\partial N}{\partial m_b}
         = (A_{\rm lo} - A_{\rm hi})
           + (m_{\min}/m_b)^{1-\alpha_1} - (m_{\rm high}/m_b)^{1-\alpha_2}.

    Returns
    -------
    dict
        Keys ``"alpha_1"``, ``"alpha_2"``, ``"m_break"``.
    """
    low_ratio = mmin / m_break
    high_ratio = m_high / m_break
    low_integral = _power_law_integral(low_ratio, alpha_1)
    high_integral = _power_law_integral(high_ratio, alpha_2)
    return {
        "alpha_1": m_break * _power_law_integral_dalpha(low_ratio, alpha_1),
        "alpha_2": -m_break * _power_law_integral_dalpha(high_ratio, alpha_2),
        "m_break": (low_integral - high_integral)
        + low_ratio ** (1.0 - alpha_1)
        - high_ratio ** (1.0 - alpha_2),
    }


def left_truncated_normal_dlog(mass, mu, sigma, mmin):
    r"""``d ln N_lt / d(mu, sigma)`` for the left-truncated normal.

    With :math:`z = (m - \mu)/\sigma`, :math:`t = (m_{\min} - \mu)/\sigma` and
    :math:`D = 1 - \Phi(t)`,

    .. math::

       \frac{\partial \ln N_{\rm lt}}{\partial\mu} = \frac{z}{\sigma}
           - \frac{\phi(t)}{\sigma D},
       \qquad
       \frac{\partial \ln N_{\rm lt}}{\partial\sigma} = \frac{z^2 - 1}{\sigma}
           - \frac{t\,\phi(t)}{\sigma D}.

    The :math:`\phi(t)/D` terms are the truncation edge responding to the peak
    moving or widening. They are not small: at :math:`\mu = 10`,
    :math:`\sigma = 1.5`, :math:`m_{\min} = 5` the truncation already removes
    0.04% of the Gaussian, and the term grows fast as a peak approaches the
    cutoff.

    Returns
    -------
    dict
        Keys ``"mu"``, ``"sigma"``, each an array shaped like ``mass``.
    """
    mass = numpy.asarray(mass, dtype=float)
    scaled = (mass - mu) / sigma
    edge = (mmin - mu) / sigma
    survival = ndtr(-edge)
    edge_gaussian = _INVERSE_ROOT_TWO_PI * numpy.exp(-0.5 * edge**2)
    ratio = edge_gaussian / (sigma * survival)
    return {
        "mu": scaled / sigma - ratio,
        "sigma": (scaled**2 - 1.0) / sigma - edge * ratio,
    }


def smoothing_ddelta_m(mass, mmin, mmax, delta_m):
    r"""``dS/d(delta_m)`` for the Planck taper.

    With :math:`x = m - m_{\min}`, :math:`d = \delta_m` and
    :math:`f = d/x + d/(x - d)` so that :math:`S = (1 + e^{f})^{-1}`,

    .. math::

       \frac{\partial f}{\partial d} = \frac{1}{x} + \frac{x}{(x - d)^2},
       \qquad
       \frac{\partial S}{\partial d} = -S(1 - S)\,\frac{\partial f}{\partial d}.

    Zero outside the ramp, where :math:`S(1-S)` vanishes anyway.

    Returns
    -------
    numpy.ndarray
    """
    mass = numpy.asarray(mass, dtype=float)
    derivative = numpy.zeros(numpy.shape(mass))
    if delta_m == 0:
        return derivative
    offset = mass - mmin
    ramp = (offset > 0.0) & (offset < delta_m) & (mass <= mmax)
    if not numpy.any(ramp):
        return derivative
    inside = offset[ramp]
    window = smoothing(mass[ramp], mmin=mmin, mmax=mmax, delta_m=delta_m)
    dexponent = 1.0 / inside + inside / (inside - delta_m) ** 2
    derivative[ramp] = -window * (1.0 - window) * dexponent
    return derivative


def smoothing_dmass(mass, mmin, mmax, delta_m):
    r"""``dS/dm`` for the Planck taper, at fixed hyper-parameters.

    Same window as :func:`smoothing_ddelta_m`, differentiated the other way:
    :math:`\partial f/\partial m = -\big(d/x^2 + d/(x-d)^2\big)`, so

    .. math::

       \frac{\partial S}{\partial m}
         = S(1-S)\left(\frac{d}{x^2} + \frac{d}{(x-d)^2}\right).

    Returns
    -------
    numpy.ndarray
    """
    mass = numpy.asarray(mass, dtype=float)
    derivative = numpy.zeros(numpy.shape(mass))
    if delta_m == 0:
        return derivative
    offset = mass - mmin
    ramp = (offset > 0.0) & (offset < delta_m) & (mass <= mmax)
    if not numpy.any(ramp):
        return derivative
    inside = offset[ramp]
    window = smoothing(mass[ramp], mmin=mmin, mmax=mmax, delta_m=delta_m)
    dexponent = delta_m / inside**2 + delta_m / (inside - delta_m) ** 2
    derivative[ramp] = window * (1.0 - window) * dexponent
    return derivative


class BrokenPowerLawTwoPeakMass(PopulationModel):
    r"""BGP primary mass and power-law mass ratio, with analytic scores.

    Usage
    -----
    >>> model = BrokenPowerLawTwoPeakMass(mmin=5.0, m_high=100.0)
    >>> model.log_prob({"mass_1_source": m1, "mass_ratio": q}, model.fiducial)
    >>> model.score({"mass_1_source": m1, "mass_ratio": q})

    Attributes
    ----------
    parameter_names : list of str
        The eleven free hyper-parameters, using
        :class:`GWForge.population.mass.Mass`'s spelling so that a generation
        config and a Fisher config read the same.
    mmin, m_high, maximum_mass : float
        Support parameters, fixed by construction. See the module docstring.
    """

    parameter_names = [
        "alpha_1",
        "alpha_2",
        "m_break",
        "lam_0",
        "lam_1",
        "mpp_1",
        "sigpp_1",
        "mpp_2",
        "sigpp_2",
        "delta_m",
        "beta",
    ]
    event_keys = ["mass_1_source"]

    def __init__(
        self,
        mmin=BGP_PARAMETERS["mmin"],
        m_high=BGP_PARAMETERS["m_high"],
        maximum_mass=BGP_PARAMETERS["maximum_mass"],
        mmin_2=BGP_PARAMETERS["mmin_2"],
        delta_m_2=BGP_PARAMETERS["delta_m_2"],
        fiducial=None,
        nodes=MASS_GRID_NODES,
    ):
        """
        Parameters
        ----------
        mmin : float
            Lower edge of the population, and the foot of the Planck taper.
        m_high : float
            Upper edge of the broken power law. The Gaussian peaks are only
            left-truncated and extend past it.
        maximum_mass : float
            Upper edge of the support and of the normalisation grid. Matches
            ``Mass``'s optional ``maximum_mass``.
        mmin_2, delta_m_2 : float
            The secondary's own Planck taper, independent of the primary's.
        fiducial : dict or None
            Defaults to the GWTC-5.0 medians, :data:`GWForge.population.mass.BGP_PARAMETERS`.
        nodes : int
            Nodes on the shared mass grid.
        """
        self.mmin = float(mmin)
        # GWTC-5.0 Tab. 5 gives the secondary its own taper, (m2_low, delta_m_2),
        # independent of the primary's. Both are support parameters here, like
        # mmin itself: they move the edge of the distribution, so their scores
        # carry Leibniz boundary terms and their centred scores vanish anyway.
        # Both default to the GWTC-5.0 medians. Pass mmin_2=mmin, delta_m_2=delta_m
        # for the single-taper model the other smoothed mass distributions use.
        self.mmin_2 = float(mmin_2)
        self.m_high = float(m_high)
        self.maximum_mass = float(maximum_mass)
        self.fiducial = dict(
            {name: BGP_PARAMETERS[name] for name in self.parameter_names}
            if fiducial is None
            else fiducial
        )
        self.delta_m_2 = float(delta_m_2)
        # The grid must reach down to the *secondary's* edge: p(q | m_1) turns
        # on there, below the primary's.
        self.grid = mass_grid(min(self.mmin, self.mmin_2), self.maximum_mass, nodes)

    # -- primary mass ----------------------------------------------------

    def _components(self, mass, parameters):
        """The three normalised mixture components and the taper, at ``mass``."""
        power_law = broken_power_law(
            mass,
            alpha_1=parameters["alpha_1"],
            alpha_2=parameters["alpha_2"],
            m_break=parameters["m_break"],
            mmin=self.mmin,
            m_high=self.m_high,
        )
        peak_1 = left_truncated_normal(
            mass, parameters["mpp_1"], parameters["sigpp_1"], low=self.mmin
        )
        peak_2 = left_truncated_normal(
            mass, parameters["mpp_2"], parameters["sigpp_2"], low=self.mmin
        )
        window = smoothing(
            mass, mmin=self.mmin, mmax=self.maximum_mass, delta_m=parameters["delta_m"]
        )
        return power_law, peak_1, peak_2, window

    def _primary_numerator(self, mass, parameters):
        """Un-normalised :math:`\\pi(m_1)`: the mixture times the taper."""
        power_law, peak_1, peak_2, window = self._components(mass, parameters)
        weight_2 = 1.0 - parameters["lam_0"] - parameters["lam_1"]
        mixture = (
            parameters["lam_0"] * power_law
            + parameters["lam_1"] * peak_1
            + weight_2 * peak_2
        )
        return mixture * window

    def _primary_numerator_derivatives(self, mass, parameters):
        """``d(numerator)/dLambda`` at ``mass``, as a dict of arrays."""
        power_law, peak_1, peak_2, window = self._components(mass, parameters)
        weight_0 = parameters["lam_0"]
        weight_1 = parameters["lam_1"]
        weight_2 = 1.0 - weight_0 - weight_1
        mixture = weight_0 * power_law + weight_1 * peak_1 + weight_2 * peak_2

        norm = broken_power_law_norm(
            parameters["alpha_1"],
            parameters["alpha_2"],
            parameters["m_break"],
            self.mmin,
            self.m_high,
        )
        dnorm = broken_power_law_norm_derivatives(
            parameters["alpha_1"],
            parameters["alpha_2"],
            parameters["m_break"],
            self.mmin,
            self.m_high,
        )
        # d(shape)/dLambda for the unnormalised broken power law, recovered from
        # the normalised one: shape = power_law * norm.
        mass = numpy.asarray(mass, dtype=float)
        m_break = parameters["m_break"]
        low = (mass >= self.mmin) & (mass < m_break)
        high = (mass >= m_break) & (mass < self.m_high)
        log_ratio = numpy.zeros(numpy.shape(mass))
        support = low | high
        log_ratio[support] = numpy.log(mass[support] / m_break)
        shape = power_law * norm
        dshape = {
            "alpha_1": numpy.where(low, -log_ratio * shape, 0.0),
            "alpha_2": numpy.where(high, -log_ratio * shape, 0.0),
            "m_break": numpy.where(low, parameters["alpha_1"] * shape / m_break, 0.0)
            + numpy.where(high, parameters["alpha_2"] * shape / m_break, 0.0),
        }

        derivatives = {}
        for name in ("alpha_1", "alpha_2", "m_break"):
            # Quotient rule on shape / norm.
            derivatives[name] = (
                weight_0
                * (dshape[name] / norm - shape * dnorm[name] / norm**2)
                * window
            )
        derivatives["lam_0"] = (power_law - peak_2) * window
        derivatives["lam_1"] = (peak_1 - peak_2) * window

        first = left_truncated_normal_dlog(
            mass, parameters["mpp_1"], parameters["sigpp_1"], self.mmin
        )
        second = left_truncated_normal_dlog(
            mass, parameters["mpp_2"], parameters["sigpp_2"], self.mmin
        )
        derivatives["mpp_1"] = weight_1 * peak_1 * first["mu"] * window
        derivatives["sigpp_1"] = weight_1 * peak_1 * first["sigma"] * window
        derivatives["mpp_2"] = weight_2 * peak_2 * second["mu"] * window
        derivatives["sigpp_2"] = weight_2 * peak_2 * second["sigma"] * window
        derivatives["delta_m"] = mixture * smoothing_ddelta_m(
            mass, self.mmin, self.maximum_mass, parameters["delta_m"]
        )
        derivatives["beta"] = numpy.zeros(numpy.shape(mass))
        return derivatives

    # -- mass ratio ------------------------------------------------------

    def _cumulative_secondary(self, parameters):
        r"""The cumulative secondary-mass integral and its derivatives.

        .. math::

           G(M) = \int_{m_{\min}}^{M} m_2^{\beta}\, S(m_2)\, dm_2,

        so that :math:`C(m_1) = m_1^{-(\beta+1)} G(m_1)` is the mass-ratio
        normalisation. Substituting :math:`m_2 = m_1 q` turns one integral per
        event into one cumulative integral on the shared grid, interpolated --
        which is also how ``_smoothed_mass._norm_p_q`` does it, so the two agree.

        Returns
        -------
        dict
            ``{"value": ..., "beta": ..., "delta_m": ...}``, each a callable
            cubic-Hermite interpolant of the corresponding cumulative integral.
        """
        grid = self.grid
        window = smoothing(
            grid, mmin=self.mmin_2, mmax=self.maximum_mass, delta_m=self.delta_m_2
        )
        power = grid ** parameters["beta"]
        # log(grid) is safe: the grid starts at mmin > 0.
        # Only beta is a free parameter of G now: the secondary taper is fixed,
        # so delta_m acts on the primary alone.
        integrands = {
            "value": power * window,
            "beta": power * numpy.log(grid) * window,
        }
        # Cubic Hermite rather than linear interpolation, because the slope of
        # each of these integrals is known exactly -- it is the integrand. Just
        # above mmin, G is a very small number rising very fast, and a linear
        # interpolant there is wrong enough that d ln p / d m_1 disagreed with a
        # finite difference by 4.7e-2 (against a column scale of order 3). With
        # the exact slopes imposed at every node that falls to 6.7e-5 at the very
        # foot of the taper, and to 6e-11 by m_1 = 20 Msun.
        return {
            name: CubicHermiteSpline(
                grid,
                cumulative_trapezoid(integrand, grid, initial=0.0),
                integrand,
                extrapolate=False,
            )
            for name, integrand in integrands.items()
        }

    def dlog_prob_dmasses(self, events, parameters):
        r"""``(d ln p / d m_1, d ln p / d m_2)`` at fixed hyper-parameters.

        Not needed by the Fisher itself -- these are derivatives with respect to
        the *event* parameters, not the hyper-parameters. They are needed by
        :mod:`GWForge.population_fisher.model`, where a trial cosmology
        moves the redshift assigned to a measured distance and so moves the
        source-frame masses:
        :math:`\partial_\Lambda m^{\rm src} = -m^{\rm src}/(1+z)\cdot\partial_\Lambda z`.

        .. math::

           \frac{\partial \ln p}{\partial m_1}
             = \frac{\partial \ln(\text{mixture})}{\partial m_1}
             + \frac{\partial \ln S(m_1)}{\partial m_1}
             - \frac{m_1^{\beta} S(m_1)}{G(m_1)},

        the last term because :math:`m_1` is the upper limit of the mass-ratio
        normalisation :math:`G`, and

        .. math::

           \frac{\partial \ln p}{\partial m_2}
             = \frac{\beta}{m_2} + \frac{\partial \ln S(m_2)}{\partial m_2}.

        Returns
        -------
        tuple of numpy.ndarray
            ``(d ln p / d m_1, d ln p / d m_2)``, ``nan`` outside the support.
        """
        self.check_events(events)
        primary = numpy.asarray(events["mass_1_source"], dtype=float)
        secondary = self._secondary_mass(events)
        inside = self._support(primary, secondary)
        safe_primary = numpy.where(inside, primary, self.mmin + 1.0)
        safe_secondary = numpy.where(inside, secondary, self.mmin_2 + 0.5)

        power_law, peak_1, peak_2, window = self._components(safe_primary, parameters)
        weight_0 = parameters["lam_0"]
        weight_1 = parameters["lam_1"]
        weight_2 = 1.0 - weight_0 - weight_1
        mixture = weight_0 * power_law + weight_1 * peak_1 + weight_2 * peak_2

        m_break = parameters["m_break"]
        low = (safe_primary >= self.mmin) & (safe_primary < m_break)
        high = (safe_primary >= m_break) & (safe_primary < self.m_high)
        dpower_law = -power_law / safe_primary * (
            numpy.where(low, parameters["alpha_1"], 0.0)
            + numpy.where(high, parameters["alpha_2"], 0.0)
        )
        dpeak_1 = -peak_1 * (safe_primary - parameters["mpp_1"]) / parameters["sigpp_1"] ** 2
        dpeak_2 = -peak_2 * (safe_primary - parameters["mpp_2"]) / parameters["sigpp_2"] ** 2
        dmixture = weight_0 * dpower_law + weight_1 * dpeak_1 + weight_2 * dpeak_2

        cumulative = self._cumulative_secondary(parameters)
        secondary_norm = cumulative["value"](safe_primary)
        secondary_window = smoothing(
            safe_secondary,
            mmin=self.mmin_2,
            mmax=self.maximum_mass,
            delta_m=self.delta_m_2,
        )

        with numpy.errstate(divide="ignore", invalid="ignore"):
            primary_derivative = (
                dmixture / mixture
                + smoothing_dmass(
                    safe_primary, self.mmin, self.maximum_mass, parameters["delta_m"]
                )
                / window
                # m_1 is the upper limit of the mass-ratio normalisation G, so
                # it enters through the *secondary's* taper evaluated there.
                - safe_primary ** parameters["beta"]
                * smoothing(
                    safe_primary,
                    mmin=self.mmin_2,
                    mmax=self.maximum_mass,
                    delta_m=self.delta_m_2,
                )
                / secondary_norm
            )
            secondary_derivative = (
                parameters["beta"] / safe_secondary
                + smoothing_dmass(
                    safe_secondary, self.mmin_2, self.maximum_mass, self.delta_m_2
                )
                / secondary_window
            )
        return (
            numpy.where(inside, primary_derivative, numpy.nan),
            numpy.where(inside, secondary_derivative, numpy.nan),
        )

    @staticmethod
    def _secondary_mass(events):
        """``m_2`` from either ``mass_2_source`` or ``mass_ratio``."""
        primary = numpy.asarray(events["mass_1_source"], dtype=float)
        if "mass_2_source" in events:
            return numpy.asarray(events["mass_2_source"], dtype=float)
        return numpy.asarray(events["mass_ratio"], dtype=float) * primary

    # -- interface -------------------------------------------------------

    def _support(self, primary, secondary):
        return (
            (primary >= self.mmin)
            & (primary <= self.maximum_mass)
            & (secondary >= self.mmin_2)
            & (secondary <= primary)
        )

    def log_prob(self, events, parameters):
        """See :meth:`GWForge.population_fisher.base.PopulationModel.log_prob`."""
        self.check_events(events)
        primary = numpy.asarray(events["mass_1_source"], dtype=float)
        secondary = self._secondary_mass(events)
        inside = self._support(primary, secondary)
        safe_primary = numpy.where(inside, primary, self.mmin + 1.0)
        safe_secondary = numpy.where(inside, secondary, self.mmin_2 + 0.5)

        numerator = self._primary_numerator(safe_primary, parameters)
        normalisation = numpy.trapezoid(
            self._primary_numerator(self.grid, parameters), self.grid
        )

        cumulative = self._cumulative_secondary(parameters)
        secondary_norm = cumulative["value"](safe_primary)
        secondary_window = smoothing(
            safe_secondary,
            mmin=self.mmin_2,
            mmax=self.maximum_mass,
            delta_m=self.delta_m_2,
        )

        with numpy.errstate(divide="ignore", invalid="ignore"):
            value = (
                numpy.log(numerator)
                - numpy.log(normalisation)
                + parameters["beta"] * numpy.log(safe_secondary)
                + numpy.log(secondary_window)
                - numpy.log(secondary_norm)
            )
        return numpy.where(inside & numpy.isfinite(value), value, -numpy.inf)

    def analytic_score(self, events, parameters):
        r"""Closed-form score, ``(n_events, 11)``.

        Each column is the log-derivative of the numerator minus that of the
        normalisation,
        :math:`\partial_\Lambda \ln(\text{numerator}) - \partial_\Lambda \ln Z`,
        for the primary mass, plus the mass-ratio
        contribution for the two parameters (``beta`` and ``delta_m``) that
        appear in both. :math:`\partial_\Lambda Z` is the trapezoid of
        :math:`\partial_\Lambda(\text{numerator})` on the *same* nodes that
        normalise ``log_prob``, which is exact because the nodes do not move.
        """
        self.check_events(events)
        primary = numpy.asarray(events["mass_1_source"], dtype=float)
        secondary = self._secondary_mass(events)
        inside = self._support(primary, secondary)
        safe_primary = numpy.where(inside, primary, self.mmin + 1.0)
        safe_secondary = numpy.where(inside, secondary, self.mmin_2 + 0.5)

        numerator = self._primary_numerator(safe_primary, parameters)
        event_derivatives = self._primary_numerator_derivatives(
            safe_primary, parameters
        )
        grid_numerator = self._primary_numerator(self.grid, parameters)
        grid_derivatives = self._primary_numerator_derivatives(self.grid, parameters)
        normalisation = numpy.trapezoid(grid_numerator, self.grid)

        cumulative = self._cumulative_secondary(parameters)
        secondary_norm = cumulative["value"](safe_primary)

        columns = []
        for name in self.parameter_names:
            with numpy.errstate(divide="ignore", invalid="ignore"):
                column = event_derivatives[name] / numerator - numpy.trapezoid(
                    grid_derivatives[name], self.grid
                ) / normalisation
            if name == "beta":
                column = (
                    column
                    + numpy.log(safe_secondary)
                    - cumulative["beta"](safe_primary) / secondary_norm
                )
            # delta_m acts on the primary alone: the secondary has its own,
            # fixed taper, so p(q | m_1) does not depend on it and there is
            # nothing to add to the column computed above.
            columns.append(column)

        score = numpy.column_stack(columns)
        score[~inside, :] = numpy.nan
        score[~numpy.all(numpy.isfinite(score), axis=1), :] = numpy.nan
        return score


# =========================================================================
# Spin: the Default BBH model
# =========================================================================

# Upper edge of the spin-magnitude support, pinned by a delta-function prior.
DEFAULT_MAXIMUM_SPIN = 1.0

# Lower edge of the cos-tilt support. The paper leaves ``t_min`` free
# (Tab. 6); this variant of the analysis pins it, and so does GWForge.
DEFAULT_MINIMUM_COSINE_TILT = -1.0

class DefaultSpin(PopulationModel):
    r"""GWTC-5.0 Default BBH spins: truncated-Gaussian magnitudes, mixture tilts.

    Usage
    -----
    >>> model = DefaultSpin()
    >>> model.log_prob({"a_1": a1, "a_2": a2, "tilt_1": t1, "tilt_2": t2}, model.fiducial)

    Attributes
    ----------
    parameter_names : list of str
        ``["mu_chi", "sigma_chi", "mu_t", "sigma_t", "xi_spin"]`` -- the names
        ``GWForge.population.spin.Spin`` takes in ``spin-parameters``. The
        posterior files call ``mu_t``/``sigma_t`` ``mu_spin``/``sigma_spin``.
    """

    parameter_names = ["mu_chi", "sigma_chi", "mu_t", "sigma_t", "xi_spin"]
    event_keys = ["a_1", "a_2"]

    def __init__(
        self,
        fiducial=None,
        maximum_spin=DEFAULT_MAXIMUM_SPIN,
        minimum_cosine_tilt=DEFAULT_MINIMUM_COSINE_TILT,
    ):
        """
        Parameters
        ----------
        fiducial : dict or None
            Defaults to the GWTC-5.0 medians,
            :data:`GWForge.population.spin.DEFAULT_BBH_SPIN_PARAMETERS`.
        maximum_spin : float
            ``amax``; pinned at 1 in the analysis.
        minimum_cosine_tilt : float
            ``t_min``; pinned at -1 in this variant.
        """
        self.fiducial = dict(
            {name: DEFAULT_BBH_SPIN_PARAMETERS[name] for name in self.parameter_names}
            if fiducial is None
            else fiducial
        )
        self.maximum_spin = float(maximum_spin)
        self.minimum_cosine_tilt = float(minimum_cosine_tilt)
        # Density of the isotropic component *of the pair*, not of one tilt.
        self.isotropic_density = 1.0 / (1.0 - self.minimum_cosine_tilt) ** 2

    # -- helpers ---------------------------------------------------------

    @staticmethod
    def _cosine_tilts(events):
        """``(cos tilt_1, cos tilt_2)``, from either name."""
        cosines = []
        for index in (1, 2):
            cosine_key = "cos_tilt_{}".format(index)
            if cosine_key in events:
                cosines.append(numpy.asarray(events[cosine_key], dtype=float))
            else:
                cosines.append(
                    numpy.cos(
                        numpy.asarray(events["tilt_{}".format(index)], dtype=float)
                    )
                )
        return cosines

    def _magnitudes(self, events):
        return [
            numpy.asarray(events["a_{}".format(index)], dtype=float)
            for index in (1, 2)
        ]

    def check_events(self, events):
        """Raise if the magnitudes or either spelling of the tilts is missing."""
        missing = [key for key in ("a_1", "a_2") if key not in events]
        for index in (1, 2):
            if not (
                "cos_tilt_{}".format(index) in events
                or "tilt_{}".format(index) in events
            ):
                missing.append("tilt_{} (or cos_tilt_{})".format(index, index))
        if missing:
            raise KeyError(
                "DefaultSpin needs event key(s) {}; got {}.".format(
                    missing, sorted(events)
                )
            )

    def _magnitude_terms(self, events, parameters):
        """Truncated-normal terms for both spin magnitudes."""
        return [
            truncated_normal(
                magnitude,
                parameters["mu_chi"],
                parameters["sigma_chi"],
                0.0,
                self.maximum_spin,
            )
            for magnitude in self._magnitudes(events)
        ]

    def _tilt_terms(self, events, parameters):
        """Truncated-normal terms for both cosine tilts."""
        return [
            truncated_normal(
                cosine,
                parameters["mu_t"],
                parameters["sigma_t"],
                self.minimum_cosine_tilt,
                1.0,
            )
            for cosine in self._cosine_tilts(events)
        ]

    # -- interface -------------------------------------------------------

    def log_prob(self, events, parameters):
        """See :meth:`GWForge.population_fisher.base.PopulationModel.log_prob`."""
        self.check_events(events)
        magnitudes = self._magnitude_terms(events, parameters)
        tilts = self._tilt_terms(events, parameters)

        # The tilt term is the density of the *pair*, not a product of two
        # per-component densities. See the module docstring.
        joint_tilt = (
            parameters["xi_spin"] * tilts[0]["density"] * tilts[1]["density"]
            + (1.0 - parameters["xi_spin"]) * self.isotropic_density
        )
        inside = (
            (magnitudes[0]["density"] > 0.0)
            & (magnitudes[1]["density"] > 0.0)
            & (joint_tilt > 0.0)
        )
        with numpy.errstate(divide="ignore", invalid="ignore"):
            total = (
                numpy.log(numpy.where(inside, magnitudes[0]["density"], 1.0))
                + numpy.log(numpy.where(inside, magnitudes[1]["density"], 1.0))
                + numpy.log(numpy.where(inside, joint_tilt, 1.0))
            )
        return numpy.where(inside, total, -numpy.inf)

    def analytic_score(self, events, parameters):
        r"""Closed-form score, ``(n_events, 5)``.

        The magnitude columns are the sum of the two components'
        :math:`\partial \ln N/\partial(\mu_\chi, \sigma_\chi)`. The tilt columns
        differentiate the *joint* mixture, so each picks up both terms of the
        product rule:

        .. math::

           \frac{\partial \ln \pi}{\partial \mu_t}
             = \frac{\xi\,N_1 N_2\left(\partial_{\mu_t}\ln N_1
                                     + \partial_{\mu_t}\ln N_2\right)}{\pi},
           \qquad
           \frac{\partial \ln \pi}{\partial \xi}
             = \frac{N_1 N_2 - (1-t_{\min})^{-2}}{\pi}.
        """
        self.check_events(events)
        magnitudes = self._magnitude_terms(events, parameters)
        tilts = self._tilt_terms(events, parameters)
        xi_spin = parameters["xi_spin"]

        product = tilts[0]["density"] * tilts[1]["density"]
        joint_tilt = xi_spin * product + (1.0 - xi_spin) * self.isotropic_density

        inside = (
            (magnitudes[0]["density"] > 0.0)
            & (magnitudes[1]["density"] > 0.0)
            & (joint_tilt > 0.0)
        )
        safe = numpy.where(inside, joint_tilt, 1.0)

        columns = numpy.column_stack(
            [
                magnitudes[0]["mu"] + magnitudes[1]["mu"],
                magnitudes[0]["sigma"] + magnitudes[1]["sigma"],
                xi_spin * product * (tilts[0]["mu"] + tilts[1]["mu"]) / safe,
                xi_spin * product * (tilts[0]["sigma"] + tilts[1]["sigma"]) / safe,
                (product - self.isotropic_density) / safe,
            ]
        )
        columns[~inside, :] = numpy.nan
        return columns


# =========================================================================
# Redshift: Madau-Dickinson
# =========================================================================

# Trapezoidal nodes for the normalisation Z and the expectations in the score.
# dV_c/dz peaks broadly around z ~ 2 and psi is smooth, so this is far more
# than convergence needs; measured drift in sigma between 2000 and 8000 nodes
# is below 1e-6 relative.
NORMALISATION_NODES = 4000


class MadauDickinsonRedshift(PopulationModel):
    r"""Merger-rate density following the Madau-Dickinson shape.

    Usage
    -----
    >>> model = MadauDickinsonRedshift(maximum_redshift=10.0)
    >>> model.log_prob({"redshift": z}, model.fiducial)
    >>> model.score({"redshift": z})

    Attributes
    ----------
    parameter_names : list of str
        ``["gamma", "kappa", "z_peak"]`` -- the names
        :class:`GWForge.population.redshift.Redshift` uses, so a generation
        config and a Fisher config read the same.
    """

    parameter_names = ["gamma", "kappa", "z_peak"]
    event_keys = ["redshift"]

    def __init__(
        self,
        maximum_redshift=10.0,
        cosmology=None,
        fiducial=None,
        nodes=NORMALISATION_NODES,
    ):
        """
        Parameters
        ----------
        maximum_redshift : float
            Upper edge of the support. Must match the ``maximum-redshift`` the
            population was generated with, or the normalisation is of a
            different distribution.
        cosmology : FlatwCDM or astropy.cosmology.FLRW or None
            Supplies ``dV_c/dz``. Defaults to Planck18.
        fiducial : dict or None
            Defaults to the GWForge generation defaults
            ``{"gamma": 2.7, "kappa": 5.6, "z_peak": 1.9}``.
        nodes : int
            Trapezoidal nodes for the normalisation.
        """
        if cosmology is None:
            from ..cosmology import astropy_cosmology

            cosmology = astropy_cosmology("Planck18")
        self.maximum_redshift = float(maximum_redshift)
        self.cosmology = cosmology
        self.fiducial = dict(
            {"gamma": 2.7, "kappa": 5.6, "z_peak": 1.9}
            if fiducial is None
            else fiducial
        )
        # The grid starts just above zero: dV_c/dz vanishes at z = 0, and its
        # logarithm is what the normalisation integrates.
        self.nodes = numpy.linspace(1e-6, self.maximum_redshift, int(nodes))
        self._log_measure_nodes = numpy.log(
            differential_comoving_volume(self.cosmology, self.nodes)
        ) - numpy.log1p(self.nodes)

    def log_shape(self, redshift, parameters):
        r"""``ln psi(z)``: the Madau-Dickinson shape alone, without the measure.

        Public because :mod:`GWForge.population_fisher.model` rebuilds
        the measure at a trial cosmology and needs the shape on its own.

        Includes the :math:`[1 + (1+z_p)^{-\kappa}]` factor that sets
        :math:`\psi(0) = 1`. That factor is independent of :math:`z`, so it
        multiplies the numerator and the normalisation alike and cancels from
        ``log_prob`` identically -- which is why :meth:`shape_gradient` omits it.
        """
        return numpy.log(
            madau_dickinson_psi_of_z(
                redshift,
                gamma=parameters["gamma"],
                kappa=parameters["kappa"],
                z_peak=parameters["z_peak"],
            )
        )

    def _log_normalisation(self, parameters):
        """``ln Z``, the trapezoidal integral of the numerator over the support."""
        integrand = numpy.exp(
            self.log_shape(self.nodes, parameters) + self._log_measure_nodes
        )
        return numpy.log(numpy.trapezoid(integrand, self.nodes))

    def log_prob(self, events, parameters):
        """See :meth:`GWForge.population_fisher.base.PopulationModel.log_prob`."""
        self.check_events(events)
        redshift = numpy.asarray(events["redshift"], dtype=float)
        inside = (redshift >= 0.0) & (redshift <= self.maximum_redshift)
        log_measure = numpy.where(
            inside,
            numpy.log(
                differential_comoving_volume(
                    self.cosmology, numpy.where(inside, redshift, 1.0)
                )
            )
            - numpy.log1p(numpy.where(inside, redshift, 1.0)),
            0.0,
        )
        value = (
            self.log_shape(numpy.where(inside, redshift, 1.0), parameters)
            + log_measure
            - self._log_normalisation(parameters)
        )
        return numpy.where(inside, value, -numpy.inf)

    def shape_gradient(self, redshift, parameters):
        r"""``d ln psi / dLambda``, shape ``(n, 3)``.

        Public for the same reason as :meth:`log_shape`.

        With :math:`u = (1+z)/(1+z_p)` and :math:`r = u^\kappa`,
        :math:`f = r/(1+r)`:

        .. math::

           \partial_\gamma = \ln(1+z), \quad
           \partial_\kappa = -f \ln u, \quad
           \partial_{z_p} = \frac{\kappa f}{1 + z_p}.
        """
        redshift = numpy.asarray(redshift, dtype=float)
        z_peak = parameters["z_peak"]
        kappa = parameters["kappa"]
        ratio = (1.0 + redshift) / (1.0 + z_peak)
        power = ratio**kappa
        fraction = power / (1.0 + power)
        return numpy.column_stack(
            [
                numpy.log1p(redshift),
                -fraction * numpy.log(ratio),
                kappa * fraction / (1.0 + z_peak),
            ]
        )

    def analytic_score(self, events, parameters):
        r"""Closed-form score, ``(n_events, 3)``.

        ``g(z) - E_p[g]``: the first term from the numerator, the second from
        :math:`\partial_\Lambda \ln Z`, evaluated on the normalisation grid.
        """
        self.check_events(events)
        redshift = numpy.asarray(events["redshift"], dtype=float)
        inside = (redshift >= 0.0) & (redshift <= self.maximum_redshift)
        safe = numpy.where(inside, redshift, 1.0)

        density = numpy.exp(
            self.log_shape(self.nodes, parameters) + self._log_measure_nodes
        )
        weights = density / numpy.trapezoid(density, self.nodes)
        expectation = numpy.trapezoid(
            weights[:, numpy.newaxis] * self.shape_gradient(self.nodes, parameters),
            self.nodes,
            axis=0,
        )

        score = self.shape_gradient(safe, parameters) - expectation
        score[~inside, :] = numpy.nan
        return score


# =========================================================================
# Spectral sirens
# =========================================================================

# Detector-frame observables this model is a density over.
OBSERVABLE_KEYS = ("mass_1", "mass_2", "luminosity_distance")


class SpectralSirenModel(PopulationModel):
    r"""Source-frame population observed in the detector frame, cosmology free.

    Usage
    -----
    >>> model = SpectralSirenModel(mass_model, redshift_model, spin_model)
    >>> model.score({"mass_1": m1_det, "mass_2": m2_det,
    ...              "luminosity_distance": distance})

    Attributes
    ----------
    parameter_names : list of str
        Mass, then redshift, then spin, then the cosmology parameters.
    maximum_redshift : float
        Upper edge of the redshift support, taken from the redshift model.
    """

    def __init__(
        self,
        mass_model,
        redshift_model,
        spin_model=None,
        cosmology=None,
        cosmology_parameters=COSMOLOGY_PARAMETERS,
        nodes=None,
    ):
        """
        Parameters
        ----------
        mass_model : GWForge.population_fisher.model.BrokenPowerLawTwoPeakMass
            Source-frame mass model. Must provide ``dlog_prob_dmasses``.
        redshift_model : GWForge.population_fisher.model.MadauDickinsonRedshift
            Supplies the rate *shape* and its gradient. Its own ``dV_c/dz`` is
            not used -- the measure is rebuilt at the trial cosmology.
        spin_model : GWForge.population_fisher.base.PopulationModel or None
            Optional, frame-invariant.
        cosmology : GWForge.cosmology.FlatwCDM or None
            Fiducial cosmology. Defaults to flat LCDM at Planck18's
            ``(H_0, \\Omega_{m,0} + \\Omega_{\\nu,0})``.
        cosmology_parameters : sequence of str
            Subset of :data:`GWForge.cosmology.COSMOLOGY_PARAMETERS` to expose.
        nodes : int or None
            Redshift nodes for the measure normalisation. Defaults to the
            redshift model's own grid size.
        """
        if cosmology is None:
            cosmology = FlatwCDM(H0=67.66, Om0=0.3111, w0=-1.0)
        self.mass_model = mass_model
        self.redshift_model = redshift_model
        self.spin_model = spin_model
        self.cosmology_parameters = list(cosmology_parameters)
        self.maximum_redshift = redshift_model.maximum_redshift

        blocks = [mass_model, redshift_model]
        if spin_model is not None:
            blocks.append(spin_model)
        self.blocks = blocks

        self.parameter_names = []
        for block in blocks:
            self.parameter_names.extend(block.parameter_names)
        self.parameter_names.extend(self.cosmology_parameters)
        self._cosmology_fiducial = {
            name: cosmology.parameters()[name] for name in self.cosmology_parameters
        }

        self.event_keys = list(OBSERVABLE_KEYS)
        if spin_model is not None:
            self.event_keys.extend(
                key for key in spin_model.event_keys if key not in self.event_keys
            )
        count = redshift_model.nodes.size if nodes is None else int(nodes)
        self.nodes = numpy.linspace(1e-6, self.maximum_redshift, count)

    @property
    def fiducial(self):
        """Block fiducials plus the cosmology, merged live.

        A property for the same reason as
        :attr:`GWForge.population_fisher.base.JointPopulationModel.fiducial`:
        the redshift block is always refitted after the model is built, and a
        snapshot would send the Fisher back to the pre-fit point.
        """
        merged = {}
        for block in self.blocks:
            merged.update(block.fiducial)
        merged.update(self._cosmology_fiducial)
        return merged

    # -- helpers ---------------------------------------------------------

    def _cosmology(self, parameters):
        """The trial cosmology, cosmology parameters taken from ``parameters``."""
        settings = {"H0": 67.66, "Om0": 0.3111, "w0": -1.0}
        settings.update(
            {name: float(parameters[name]) for name in self.cosmology_parameters}
        )
        return FlatwCDM(**settings)

    def _source_frame(self, events, cosmology):
        """``(z, m_1^src, m_2^src)`` at a trial cosmology."""
        distance = numpy.asarray(events["luminosity_distance"], dtype=float)
        redshift = cosmology.redshift_of_distance(distance)
        one_plus_z = 1.0 + redshift
        return (
            redshift,
            numpy.asarray(events["mass_1"], dtype=float) / one_plus_z,
            numpy.asarray(events["mass_2"], dtype=float) / one_plus_z,
        )

    def _log_redshift_numerator(self, redshift, parameters, cosmology):
        r"""``ln psi(z) + ln(dV_c/dz) - ln(1+z)``, the measure at the trial cosmology."""
        return (
            self.redshift_model.log_shape(redshift, parameters)
            + numpy.log(cosmology.differential_comoving_volume(redshift))
            - numpy.log1p(redshift)
        )

    def _redshift_weights(self, parameters, cosmology):
        """Normalised ``p(z)`` on :attr:`nodes`, and its normalisation."""
        density = numpy.exp(
            self._log_redshift_numerator(self.nodes, parameters, cosmology)
        )
        normalisation = numpy.trapezoid(density, self.nodes)
        return density / normalisation, normalisation

    # -- interface -------------------------------------------------------

    def log_prob(self, events, parameters):
        """See :meth:`GWForge.population_fisher.base.PopulationModel.log_prob`."""
        self.check_events(events)
        cosmology = self._cosmology(parameters)
        redshift, primary, secondary = self._source_frame(events, cosmology)
        inside = (redshift > 0.0) & (redshift <= self.maximum_redshift)
        safe = numpy.where(inside, redshift, 1.0)

        _, normalisation = self._redshift_weights(parameters, cosmology)
        total = (
            self.mass_model.log_prob(
                {"mass_1_source": primary, "mass_2_source": secondary}, parameters
            )
            + self._log_redshift_numerator(safe, parameters, cosmology)
            - numpy.log(normalisation)
            - 2.0 * numpy.log1p(safe)
            - numpy.log(cosmology.ddL_dz(safe))
        )
        if self.spin_model is not None:
            total = total + self.spin_model.log_prob(events, parameters)
        return numpy.where(inside, total, -numpy.inf)

    def analytic_score(self, events, parameters):
        """Closed-form score, ``(n_events, len(parameter_names))``.

        The population columns are the source-frame scores evaluated at the
        source-frame parameters the trial cosmology implies; the cosmology
        columns are the chain rule in the module docstring.
        """
        self.check_events(events)
        cosmology = self._cosmology(parameters)
        redshift, primary, secondary = self._source_frame(events, cosmology)
        inside = (redshift > 0.0) & (redshift <= self.maximum_redshift)
        safe = numpy.where(inside, redshift, 1.0)
        source_events = {"mass_1_source": primary, "mass_2_source": secondary}

        columns = [self.mass_model.analytic_score(source_events, parameters)]

        # Redshift shape: g(z) - E_p[g], the expectation taken under the trial
        # cosmology's measure rather than the redshift model's fixed one.
        weights, _ = self._redshift_weights(parameters, cosmology)
        gradient_nodes = self.redshift_model.shape_gradient(self.nodes, parameters)
        expectation = numpy.trapezoid(
            weights[:, numpy.newaxis] * gradient_nodes, self.nodes, axis=0
        )
        columns.append(
            self.redshift_model.shape_gradient(safe, parameters) - expectation
        )

        if self.spin_model is not None:
            columns.append(self.spin_model.analytic_score(events, parameters))

        columns.append(self._cosmology_score(safe, primary, secondary, parameters, cosmology))

        score = numpy.hstack(columns)
        score[~inside, :] = numpy.nan
        return score

    def _cosmology_score(self, redshift, primary, secondary, parameters, cosmology):
        """The ``(H0, Om0, w0)`` columns; see the chain rule in the module docstring."""
        one_plus_z = 1.0 + redshift
        efunc = cosmology.efunc(redshift)
        comoving = cosmology.comoving_distance(redshift)
        ddL_dz = cosmology.ddL_dz(redshift)

        # d/dz of the terms that depend on z at fixed cosmology.
        dlog_volume_dz = (
            2.0 * cosmology.hubble_distance / (efunc * comoving)
            - cosmology.defunc_dz(redshift) / efunc
        )
        dlog_ddL_dz = cosmology.d2dL_dz2(redshift) / ddL_dz
        dlog_mass_dm1, dlog_mass_dm2 = self.mass_model.dlog_prob_dmasses(
            {"mass_1_source": primary, "mass_2_source": secondary}, parameters
        )
        shape_gradient_z = self._dlog_shape_dz(redshift, parameters)

        sensitivity = (
            dlog_mass_dm1 * (-primary / one_plus_z)
            + dlog_mass_dm2 * (-secondary / one_plus_z)
            + shape_gradient_z
            + dlog_volume_dz
            - 3.0 / one_plus_z
            - dlog_ddL_dz
        )

        derivatives = cosmology.derivatives(
            redshift, parameters=self.cosmology_parameters
        )
        dredshift = cosmology.dredshift_dparameter(
            redshift, parameters=self.cosmology_parameters
        )
        # The normalisation Z_z depends on the cosmology through dV_c/dz only,
        # and its derivative is the expectation of the explicit term over p(z).
        weights, _ = self._redshift_weights(parameters, cosmology)
        node_volume = cosmology.differential_comoving_volume(self.nodes)
        node_derivatives = cosmology.derivatives(
            self.nodes, parameters=self.cosmology_parameters
        )

        volume = cosmology.differential_comoving_volume(redshift)
        columns = []
        for name in self.cosmology_parameters:
            explicit = (
                derivatives[name]["differential_comoving_volume"] / volume
                - derivatives[name]["ddL_dz"] / ddL_dz
            )
            dlog_normalisation = numpy.trapezoid(
                weights
                * node_derivatives[name]["differential_comoving_volume"]
                / node_volume,
                self.nodes,
            )
            columns.append(
                sensitivity * dredshift[name] + explicit - dlog_normalisation
            )
        return numpy.column_stack(columns)

    @staticmethod
    def _dlog_shape_dz(redshift, parameters):
        r""":math:`\partial \ln\psi/\partial z = (\gamma - \kappa f)/(1+z)`.

        With :math:`r = ((1+z)/(1+z_p))^{\kappa}` and :math:`f = r/(1+r)`.
        """
        ratio = (1.0 + redshift) / (1.0 + parameters["z_peak"])
        power = ratio ** parameters["kappa"]
        fraction = power / (1.0 + power)
        return (parameters["gamma"] - parameters["kappa"] * fraction) / (
            1.0 + redshift
        )
