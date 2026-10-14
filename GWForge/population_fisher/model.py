r"""The population models GWForge forecasts on, with their analytic scores."""

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
    r""":math:`A(x, \alpha) = \int_x^1 t^{-\alpha} dt = (1 - x^{1-\alpha})/(1-\alpha)`."""
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
    r"""Derivative of :func:`_power_law_integral` with respect to the slope."""
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

    Attributes
    ----------
    parameter_names : list of str
        The eleven free hyper-parameters, using
        :class:`GWForge.population.mass.Mass`'s spelling so that a generation
        config and a Fisher config read the same.
    mmin, m_high, maximum_mass : float
        Support parameters, fixed by construction: hard cutoffs have a
        constant score, so their centred score vanishes and they are not
        measurable within Term I.
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
        r"""Closed-form score, ``(n_events, 11)``."""
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
        # per-component densities.
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
        r"""Closed-form score, ``(n_events, 5)``."""
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
        r"""``ln psi(z)``: the Madau-Dickinson shape alone, without the measure."""
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
        r"""``d ln psi / dLambda``, shape ``(n, 3)``."""
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
        r"""Closed-form score, ``(n_events, 3)``."""
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
        """Block fiducials plus the cosmology, merged live."""
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
        """Closed-form score, ``(n_events, len(parameter_names))``."""
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
        """The ``(H0, Om0, w0)`` columns, by the chain rule through ``z(d_L)``."""
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
        r""":math:`\partial \ln\psi/\partial z = (\gamma - \kappa f)/(1+z)`."""
        ratio = (1.0 + redshift) / (1.0 + parameters["z_peak"])
        power = ratio ** parameters["kappa"]
        fraction = power / (1.0 + power)
        return (parameters["gamma"] - parameters["kappa"] * fraction) / (
            1.0 + redshift
        )
