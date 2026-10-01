r"""Maximum-likelihood fits of a population model to a set of samples.

This exists for one necessary reason and one useful one.

**Necessary.** :class:`GWForge.population.redshift.Redshift` always convolves
the Madau-Dickinson star-formation rate with a formation-to-merger time-delay
distribution, so a population it generates is *not* distributed as
:math:`\psi_{\rm MD}(z)\,dV_c/dz/(1+z)`. Running a Fisher forecast for the
direct Madau-Dickinson form on such a catalogue at the generation
:math:`(\gamma, \kappa, z_p)` would be a forecast about the wrong fiducial
point. Fitting the direct form to the injected redshifts gives the *effective*
parameters that population actually corresponds to, and those are what the
Fisher should be evaluated at.

**Useful.** The mass and spin blocks need no such fit: ``gwforge_population``
draws them straight from the BGP and ``Default`` densities, so the fiducials are
the config values. Running the fit on them anyway is therefore a test with
teeth -- it must return the config values to within Monte-Carlo error, and if it
does not, the log density here disagrees with the sampler over there.

Both use the analytic gradient, so an eleven-parameter BGP fit to 30,000 samples
is seconds rather than minutes.

Simplex parameters
------------------

``lam_0`` and ``lam_1`` are mixture weights and must satisfy
:math:`\lambda_0, \lambda_1 \ge 0` and :math:`\lambda_0 + \lambda_1 \le 1`; a
box constraint cannot express the last one. They are therefore fitted through a
softmax of two unconstrained logits, which stays in the simplex by construction
and has no boundary for the optimiser to stick to. The chain rule for the
gradient is in :class:`_SimplexReparameterisation`.
"""

import logging

import numpy
from scipy.optimize import minimize

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Mixture weights fitted through a softmax rather than a box.
SIMPLEX_PARAMETERS = ("lam_0", "lam_1")

# L-BFGS-B settings. ``ftol`` is tight because the objective is a mean log
# density -- an O(1) number whose interesting variation is in the fifth digit.
FIT_OPTIONS = {"maxiter": 3000, "ftol": 1e-12, "gtol": 1e-10}

# Objective value returned when a trial point puts samples outside the model's
# support. Large and **finite**, not ``inf``: L-BFGS-B's line search aborts on a
# non-finite value and reports convergence, so an ``inf`` here does not push the
# optimiser back inside the feasible region -- it stops it dead at the starting
# point, which looks exactly like a perfect fit.
#
# This is not hypothetical. The Beta spin model needs
# a constraint no box can express -- the Beta models need
# ``sigma_squared_chi < mu_chi (1 - mu_chi)``, which no box constraint can
# express, so L-BFGS-B's first trial step leaves the region every time. With
# ``inf`` the spin block returned its starting values to twelve significant
# figures from any start; with this it converges normally.
OUT_OF_SUPPORT_PENALTY = 1e8

# Default box constraints, used when the caller gives none. Wide enough not to
# bind at any sensible fiducial, narrow enough to keep the optimiser out of
# regions where the density is undefined.
DEFAULT_BOUNDS = {
    "gamma": (-2.0, 12.0),
    "kappa": (0.1, 20.0),
    "z_peak": (0.1, 8.0),
    "alpha_1": (-3.0, 12.0),
    "alpha_2": (-3.0, 15.0),
    "m_break": (6.0, 95.0),
    "mpp_1": (5.0, 120.0),
    "sigpp_1": (0.3, 30.0),
    "mpp_2": (5.0, 120.0),
    "sigpp_2": (0.3, 30.0),
    "delta_m": (0.01, 20.0),
    "beta": (-4.0, 10.0),
    # GWTC-5.0 Tab. 6 priors for the Default BBH spin model. Note sigma_chi is
    # a standard deviation; the variance-parameterised sigma_squared_chi belongs
    # to the Beta models.
    "mu_chi": (0.0, 1.0),
    "sigma_chi": (0.005, 1.0),
    "mu_t": (-1.0, 1.0),
    "sigma_t": (0.01, 4.0),
    "xi_spin": (0.0, 1.0),
    "sigma_squared_chi": (1e-4, 0.2),
}


class _SimplexReparameterisation:
    r"""Softmax map from two logits to ``(lam_0, lam_1)``.

    With :math:`w = (v_0, v_1, 0)` and :math:`\lambda_i = e^{w_i}/\sum_j e^{w_j}`,

    .. math::

       \frac{\partial \lambda_i}{\partial v_j}
         = \lambda_i(\delta_{ij} - \lambda_j),

    which is the standard softmax Jacobian restricted to the two free logits.
    Fixing the third logit at zero removes the map's one redundant direction.
    """

    names = SIMPLEX_PARAMETERS

    @staticmethod
    def to_weights(logits):
        """Two logits to ``(lam_0, lam_1, lam_2)``."""
        exponentials = numpy.exp(
            numpy.append(numpy.asarray(logits, dtype=float), 0.0)
            - numpy.max(numpy.append(numpy.asarray(logits, dtype=float), 0.0))
        )
        return exponentials / exponentials.sum()

    @staticmethod
    def to_logits(lam_0, lam_1):
        """``(lam_0, lam_1)`` back to two logits."""
        lam_2 = 1.0 - lam_0 - lam_1
        if lam_2 <= 0.0:
            raise ValueError(
                "lam_0 + lam_1 must be below 1 for the third mixture weight to "
                "exist; got {} + {}.".format(lam_0, lam_1)
            )
        return numpy.array([numpy.log(lam_0 / lam_2), numpy.log(lam_1 / lam_2)])

    @staticmethod
    def jacobian(weights):
        """``d(lam_0, lam_1)/d(logit_0, logit_1)``, a ``(2, 2)`` array."""
        first_two = weights[:2]
        return numpy.diag(first_two) - numpy.outer(first_two, first_two)


def fit_model(
    model,
    samples,
    free_parameters=None,
    start=None,
    bounds=None,
    fixed_parameters=None,
):
    """Fit a population model to samples by maximum likelihood.

    Parameters
    ----------
    model : GWForge.population_fisher.base.PopulationModel
    samples : dict
        Event parameters to fit to. For a fiducial, these should be *all*
        injections, not just the detected ones: the population model describes
        the astrophysical population, and fitting to the detected subset would
        absorb the selection function into the hyper-parameters.
    free_parameters : sequence of str or None
        Defaults to every hyper-parameter not in ``fixed_parameters``.
    start : dict or None
        Starting point. Defaults to ``model.fiducial``.
    bounds : dict or None
        Name to ``(low, high)``. Merged onto :data:`DEFAULT_BOUNDS`.
    fixed_parameters : dict or None

    Returns
    -------
    dict
        The fitted values, merged onto the starting point, with the keys
        ``"_success"``, ``"_message"`` and ``"_mean_log_prob"`` describing the
        optimisation.
    """
    fixed_parameters = dict(fixed_parameters or {})
    start = dict(model.fiducial if start is None else start)
    start.update(fixed_parameters)
    if free_parameters is None:
        free_parameters = [
            name for name in model.parameter_names if name not in fixed_parameters
        ]
    free_parameters = list(free_parameters)

    limits = dict(DEFAULT_BOUNDS)
    limits.update(bounds or {})

    simplex = [name for name in SIMPLEX_PARAMETERS if name in free_parameters]
    use_simplex = len(simplex) == len(SIMPLEX_PARAMETERS)
    if simplex and not use_simplex:
        raise ValueError(
            "Fit {} together or not at all: fitting one mixture weight without "
            "the other leaves the simplex constraint unexpressible.".format(
                list(SIMPLEX_PARAMETERS)
            )
        )
    plain = [name for name in free_parameters if name not in SIMPLEX_PARAMETERS]

    def unpack(vector):
        """Optimiser vector to a full hyper-parameter dict."""
        trial = dict(start)
        trial.update(dict(zip(plain, vector[: len(plain)])))
        if use_simplex:
            weights = _SimplexReparameterisation.to_weights(vector[len(plain) :])
            trial["lam_0"] = weights[0]
            trial["lam_1"] = weights[1]
        return trial

    def objective(vector):
        trial = unpack(vector)
        values = model.log_prob(samples, trial)
        if not numpy.all(numpy.isfinite(values)):
            return OUT_OF_SUPPORT_PENALTY
        return -float(numpy.mean(values))

    def gradient(vector):
        trial = unpack(vector)
        scores = model.score(samples, trial, method="auto")
        finite = numpy.all(numpy.isfinite(scores), axis=1)
        if not numpy.any(finite):
            return numpy.zeros(len(vector))
        # Average over the samples that are still in support. Outside the
        # feasible region the objective is the flat penalty above and the line
        # search is driving the step down anyway; this at least keeps the
        # direction meaningful on the boundary.
        mean_score = -numpy.mean(scores[finite], axis=0)
        indices = {name: i for i, name in enumerate(model.parameter_names)}
        columns = [mean_score[indices[name]] for name in plain]
        if use_simplex:
            weights = _SimplexReparameterisation.to_weights(vector[len(plain) :])
            weight_gradient = numpy.array(
                [mean_score[indices[name]] for name in SIMPLEX_PARAMETERS]
            )
            columns.extend(
                _SimplexReparameterisation.jacobian(weights).T @ weight_gradient
            )
        return numpy.array(columns)

    initial = [float(start[name]) for name in plain]
    box = [limits.get(name, (-numpy.inf, numpy.inf)) for name in plain]
    if use_simplex:
        initial.extend(
            _SimplexReparameterisation.to_logits(start["lam_0"], start["lam_1"])
        )
        box.extend([(-20.0, 20.0)] * len(SIMPLEX_PARAMETERS))

    # A single out-of-support sample makes the objective -inf *everywhere*, so
    # L-BFGS-B cannot take a first step and returns the starting point -- which
    # looks exactly like a perfect fit. Drop them at the start, loudly, so a
    # model that disagrees with the catalogue shows up as a warning rather than
    # as suspiciously exact agreement.
    samples, dropped = _restrict_to_support(model, samples, start)
    total = len(numpy.atleast_1d(samples[model.event_keys[0]])) + dropped
    if dropped:
        logging.warning(
            "%d of %d samples (%.2f%%) are outside %s's support at the starting "
            "point and were dropped from the fit. A large fraction means the "
            "model and the catalogue describe different populations.",
            dropped,
            total,
            100.0 * dropped / total,
            type(model).__name__,
        )

    logging.info(
        "Fitting %s to %d samples over %s",
        type(model).__name__,
        len(numpy.atleast_1d(samples[model.event_keys[0]])),
        free_parameters,
    )
    result = minimize(
        objective,
        numpy.array(initial, dtype=float),
        jac=gradient,
        method="L-BFGS-B",
        bounds=box,
        options=FIT_OPTIONS,
    )
    fitted = unpack(result.x)
    if not result.success:
        logging.warning("Fit did not converge cleanly: %s", result.message)
    fitted["_success"] = bool(result.success)
    fitted["_message"] = str(result.message)
    fitted["_mean_log_prob"] = -float(result.fun)
    fitted["_n_dropped"] = dropped
    fitted["_n_samples"] = total - dropped
    logging.info(
        "Fitted: %s",
        ", ".join("{}={:.5g}".format(name, fitted[name]) for name in free_parameters),
    )
    return fitted


def _restrict_to_support(model, samples, parameters):
    """Drop samples the model gives zero probability at ``parameters``.

    Returns
    -------
    tuple
        ``(samples, n_dropped)``. ``samples`` is the original object when
        nothing was dropped.
    """
    values = model.log_prob(samples, parameters)
    finite = numpy.isfinite(values)
    dropped = int((~finite).sum())
    if not dropped:
        return samples, 0
    return {
        key: numpy.asarray(value)[finite]
        for key, value in samples.items()
        if numpy.ndim(value) == 1 and len(numpy.atleast_1d(value)) == len(finite)
    }, dropped


def compare_to(fitted, reference, names):
    """A plain-text table of fitted values against the values used to generate.

    Parameters
    ----------
    fitted : dict
    reference : dict
    names : sequence of str

    Returns
    -------
    str
    """
    lines = [
        "  {:>18s}   {:>12s}   {:>12s}   {:>10s}".format(
            "Parameter", "Generated", "Fitted", "ratio"
        ),
        "  " + "-" * 60,
    ]
    for name in names:
        generated = float(reference[name])
        value = float(fitted[name])
        ratio = value / generated if abs(generated) > 1e-30 else float("nan")
        lines.append(
            "  {:>18s}   {:>12.5g}   {:>12.5g}   {:>10.4f}".format(
                name, generated, value, ratio
            )
        )
    return "\n".join(lines)
