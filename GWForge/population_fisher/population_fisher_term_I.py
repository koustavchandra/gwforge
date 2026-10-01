r"""The Term-I population Fisher matrix.

Gair, Ghosh et al. (`arXiv:2205.07893 <https://arxiv.org/abs/2205.07893>`_)
expand the hyper-parameter Fisher into five terms. This module implements the
first, which dominates whenever per-event measurement uncertainties are small
compared with the width of the population distribution -- the regime XG
detectors put most sources in.

For :math:`N_{\rm det}` detected events :math:`\theta_k`,

.. math::

   s_k^i = \frac{\partial \ln p(\theta_k \mid \Lambda)}{\partial \Lambda^i},
   \qquad
   \langle s \rangle^i = \frac{1}{N_{\rm det}}\sum_k s_k^i,
   \qquad
   \Gamma^{ij} = \sum_k \big(s_k^i - \langle s\rangle^i\big)
                        \big(s_k^j - \langle s\rangle^j\big),

and :math:`\Sigma = \Gamma^{-1}`.

Where the selection function went
---------------------------------

Not away -- it is never *evaluated*, which is a different thing. Term I is an
expectation under the **detected** density, so :math:`p_{\rm det}` is in it
twice over: as the measure, and through the detectable fraction
:math:`\alpha(\Lambda)`. Two facts collapse both into the centred sum above.

First, :math:`p_{\rm det}(\theta)` carries no hyper-parameter, so it drops from
the score. Second, the detected catalogue is itself a fair Monte-Carlo draw from
:math:`p_{\rm det}(\theta) p(\theta|\Lambda)/\alpha(\Lambda)`, so

.. math::

   \frac{\partial \ln \alpha}{\partial \Lambda^i}
     = \int p_{\rm det}(\theta)\,\frac{p(\theta|\Lambda)}{\alpha(\Lambda)}\,
       \frac{\partial \ln p}{\partial \Lambda^i}\,d\theta
     = \mathbb{E}_{\rm det}\!\left[s^i\right],

which the sample mean estimates. So there is no separate selection integral to
evaluate, no injection reweighting and no rate normalisation -- but the
selection is still doing work, carried entirely by *which events are in the
sum*. The SNR threshold is an input to the forecast, not a preprocessing step:
change it and :math:`\Gamma` changes. The absolute scale is set by how many
events were detected, and every :math:`\sigma` therefore scales as
:math:`1/\sqrt{T_{\rm obs}}`.

``validation/population_detection_fraction.py`` reports what a given threshold
removes. :math:`p_{\rm det}` reappears explicitly at higher order -- as
:math:`\partial_\theta \ln p_{\rm det}` in Gair's :math:`\Gamma_4` -- which is
one reason this module stops at Term I.

Two consequences worth knowing. First, any :math:`\Lambda`-dependent but
*event-independent* normalisation cancels identically under centring, so a model
that forgets one still gives the right Fisher (it will not give the right
``log_prob``, which is why the models here compute it anyway). Second, the
identity is exact only if :math:`p_{\rm det}` carries no hyper-parameter -- the
reason :mod:`GWForge.population_fisher.model` builds its density in the
detector frame rather than pushing events to the source frame at a trial
cosmology.

Degenerate directions
---------------------

A parameter whose score is the *same number* for every event -- a hard mass
cutoff, in the interior of a power law -- has an identically zero centred score,
so its Fisher row and column vanish and the matrix is singular. That is a
statement about the model, not a numerical accident: within Term I such a
parameter is unmeasurable. :func:`population_fisher` detects it, names it, and
falls back to a pseudo-inverse rather than returning a plausible-looking number.
"""

import logging

import numpy

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Above this condition number the inverse is not trustworthy in double
# precision and the pseudo-inverse is used instead, with a warning.
CONDITION_NUMBER_WARNING = 1e14

# Singular-value threshold for the pseudo-inverse fallback.
DEFAULT_RCOND = 1e-12

# A centred-score column whose root-mean-square is below this fraction of the
# largest column's is reported as degenerate.
DEGENERACY_THRESHOLD = 1e-10


class PopulationFisherResult:
    """Everything :func:`population_fisher` computed, and how to read it.

    Attributes
    ----------
    parameter_names : list of str
        Free hyper-parameters, in matrix order.
    fixed_parameters : dict
        Hyper-parameters held fixed, excluded from the matrix.
    fiducial : dict
        Values the Fisher was evaluated at, free and fixed together.
    fisher : numpy.ndarray
        ``(n, n)`` Term-I Fisher matrix.
    covariance : numpy.ndarray
        Its inverse (or pseudo-inverse).
    sigma : dict
        Name to one-sigma marginalised uncertainty.
    condition_number : float
    n_events : int
        Events with finite scores, i.e. those that entered the sum.
    n_total : int
        Events offered, before the support and finiteness cuts.
    scores : numpy.ndarray
        ``(n_events, n)`` raw, uncentred scores.
    mean_score : numpy.ndarray
        ``(n,)``; the Monte-Carlo estimate of ``d ln alpha / d Lambda``.
    degenerate : list of str
        Free parameters whose centred score vanished.
    """

    def __init__(
        self,
        parameter_names,
        fixed_parameters,
        fiducial,
        fisher,
        covariance,
        condition_number,
        n_events,
        n_total,
        scores,
        mean_score,
        degenerate,
    ):
        self.parameter_names = list(parameter_names)
        self.fixed_parameters = dict(fixed_parameters)
        self.fiducial = dict(fiducial)
        self.fisher = fisher
        self.covariance = covariance
        self.condition_number = condition_number
        self.n_events = int(n_events)
        self.n_total = int(n_total)
        self.scores = scores
        self.mean_score = mean_score
        self.degenerate = list(degenerate)
        self.sigma = dict(
            zip(self.parameter_names, numpy.sqrt(numpy.abs(numpy.diag(covariance))))
        )

    @classmethod
    def from_npz(cls, path):
        """Reload a forecast that :file:`bin/gwforge_population_fisher` saved.

        Comparing two networks, or redrawing a figure in a new style, should not
        mean recomputing either forecast. The saved archive carries everything a
        plot needs: the matrix, its covariance, the parameter names, the
        fiducial values and the event counts.

        The score matrix is *not* saved -- it is one row per detected event and
        would dominate the file -- so the reloaded result has an empty
        ``scores`` and ``mean_score``. Everything the plotting path touches is
        present.

        Parameters
        ----------
        path : str
            A ``.npz`` written by the CLI.

        Returns
        -------
        PopulationFisherResult
        """
        with numpy.load(path, allow_pickle=True) as archive:
            names = [str(name) for name in archive["parameter_names"]]
            fiducial = dict(zip(names, archive["fiducial"].astype(float)))
            return cls(
                parameter_names=names,
                fixed_parameters={},
                fiducial=fiducial,
                fisher=archive["fisher"],
                covariance=archive["covariance"],
                condition_number=(
                    float(archive["condition_number"])
                    if "condition_number" in archive
                    else float("nan")
                ),
                n_events=int(archive["n_events"]),
                n_total=int(archive["n_total"]),
                scores=numpy.empty((0, len(names))),
                mean_score=numpy.zeros(len(names)),
                degenerate=[],
            )

    def correlation_matrix(self):
        """Correlation matrix implied by :attr:`covariance`."""
        deviations = numpy.sqrt(numpy.abs(numpy.diag(self.covariance)))
        with numpy.errstate(divide="ignore", invalid="ignore"):
            return self.covariance / numpy.outer(deviations, deviations)

    def scaled_to(self, n_events):
        """This forecast rescaled to a different number of detected events.

        The Term-I Fisher is a sum over events, so it is proportional to
        :math:`N_{\\rm det}` in expectation and :math:`\\sigma \\propto
        N_{\\rm det}^{-1/2}`. Use this to turn a one-month catalogue into a
        ten-year forecast without regenerating anything -- or, more usefully, to
        state plainly that a headline number is a pure :math:`\\sqrt{T}`
        extrapolation rather than a computed one.

        Parameters
        ----------
        n_events : int

        Returns
        -------
        PopulationFisherResult
        """
        factor = float(n_events) / self.n_events
        return PopulationFisherResult(
            parameter_names=self.parameter_names,
            fixed_parameters=self.fixed_parameters,
            fiducial=self.fiducial,
            fisher=self.fisher * factor,
            covariance=self.covariance / factor,
            condition_number=self.condition_number,
            n_events=n_events,
            n_total=int(round(self.n_total * factor)),
            scores=self.scores,
            mean_score=self.mean_score,
            degenerate=self.degenerate,
        )

    def summary(self):
        """A plain-text table of the forecast.

        Returns
        -------
        str
        """
        rule = "-" * 74
        lines = [
            rule,
            "Population Fisher (Gair+2022 Term I)",
            rule,
            "  Events used      : {} / {} offered".format(self.n_events, self.n_total),
            "  Condition number : {:.4g}".format(self.condition_number),
        ]
        if self.fixed_parameters:
            lines.append(
                "  Fixed            : "
                + ", ".join(
                    "{}={:g}".format(name, value)
                    for name, value in sorted(self.fixed_parameters.items())
                )
            )
        if self.degenerate:
            lines.append("  Degenerate       : " + ", ".join(self.degenerate))
        lines += [
            "",
            "  {:>18s}   {:>12s}   {:>12s}   {:>11s}".format(
                "Parameter", "Fiducial", "sigma", "sigma/|fid|"
            ),
            "  " + "-" * 62,
        ]
        for index, name in enumerate(self.parameter_names):
            fiducial = float(self.fiducial[name])
            deviation = self.sigma[name]
            relative = (
                deviation / abs(fiducial) if abs(fiducial) > 1e-30 else float("nan")
            )
            lines.append(
                "  {:>18s}   {:>12.5g}   {:>12.5g}   {:>11.4g}".format(
                    name, fiducial, deviation, relative
                )
            )
        lines.append(rule)
        return "\n".join(lines)

    def __repr__(self):
        return "PopulationFisherResult({} parameters, {} events)".format(
            len(self.parameter_names), self.n_events
        )


def covariance_from_fisher(fisher, names=None, rcond=DEFAULT_RCOND):
    """Invert a Fisher matrix, falling back to a pseudo-inverse when it is singular.

    The matrix is normalised by its diagonal before inversion -- the same trick
    :func:`GWForge.fisher.matrix.covariance` uses -- because hyper-parameters
    differ by many orders of magnitude in scale (``H0`` near 68, ``lam_1`` near
    0.05) and the raw condition number is then dominated by units rather than by
    degeneracy.

    Parameters
    ----------
    fisher : numpy.ndarray
    names : sequence of str or None
        Used only to name parameters in the warning.
    rcond : float

    Returns
    -------
    tuple
        ``(covariance, condition_number)``.
    """
    diagonal = numpy.diag(fisher).astype(float)
    scale = numpy.sqrt(numpy.where(diagonal > 0.0, diagonal, 1.0))
    normalised = fisher / numpy.outer(scale, scale)

    eigenvalues = numpy.linalg.eigvalsh(normalised)
    if eigenvalues[0] > 0.0:
        condition_number = float(eigenvalues[-1] / eigenvalues[0])
    else:
        condition_number = numpy.inf

    if numpy.isfinite(condition_number) and condition_number < 1.0 / rcond:
        inverse = numpy.linalg.inv(normalised)
    else:
        logging.warning(
            "Population Fisher is ill-conditioned (condition number %.4g); using "
            "the pseudo-inverse. Parameters: %s",
            condition_number,
            list(names) if names is not None else "unnamed",
        )
        inverse = numpy.linalg.pinv(normalised, rcond=rcond)
    return inverse / numpy.outer(scale, scale), condition_number


def population_fisher(
    model,
    events,
    free_parameters=None,
    fixed_parameters=None,
    parameters=None,
    method="auto",
    rcond=DEFAULT_RCOND,
    n_total=None,
):
    """Term-I population Fisher matrix for a detected catalogue.

    Parameters
    ----------
    model : GWForge.population_fisher.base.PopulationModel
    events : dict
        Source parameters of the **detected** events. They must be drawn from
        the population at ``parameters``, or the centred-score identity that
        supplies the selection term does not hold.
    free_parameters : sequence of str or None
        Hyper-parameters to vary. Defaults to every name the model has that is
        not in ``fixed_parameters``.
    fixed_parameters : dict or None
        Name to value, held fixed and excluded from the matrix.
    parameters : dict or None
        Fiducial point. Defaults to ``model.fiducial``.
    method : str
        Passed to :meth:`~GWForge.population_fisher.base.PopulationModel.score`.
    rcond : float
        Pseudo-inverse threshold.
    n_total : int or None
        Events offered before cuts, for reporting. Defaults to ``len(events)``.

    Returns
    -------
    PopulationFisherResult

    Raises
    ------
    ValueError
        If no parameters are free, or no event has a finite score.
    """
    fixed_parameters = dict(fixed_parameters or {})
    fiducial = dict(model.fiducial if parameters is None else parameters)
    fiducial.update(fixed_parameters)

    if free_parameters is None:
        free_parameters = [
            name for name in model.parameter_names if name not in fixed_parameters
        ]
    free_parameters = list(free_parameters)
    unknown = [name for name in free_parameters if name not in model.parameter_names]
    if unknown:
        raise ValueError(
            "{} is not a hyper-parameter of {}. Known: {}".format(
                unknown, type(model).__name__, model.parameter_names
            )
        )
    overlap = [name for name in free_parameters if name in fixed_parameters]
    if overlap:
        raise ValueError("{} cannot be both free and fixed.".format(overlap))
    if not free_parameters:
        raise ValueError("Every hyper-parameter is fixed; there is nothing to compute.")

    # Anything the model knows about but that was not asked for is held at its
    # fiducial value -- which, after a fit, is the *fitted* value, not whatever
    # a config file said. Record it all as fixed so the summary states what was
    # pinned rather than leaving it to be inferred from what is missing.
    fixed_parameters = {
        name: fiducial[name]
        for name in model.parameter_names
        if name not in free_parameters
    }

    logging.info(
        "Population Fisher: model=%s, free=%s, fixed=%s",
        type(model).__name__,
        free_parameters,
        sorted(fixed_parameters),
    )

    all_scores = model.score(events, fiducial, method=method)
    columns = [model.parameter_names.index(name) for name in free_parameters]
    scores = numpy.asarray(all_scores, dtype=float)[:, columns]

    offered = len(scores) if n_total is None else int(n_total)
    finite = numpy.all(numpy.isfinite(scores), axis=1)
    dropped = len(scores) - int(finite.sum())
    if dropped:
        logging.warning(
            "%d of %d events dropped: outside the model's support, or a "
            "non-finite score.",
            dropped,
            len(scores),
        )
    scores = scores[finite]
    if len(scores) == 0:
        raise ValueError(
            "No event has a finite score. Check that the model's support "
            "matches the catalogue -- most often mmin is above the lightest "
            "source, or maximum_redshift is below the most distant one."
        )

    mean_score = scores.mean(axis=0)
    logging.info(
        "d ln alpha / d Lambda (selection term): %s",
        ", ".join(
            "{}={:.4g}".format(name, value)
            for name, value in zip(free_parameters, mean_score)
        ),
    )

    centred = scores - mean_score
    fisher = centred.T @ centred
    fisher = 0.5 * (fisher + fisher.T)

    magnitudes = numpy.sqrt(numpy.mean(centred**2, axis=0))
    largest = magnitudes.max()
    degenerate = [
        name
        for name, magnitude in zip(free_parameters, magnitudes)
        if largest > 0.0 and magnitude <= DEGENERACY_THRESHOLD * largest
    ]
    if degenerate:
        logging.warning(
            "Centred score vanishes for %s: Term I cannot measure these. They "
            "are event-independent directions -- usually a hard cutoff in the "
            "interior of a power law. Pin them via fixed_parameters.",
            degenerate,
        )

    covariance, condition_number = covariance_from_fisher(
        fisher, names=free_parameters, rcond=rcond
    )
    if condition_number > CONDITION_NUMBER_WARNING:
        logging.warning(
            "Condition number %.4g: some directions are close to degenerate and "
            "their sigmas should not be quoted without checking the correlation "
            "matrix.",
            condition_number,
        )

    return PopulationFisherResult(
        parameter_names=free_parameters,
        fixed_parameters=fixed_parameters,
        fiducial=fiducial,
        fisher=fisher,
        covariance=covariance,
        condition_number=condition_number,
        n_events=len(scores),
        n_total=offered,
        scores=scores,
        mean_score=mean_score,
        degenerate=degenerate,
    )
