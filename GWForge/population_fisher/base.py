r"""Base class for population models, and how their scores get computed.

A population model here is a normalised probability density over the *source*
parameters of a single event, controlled by hyper-parameters :math:`\Lambda`.
The population Fisher only ever needs two things from it: the log density, and
the **score**

.. math::

   s^i(\theta) = \frac{\partial \ln p(\theta \mid \Lambda)}{\partial \Lambda^i}.

Derivative policy
-----------------

Two tiers, in strict order of preference, resolved by :meth:`PopulationModel.score`:

``analytic``
    Closed form, written out by hand. This is the production path for every
    model GWForge ships, cosmology included -- see
    :mod:`GWForge.population_fisher.model` and :mod:`GWForge.cosmology`.

``finite-difference``
    Centred differences. Present as the **test oracle** -- every analytic score
    is checked against it -- and as the fallback for a model added later whose
    derivatives nobody has written out yet. It is not a production path. A
    finite difference on a quantity that is itself a quadrature has an error
    floor set by the quadrature, not by the step, which is exactly the trap the
    previous generation of this code fell into: shrinking the step from 1e-4 to
    1e-5 made the cosmology columns *two orders of magnitude worse*.

There is deliberately no autodiff tier. Every model here has closed-form scores,
so JAX would have been a dependency that the production path never executed.

Support
-------

Events outside a model's support get ``-inf`` from ``log_prob`` and ``nan`` from
``score``; :func:`GWForge.population_fisher.population_fisher_term_I.population_fisher` drops those
rows and says how many. That is deliberate: silently coercing them to zero would
put an event that the model says cannot exist into the Fisher sum with weight
one.
"""

import logging

import numpy

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Score methods :meth:`PopulationModel.score` understands.
SCORE_METHODS = ("auto", "analytic", "finite-difference")

# Relative and absolute finite-difference steps. The step for parameter ``i``
# is ``max(FD_RELATIVE_STEP * |Lambda_i|, FD_ABSOLUTE_STEP)``.
FD_RELATIVE_STEP = 1e-4
FD_ABSOLUTE_STEP = 1e-8


class PopulationModel:
    """Base class: a normalised density over event parameters, given hyper-parameters.

    Subclasses must set :attr:`parameter_names`, :attr:`fiducial` and
    :attr:`event_keys`, and implement :meth:`log_prob`. Implementing
    :meth:`analytic_score` is strongly preferred but optional.

    Attributes
    ----------
    parameter_names : list of str
        Hyper-parameter names, in the order the score columns come out.
    fiducial : dict
        Default hyper-parameter values.
    event_keys : list of str
        Keys the model requires in the ``events`` dict.
    """

    parameter_names = []
    fiducial = {}
    event_keys = []

    # -- interface -------------------------------------------------------

    def log_prob(self, events, parameters):
        """Log density of each event.

        Parameters
        ----------
        events : dict
            Arrays of event parameters, all the same length.
        parameters : dict
            Hyper-parameter values; must cover :attr:`parameter_names`.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_events,)``. ``-inf`` outside the support.
        """
        raise NotImplementedError

    def analytic_score(self, events, parameters):
        """Closed-form score, or ``None`` if this model has none.

        Returns
        -------
        numpy.ndarray or None
            Shape ``(n_events, len(parameter_names))``, columns in
            :attr:`parameter_names` order.
        """
        return None

    # -- score dispatch --------------------------------------------------

    def score(self, events, parameters=None, method="auto"):
        """Score of each event with respect to every hyper-parameter.

        Parameters
        ----------
        events : dict
        parameters : dict or None
            Defaults to :attr:`fiducial`.
        method : str
            One of :data:`SCORE_METHODS`. ``"auto"`` tries analytic, then falls
            back to finite differences.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_events, len(parameter_names))``.
        """
        if method not in SCORE_METHODS:
            raise ValueError(
                "method must be one of {}, got {!r}".format(SCORE_METHODS, method)
            )
        parameters = dict(self.fiducial if parameters is None else parameters)

        if method in ("auto", "analytic"):
            values = self.analytic_score(events, parameters)
            if values is not None:
                return numpy.asarray(values, dtype=float)
            if method == "analytic":
                raise NotImplementedError(
                    "{} has no analytic_score.".format(type(self).__name__)
                )

        return self._finite_difference_score(events, parameters)

    def _finite_difference_score(self, events, parameters):
        """Score by centred finite differences. The test oracle; see the module docstring."""
        names = list(self.parameter_names)
        vector = self.to_vector(parameters)
        n_events = len(numpy.atleast_1d(events[self.event_keys[0]]))
        score = numpy.zeros((n_events, len(names)))

        for index, value in enumerate(vector):
            step = max(FD_RELATIVE_STEP * abs(value), FD_ABSOLUTE_STEP)
            upper = vector.copy()
            upper[index] += step
            lower = vector.copy()
            lower[index] -= step
            above = self.log_prob(events, self.from_vector(upper, parameters))
            below = self.log_prob(events, self.from_vector(lower, parameters))
            finite = numpy.isfinite(above) & numpy.isfinite(below)
            score[finite, index] = (above[finite] - below[finite]) / (2.0 * step)
            score[~finite, index] = numpy.nan
        return score

    # -- helpers ---------------------------------------------------------

    def to_vector(self, parameters):
        """Hyper-parameter dict to array, in :attr:`parameter_names` order."""
        return numpy.array(
            [float(parameters[name]) for name in self.parameter_names], dtype=float
        )

    def from_vector(self, vector, parameters=None):
        """Array back to a hyper-parameter dict, merged onto ``parameters``."""
        merged = dict(self.fiducial if parameters is None else parameters)
        merged.update(dict(zip(self.parameter_names, vector)))
        return merged

    def check_events(self, events):
        """Raise if ``events`` is missing a key this model needs."""
        missing = [key for key in self.event_keys if key not in events]
        if missing:
            raise KeyError(
                "{} needs event key(s) {}; got {}.".format(
                    type(self).__name__, missing, sorted(events)
                )
            )

    def _n_events(self, events):
        return len(numpy.atleast_1d(events[self.event_keys[0]]))

    def _column_index(self, name):
        return self.parameter_names.index(name)


class JointPopulationModel(PopulationModel):
    """Product of independent sub-models over disjoint event parameters.

    The log density adds and the score columns concatenate, in the order the
    sub-models were given. Duplicate hyper-parameter names across sub-models are
    an error rather than a silent tie: two blocks sharing a name almost always
    means one of them was built with the wrong keyword.

    Usage
    -----
    >>> model = JointPopulationModel([mass_model, redshift_model, spin_model])
    >>> model.score(events)
    """

    def __init__(self, models):
        """
        Parameters
        ----------
        models : sequence of PopulationModel
        """
        self.models = list(models)
        names = []
        for model in self.models:
            for name in model.parameter_names:
                if name in names:
                    raise ValueError(
                        "Hyper-parameter '{}' appears in more than one sub-model "
                        "of the joint population model.".format(name)
                    )
                names.append(name)
        self.parameter_names = names
        self.event_keys = []
        for model in self.models:
            for key in model.event_keys:
                if key not in self.event_keys:
                    self.event_keys.append(key)

    @property
    def fiducial(self):
        """Sub-model fiducials, merged live.

        A property rather than a dict copied at construction: the fiducials of
        a block can be *refitted* after the joint model is built (see
        :mod:`GWForge.population_fisher.compare`, which the redshift block
        always needs), and a snapshot taken in ``__init__`` would silently send
        the Fisher back to the pre-fit point.
        """
        merged = {}
        for model in self.models:
            merged.update(model.fiducial)
        return merged

    def log_prob(self, events, parameters):
        """Sum of the sub-model log densities."""
        total = 0.0
        for model in self.models:
            total = total + model.log_prob(events, parameters)
        return total

    def analytic_score(self, events, parameters):
        """Horizontally stacked sub-model scores, ``None`` if any sub-model lacks one."""
        columns = []
        for model in self.models:
            values = model.analytic_score(events, parameters)
            if values is None:
                return None
            columns.append(numpy.asarray(values, dtype=float))
        return numpy.hstack(columns)

    def score(self, events, parameters=None, method="auto"):
        """Per-sub-model dispatch, so one analytic block is not lost to another's fallback."""
        parameters = dict(self.fiducial if parameters is None else parameters)
        return numpy.hstack(
            [model.score(events, parameters, method=method) for model in self.models]
        )
