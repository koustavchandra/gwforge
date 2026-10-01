r"""Figures for population Fisher forecasts.

Two kinds.

The **corner** is the covariance itself: a Fisher forecast is a Gaussian, so the
way to draw it is to realise that Gaussian and hand the samples to
:func:`GWForge.plotting.corner_plot`, which is where the style, the palettes and
the contour machinery live.

The **reconstruction** plots are the more useful ones. A table of sigmas says how
well each hyper-parameter is pinned; a band on :math:`\pi(m_1)` or :math:`p(z)`
says what that means for the distribution anyone actually cares about, and it is
where a degeneracy shows up as a band that is narrow everywhere the parameters
are correlated and wide where they are not. Both overlay the injected histogram,
which is also the quickest way to catch a fiducial that was fitted badly.

Everything here draws inside ``matplotlib.rc_context`` so that importing this
module, or calling anything in it, leaves the caller's matplotlib configuration
alone.
"""

import logging

import matplotlib
import numpy
import pylab

from ..fisher.plot import samples_from_covariance
from ..plotting import (
    LVK_COLOUR,
    XG_COLOUR,
    corner_plot as _corner_plot,
    labels_for,
    new_rcParams,
    palette,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Credible level for the reconstruction bands.
CREDIBLE_LEVEL = 0.90

# Hyper-parameter draws behind each band.
BAND_DRAWS = 2000

# Samples behind each corner plot. Enough that the 90% contour is smooth
# without the histogram in each panel being noisy.
CORNER_SAMPLES = 60000

# Enclosed-probability contour levels, matching the reference figures.
CONTOUR_LEVELS = (0.67, 0.90)

__all__ = [
    "CONTOUR_LEVELS",
    "CREDIBLE_LEVEL",
    "LVK_COLOUR",
    "XG_COLOUR",
    "corner_plot",
    "labels_for",
    "plot_mass_distribution",
    "plot_redshift_distribution",
    "samples_from_covariance",
]


def corner_plot(
    results,
    save=None,
    title=None,
    seed=None,
    parameters=None,
    colours=None,
    palette_name="tealrose",
):
    """Corner plot of one or more population Fisher forecasts, overlaid.

    Parameters
    ----------
    results : sequence of tuple
        ``(label, PopulationFisherResult)`` pairs. Every result must expose the
        same free parameters, so that overlaying them compares like with like.
    save : str or None
        Path to write to. The figure is closed after saving.
    title : str or None
    seed : int or None
        Each dataset is drawn with a **freshly seeded** generator, so they share
        the same underlying standard normals and the contours nest visually
        instead of jittering independently against each other.
    parameters : sequence of str or None
        Restrict the plot to these. Drawn from the *marginal* of the full
        covariance, so the contours are marginalised, not the conditional ones a
        re-fit over a subset would give. Defaults to all of them, which is
        unreadable beyond about eight.
    colours : sequence of str or None
        One per dataset. Defaults to the palette's cycle, which starts
        teal/rose.
    palette_name : str
        ``"tealrose"`` or ``"okabe-ito"``.

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If the results do not share their parameter list, or a requested
        parameter is not among them.
    """
    names = results[0][1].parameter_names
    for label, result in results:
        if result.parameter_names != names:
            raise ValueError(
                "'{}' has parameters {} but the first result has {}; overlaying "
                "them would put different quantities on the same axis.".format(
                    label, result.parameter_names, names
                )
            )
    if parameters is None:
        parameters = names
    missing = [name for name in parameters if name not in names]
    if missing:
        raise ValueError(
            "Cannot plot {}: not among the free parameters {}.".format(missing, names)
        )
    # Marginalising is slicing the *covariance*; slicing the Fisher instead would
    # condition on the parameters left out and give contours that are too tight.
    columns = [names.index(name) for name in parameters]
    centre = numpy.array([results[0][1].fiducial[name] for name in parameters])

    draws, widest = [], numpy.zeros(len(parameters))
    for _, result in results:
        sample = samples_from_covariance(
            [result.fiducial[name] for name in parameters],
            result.covariance[numpy.ix_(columns, columns)],
            size=CORNER_SAMPLES,
            seed=seed,
        )
        draws.append(sample)
        widest = numpy.maximum(widest, sample.std(axis=0))
    # One shared range, set by the *widest* dataset, so the panels line up and
    # the tightest contour is not clipped to a dot.
    span = [
        (centre[i] - 4.0 * widest[i], centre[i] + 4.0 * widest[i])
        for i in range(len(parameters))
    ]

    if colours is None:
        colours = palette(palette_name)[0]

    figure = None
    for index, ((label, _), sample) in enumerate(zip(results, draws)):
        figure = _corner_plot(
            sample,
            labels=labels_for(parameters),
            levels=CONTOUR_LEVELS,
            span=span,
            truths=centre,
            colour=colours[index % len(colours)],
            label=label,
            show_titles=(index == 0),
            palette_name=palette_name,
            fig=figure,
        )
    if title:
        figure.suptitle(title, y=1.02)
    if save:
        figure.savefig(save, bbox_inches="tight")
        logging.info("Wrote %s", save)
        pylab.close(figure)
    return figure


def _band(curves):
    """Median and symmetric credible band of a stack of curves."""
    lower = 0.5 * (1.0 - CREDIBLE_LEVEL)
    return (
        numpy.median(curves, axis=0),
        numpy.quantile(curves, lower, axis=0),
        numpy.quantile(curves, 1.0 - lower, axis=0),
    )


def _draw_parameters(result, size, seed):
    """Hyper-parameter draws from the Fisher covariance, as dicts."""
    mean = [result.fiducial[name] for name in result.parameter_names]
    samples = samples_from_covariance(mean, result.covariance, size=size, seed=seed)
    for row in samples:
        trial = dict(result.fiducial)
        trial.update(dict(zip(result.parameter_names, row)))
        yield trial


def _reconstruction(model, result, evaluate, size, seed):
    """Stack of curves from Fisher draws, dropping any that leave the support."""
    curves = []
    for trial in _draw_parameters(result, size, seed):
        try:
            values = evaluate(model, trial)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            continue
        if numpy.all(numpy.isfinite(values)):
            curves.append(values)
    if not curves:
        raise ValueError(
            "Every hyper-parameter draw left the model's support. The forecast "
            "is too weak to reconstruct a distribution from -- check the "
            "condition number."
        )
    dropped = size - len(curves)
    if dropped:
        logging.warning(
            "%d of %d draws left the model's support and were dropped from the "
            "band. A large fraction here means the Gaussian the Fisher matrix "
            "describes reaches into unphysical parameters.",
            dropped,
            size,
        )
    return numpy.array(curves)


def plot_mass_distribution(
    model,
    result,
    injections=None,
    save=None,
    size=BAND_DRAWS,
    seed=None,
    colour=XG_COLOUR,
):
    """Reconstructed primary-mass distribution with a credible band.

    Parameters
    ----------
    model : GWForge.population_fisher.model.BrokenPowerLawTwoPeakMass
    result : GWForge.population_fisher.population_fisher_term_I.PopulationFisherResult
    injections : array_like or None
        ``mass_1_source`` of the injections, drawn as a histogram behind the
        band.
    save : str or None
    size : int
    seed : int or None
    colour : str
        Band and median colour. Defaults to the XG teal; pass
        :data:`GWForge.plotting.LVK_COLOUR` for a 2G network.

    Returns
    -------
    matplotlib.figure.Figure
    """
    grid = numpy.linspace(model.mmin, min(model.m_high, 120.0), 500)

    def evaluate(mass_model, trial):
        numerator = mass_model._primary_numerator(grid, trial)
        normalisation = numpy.trapezoid(
            mass_model._primary_numerator(mass_model.grid, trial), mass_model.grid
        )
        return numerator / normalisation

    curves = _reconstruction(model, result, evaluate, size, seed)
    median, lower, upper = _band(curves)

    with matplotlib.rc_context(new_rcParams(width="page")):
        figure, axis = pylab.subplots()
        if injections is not None:
            injections = numpy.asarray(injections)
            axis.hist(
                injections[injections <= grid[-1]],
                bins=120,
                density=True,
                histtype="stepfilled",
                color="0.8",
                edgecolor="0.5",
                label="injections",
            )
        axis.plot(
            grid,
            evaluate(model, result.fiducial),
            color="k",
            linewidth=1.2,
            label="fiducial",
        )
        axis.plot(grid, median, color=colour, linewidth=2.2, label="median")
        axis.fill_between(
            grid,
            lower,
            upper,
            color=colour,
            alpha=0.2,
            linewidth=0,
            label="{:.0f}% credible region".format(100 * CREDIBLE_LEVEL),
        )
        # Log, unlike the reference's linear fit-check panel: the BGP mass
        # function spans four decades between the peaks and the m_high cutoff,
        # and a linear axis shows only the low-mass peak.
        axis.set_yscale("log")
        axis.set_xlabel(r"$m_1^{\mathrm{src}}\ [M_{\odot}]$")
        axis.set_ylabel(r"$\pi(m_1)$")
        axis.set_xlim(grid.min(), grid.max())
        axis.set_ylim(max(median.max() * 1e-5, 1e-8), median.max() * 3)
        axis.legend(frameon=False)

    if save:
        figure.savefig(save, bbox_inches="tight")
        logging.info("Wrote %s", save)
        pylab.close(figure)
    return figure


def plot_redshift_distribution(
    model,
    result,
    injections=None,
    save=None,
    size=BAND_DRAWS,
    seed=None,
    colour=XG_COLOUR,
):
    """Reconstructed merger-rate redshift distribution with a credible band.

    Parameters
    ----------
    model : GWForge.population_fisher.model.MadauDickinsonRedshift
    result : GWForge.population_fisher.population_fisher_term_I.PopulationFisherResult
    injections : array_like or None
        Injected redshifts, drawn as a histogram behind the band.
    save : str or None
    size : int
    seed : int or None
    colour : str
        Band and median colour. Defaults to the XG teal; pass
        :data:`GWForge.plotting.LVK_COLOUR` for a 2G network.

    Returns
    -------
    matplotlib.figure.Figure
    """
    grid = numpy.linspace(1e-3, model.maximum_redshift, 400)

    def evaluate(redshift_model, trial):
        return numpy.exp(redshift_model.log_prob({"redshift": grid}, trial))

    curves = _reconstruction(model, result, evaluate, size, seed)
    median, lower, upper = _band(curves)

    with matplotlib.rc_context(new_rcParams(width="page")):
        figure, axis = pylab.subplots()
        if injections is not None:
            axis.hist(
                injections,
                bins=120,
                density=True,
                histtype="stepfilled",
                color="0.8",
                edgecolor="0.5",
                label="injections",
            )
        axis.plot(
            grid,
            evaluate(model, result.fiducial),
            color="k",
            linewidth=1.2,
            label="fiducial",
        )
        axis.plot(grid, median, color=colour, linewidth=2.2, label="median")
        axis.fill_between(
            grid,
            lower,
            upper,
            color=colour,
            alpha=0.2,
            linewidth=0,
            label="{:.0f}% credible region".format(100 * CREDIBLE_LEVEL),
        )
        axis.set_xlabel("$z$")
        axis.set_ylabel("$p(z)$")
        axis.set_xlim(grid.min(), grid.max())
        axis.set_ylim(bottom=0.0)
        axis.legend(frameon=False)

    if save:
        figure.savefig(save, bbox_inches="tight")
        logging.info("Wrote %s", save)
        pylab.close(figure)
    return figure
