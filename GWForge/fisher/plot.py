"""Corner plots for Fisher and DALI results.

A Fisher forecast is a covariance matrix, not a set of samples, so the way to
draw it is to realise the multivariate normal it describes and hand that to the
same corner machinery everything else uses. That also makes a Fisher result and
a sampled DALI posterior directly overlayable, which is the point.
"""

import logging

import numpy
import pylab

from ..utils import GWLATEX_LABELS

#: Contour levels at one, two and three sigma in two dimensions.
CONTOUR_LEVELS = (
    1 - numpy.exp(-0.5),
    1 - numpy.exp(-2),
    1 - numpy.exp(-9 / 2.0),
)

#: Colour cycle for overlaid datasets, colourblind-safe and distinguishable in
#: greyscale.
OVERLAY_COLOURS = ("#ca0020", "#0571b0", "#4daf4a", "#984ea3")


def samples_from_covariance(mean, covariance, size=200000, seed=None):
    """Draw from the multivariate normal a Fisher matrix describes.

    Parameters
    ----------
    mean : array_like
        Fiducial parameter values.
    covariance : array_like
        ``(n, n)`` covariance matrix.
    size : int
        Number of samples.
    seed : int or None
        Seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        ``(size, n)`` samples.

    Notes
    -----
    A Fisher covariance can pick up tiny negative eigenvalues from rounding when
    the matrix is poorly conditioned, which makes the normal undrawable. Those
    are clipped to zero with a warning; a genuinely indefinite covariance means
    the forecast itself is not trustworthy.
    """
    mean = numpy.asarray(mean, dtype=float)
    covariance = numpy.asarray(covariance, dtype=float)
    covariance = 0.5 * (covariance + covariance.T)

    eigenvalues, eigenvectors = numpy.linalg.eigh(covariance)
    if numpy.any(eigenvalues <= 0):
        smallest = eigenvalues.min()
        logging.warning(
            "Covariance is not positive definite (smallest eigenvalue %.3e); "
            "clipping to draw samples, but the forecast is unreliable.",
            smallest,
        )
        eigenvalues = numpy.clip(eigenvalues, 0.0, None)
        covariance = eigenvectors @ numpy.diag(eigenvalues) @ eigenvectors.T

    generator = numpy.random.default_rng(seed)
    return generator.multivariate_normal(mean, covariance, size=size, method="eigh")


def labels_for(parameters):
    """LaTeX axis labels, falling back to the raw name when none is known.

    The fallback escapes underscores so matplotlib's mathtext does not read them
    as subscripts and silently mangle the label.
    """
    return [
        GWLATEX_LABELS.get(name, name.replace("_", r"\_")) for name in parameters
    ]


def overlay_corner(datasets, parameters, truths=None, save=None, title=None):
    """Overlay several sample sets on one corner plot.

    Parameters
    ----------
    datasets : sequence of tuple
        ``(label, samples)`` pairs, each ``samples`` of shape
        ``(n_samples, n_parameters)`` with columns ordered as ``parameters``.
    parameters : sequence of str
        Parameter names, used for the axis labels.
    truths : array_like or None
        Reference values to mark.
    save : str or None
        Path to write the figure to. When None the figure is returned unsaved.
    title : str or None
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If ``datasets`` is empty or a sample array has the wrong width.
    """
    from corner import corner

    if not datasets:
        raise ValueError("overlay_corner needs at least one dataset")
    count = len(parameters)
    for label, samples in datasets:
        if numpy.shape(samples)[1] != count:
            raise ValueError(
                "dataset '{}' has {} columns but {} parameters were named".format(
                    label, numpy.shape(samples)[1], count
                )
            )

    figure = pylab.figure(figsize=(3 * count, 3 * count))
    options = dict(
        bins=50,
        smooth=0.9,
        labels=labels_for(parameters),
        label_kwargs=dict(fontsize=14),
        title_kwargs=dict(fontsize=14),
        levels=CONTOUR_LEVELS,
        plot_density=False,
        plot_datapoints=False,
        fill_contours=False,
        max_n_ticks=3,
        truths=None if truths is None else numpy.asarray(truths, dtype=float),
        truth_color="black",
    )

    single = len(datasets) == 1
    if single:
        options.update(
            show_titles=True,
            quantiles=[0.16, 0.84],
            title_quantiles=[0.16, 0.5, 0.84],
        )

    for index, (label, samples) in enumerate(datasets):
        corner(
            numpy.asarray(samples),
            fig=figure,
            color=OVERLAY_COLOURS[index % len(OVERLAY_COLOURS)],
            **options,
        )

    handles = [
        pylab.Line2D([], [], color=OVERLAY_COLOURS[index % len(OVERLAY_COLOURS)], label=label)
        for index, (label, _) in enumerate(datasets)
    ]
    figure.legend(handles=handles, loc="upper right", fontsize=2 + 6 * count // 3, frameon=False)
    if title:
        figure.suptitle(title, fontsize=8 + 4 * count // 3)

    if save:
        logging.info("Saving corner plot to %s", save)
        figure.savefig(save, bbox_inches="tight", dpi=100)
    return figure
