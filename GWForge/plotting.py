r"""Corner plots and figure style for GWForge.
"""

import logging

import numpy
import pylab
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter

from .utils import GWLATEX_LABELS

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Golden ratio, the default figure aspect.
PHI = 1.618033988749895

# CARTOColors TealRose, 7-class diverging. Transcribed from
# ``palettable.cartocolors.colormaps._TEALROSE``.
TEALROSE = (
    "#009392",
    "#72aaa1",
    "#b1c7b3",
    "#f1eac8",
    "#e5b9ad",
    "#d98994",
    "#d0587e",
)

# Okabe & Ito (2008) colour-universal-design qualitative set. Black is last so
# that the first two datasets in an overlay come out orange and sky blue.
OKABE_ITO = (
    "#e69f00",
    "#56b4e9",
    "#009e73",
    "#f0e442",
    "#0072b2",
    "#d55e00",
    "#cc79a7",
    "#000000",
)

_OKABE_ITO_RAMP = ("#0072b2",
                   "#56b4e9",
                   "#009e73",
                   "#f0e442",
                   "#e69f00",
                   "#d55e00")

# Default palette name.
DEFAULT_PALETTE = "tealrose"

XG_COLOUR = TEALROSE[0]
LVK_COLOUR = TEALROSE[-1]

# Grey for truth lines, matching the reference figures.
TRUTH_COLOUR = "0.2"

DEFAULT_LEVELS = (0.5, 0.9)
DEFAULT_QUANTILES = (0.16, 0.5, 0.84)
DEFAULT_SPAN_SIGMA = 4.0

HISTOGRAM_BINS = 60
SMOOTHING_BINS = 1.0

POPULATION_LATEX_LABELS = dict(GWLATEX_LABELS)
POPULATION_LATEX_LABELS.update(
    {
        # BGP mass model
        "alpha_1": r"$\alpha_{1}$",
        "alpha_2": r"$\alpha_{2}$",
        "m_break": r"$m_{\mathrm{break}}\ [M_{\odot}]$",
        "mmin": r"$m_{\mathrm{min}}\ [M_{\odot}]$",
        "m_high": r"$m_{\mathrm{high}}\ [M_{\odot}]$",
        "maximum_mass": r"$m_{\mathrm{max}}\ [M_{\odot}]$",
        "lam_0": r"$\lambda_{0}$",
        "lam_1": r"$\lambda_{1}$",
        "mpp_1": r"$\mu_{1}\ [M_{\odot}]$",
        "sigpp_1": r"$\sigma_{1}\ [M_{\odot}]$",
        "mpp_2": r"$\mu_{2}\ [M_{\odot}]$",
        "sigpp_2": r"$\sigma_{2}\ [M_{\odot}]$",
        "delta_m": r"$\delta_{m}\ [M_{\odot}]$",
        "beta": r"$\beta_{q}$",
        # Madau-Dickinson merger rate
        "gamma": r"$\gamma$",
        "kappa": r"$\kappa$",
        "z_peak": r"$z_{\mathrm{p}}$",
        # Default spin model
        "mu_chi": r"$\mu_{\chi}$",
        "sigma_chi": r"$\sigma_{\chi}$",
        "mu_t": r"$\mu_{t}$",
        "sigma_t": r"$\sigma_{t}$",
        "xi_spin": r"$\xi$",
        # The Beta models still carry a variance.
        "sigma_squared_chi": r"$\sigma^{2}_{\chi}$",
        # BGP's independent secondary taper.
        "mmin_2": r"$m_{2,\mathrm{low}}\ [M_{\odot}]$",
        "delta_m_2": r"$\delta_{m,2}\ [M_{\odot}]$",
        # Cosmology
        "H0": r"$H_{0}\ [\mathrm{km/s/Mpc}]$",
        "Om0": r"$\Omega_{m,0}$",
        "w0": r"$w_{0}$",
        # Event parameters, for catalogue scatter plots
        "redshift": r"$z$",
        "mass_1_source": r"$m_{1}^{\mathrm{src}}\ [M_{\odot}]$",
        "mass_2_source": r"$m_{2}^{\mathrm{src}}\ [M_{\odot}]$",
        "network_snr": r"$\rho_{\mathrm{net}}$",
    }
)


def new_rcParams(width="column", aspect_ratio=PHI):
    """Publication ready plot
    Parameters
    ----------
    width : str
        ``"column"`` (246 pt) or ``"page"`` (510 pt), each doubled so the saved
        figure has room for readable fonts when scaled down.
    aspect_ratio : float
        Width divided by height. Defaults to the golden ratio.

    Returns
    -------
    dict
        Suitable for ``matplotlib.rc_context`` -- which is how it should be
        used; see the note in the module docstring about
        :mod:`GWForge.utils` restyling matplotlib on import.
    """
    scale_factor = 2
    if width == "column":
        figure_width_points = scale_factor * 246.0
        font_size = scale_factor * 7.96
    elif width == "page":
        figure_width_points = scale_factor * 510.0
        font_size = scale_factor * 9.0
    else:
        raise ValueError("width must be 'column' or 'page', got {!r}".format(width))

    inches_per_point = 1.0 / 72.27
    figure_width = figure_width_points * inches_per_point
    return {
        "figure.figsize": (figure_width, figure_width / aspect_ratio),
        "font.size": font_size,
        "text.usetex": False,
        "font.family": "stixgeneral",
        "mathtext.fontset": "stix",
        "axes.labelsize": "medium",
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "legend.fontsize": "medium",
        "legend.handlelength": 1.5,
        "grid.linestyle": "--",
        "grid.color": "#bbbbbb",
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "savefig.format": "pdf",
    }


def palette(name=DEFAULT_PALETTE):
    """Discrete colour cycle and continuous colormap for a named palette.

    Parameters
    ----------
    name : str
        ``"tealrose"`` or ``"okabe-ito"``.

    Returns
    -------
    tuple
        ``(colours, colormap)``. ``colours`` is a tuple of hex strings to cycle
        through for overlaid datasets; ``colormap`` is a
        ``matplotlib.colors.LinearSegmentedColormap`` for
        :func:`corner_plot`'s ``colour_by`` scatter.

    Raises
    ------
    ValueError
        For an unknown palette name.
    """
    key = str(name).lower().replace("_", "-")
    if key == "tealrose":
        colours = (
            TEALROSE[0],
            TEALROSE[6],
            TEALROSE[1],
            TEALROSE[5],
            TEALROSE[2],
            TEALROSE[4],
            TEALROSE[3],
        )
        colormap = LinearSegmentedColormap.from_list("tealrose", TEALROSE, N=256)
    elif key == "okabe-ito":
        colours = OKABE_ITO
        colormap = LinearSegmentedColormap.from_list(
            "okabe-ito", _OKABE_ITO_RAMP, N=256
        )
    else:
        raise ValueError(
            "Unknown palette {!r}; choose 'tealrose' or 'okabe-ito'.".format(name)
        )
    return colours, colormap


def labels_for(parameters):
    """LaTeX labels for parameter names, falling back to the escaped raw name.

    The fallback escapes underscores so mathtext does not read them as
    subscripts and silently mangle the label.

    Parameters
    ----------
    parameters : sequence of str

    Returns
    -------
    list of str
    """
    return [
        POPULATION_LATEX_LABELS.get(name, name.replace("_", r"\_"))
        for name in parameters
    ]


def _quantile(values, quantiles, weights=None):
    """Weighted quantiles, matching zeus's ``_quantile``.

    ``numpy.quantile`` has no weights argument, so the weighted case sorts,
    builds the cumulative weight as a CDF and interpolates.

    Parameters
    ----------
    values : array_like
    quantiles : array_like
        In ``[0, 1]``.
    weights : array_like or None

    Returns
    -------
    numpy.ndarray
    """
    values = numpy.asarray(values, dtype=float)
    quantiles = numpy.asarray(quantiles, dtype=float)
    if weights is None:
        return numpy.percentile(values, 100.0 * quantiles)
    weights = numpy.asarray(weights, dtype=float)
    order = numpy.argsort(values)
    cumulative = numpy.cumsum(weights[order])
    cumulative /= cumulative[-1]
    return numpy.interp(quantiles, cumulative, values[order])


def enclosed_levels(density, levels):
    """Density values enclosing the requested fractions of the total probability.

    Sort the histogram, walk down the cumulative sum, and read off where each
    fraction is reached. This is what makes a "90% contour" contain 90% of the
    samples rather than trace an arbitrary isopleth of the density.

    Parameters
    ----------
    density : numpy.ndarray
        Histogram or smoothed density, any shape.
    levels : sequence of float
        Enclosed probabilities in ``(0, 1)``.

    Returns
    -------
    numpy.ndarray
        Contour values, sorted ascending as ``matplotlib.contour`` requires.
    """
    flat = numpy.sort(numpy.asarray(density, dtype=float).ravel())[::-1]
    cumulative = numpy.cumsum(flat)
    cumulative /= cumulative[-1]
    values = []
    for level in levels:
        index = numpy.searchsorted(cumulative, level)
        values.append(flat[min(index, len(flat) - 1)])
    # Duplicates arise when two levels fall inside one histogram step; nudge
    # them apart so contour() does not drop one silently.
    values = numpy.sort(numpy.unique(values))
    return values


def _panel_density(x_values, y_values, weights, span_x, span_y, bins, smooth):
    """Smoothed 2-D density on a fixed span, plus its bin centres."""
    density, x_edges, y_edges = numpy.histogram2d(
        x_values, y_values, bins=bins, range=[span_x, span_y], weights=weights
    )
    if smooth:
        density = gaussian_filter(density, smooth)
    return (
        density.T,
        0.5 * (x_edges[1:] + x_edges[:-1]),
        0.5 * (y_edges[1:] + y_edges[:-1]),
    )


def _default_span(samples, sigma=DEFAULT_SPAN_SIGMA):
    """Mean +/- ``sigma`` sample standard deviations, per column."""
    samples = numpy.atleast_2d(samples)
    centre = samples.mean(axis=0)
    width = sigma * samples.std(axis=0)
    return [(centre[i] - width[i], centre[i] + width[i]) for i in range(len(centre))]


def corner_plot(
    samples,
    labels=None,
    weights=None,
    levels=DEFAULT_LEVELS,
    span=None,
    quantiles=DEFAULT_QUANTILES,
    truths=None,
    colour=None,
    alpha=0.5,
    linewidth=1.5,
    fill=None,
    colour_by=None,
    colour_label=None,
    palette_name=DEFAULT_PALETTE,
    scatter_size=2.0,
    show_titles=True,
    title_fmt=".3f",
    label=None,
    fig=None,
    size=None,
    bins=HISTOGRAM_BINS,
    smooth=SMOOTHING_BINS,
):
    r"""Corner plot: 1-D marginals on the diagonal, 2-D contours below it.

    Call once per dataset, passing the returned figure back as ``fig`` to
    overlay a second. Pass the same ``span`` to both so the panels line up --
    ranges are computed per call otherwise.

    Parameters
    ----------
    samples : array_like
        ``(n_samples, n_dimensions)``.
    labels : sequence of str or None
        Axis labels, already LaTeX. See :func:`labels_for`.
    weights : array_like or None
        Per-sample weights.
    levels : sequence of float
        **Enclosed probabilities**, not density values. See
        :func:`enclosed_levels`.
    span : sequence of tuple or None
        ``(low, high)`` per dimension. Defaults to mean +/- 4 sample standard
        deviations.
    quantiles : sequence of float
        Three quantiles for the diagonal titles: lower, central, upper.
    truths : array_like or None
        Reference values, drawn as lines.
    colour : str or None
        Contour colour. Defaults to the first entry of ``palette_name``.
    alpha : float
        Fill opacity. The outline is always drawn at full opacity so overlaid
        contours stay legible.
    linewidth : float
    fill : bool or None
        Fill the contours as well as outlining them. ``None`` means *yes,
        unless* ``colour_by`` is given -- a filled contour drawn over a scatter
        hides exactly the points the scatter exists to show.
    colour_by : array_like or None
        Length ``n_samples``. When given, each lower-triangle panel scatters the
        samples coloured by this on the palette's continuous map, with one
        shared colourbar. This is the pycbc posterior idiom -- events coloured
        by redshift or by network SNR.
    colour_label : str or None
        Colourbar label.
    palette_name : str
        ``"tealrose"`` or ``"okabe-ito"``.
    scatter_size : float
        Marker size for the ``colour_by`` scatter.
    show_titles : bool
        Quantile titles over the diagonal panels. Pass ``False`` on overlays so
        the first dataset's titles are not overwritten.
    title_fmt : str
    label : str or None
        Legend entry for this dataset.
    fig : matplotlib.figure.Figure or None
        Existing corner to draw onto.
    size : tuple or None
        Figure size in inches. Defaults to ``(3 n, 3 n)``.
    bins : int
        Histogram bins per axis.
    smooth : float
        Gaussian smoothing in units of bins. ``0`` disables it.

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If ``samples`` is not two-dimensional, or if ``labels``, ``truths``,
        ``span`` or ``colour_by`` do not match its shape, or if ``fig`` has a
        different number of panels.
    """
    samples = numpy.atleast_2d(numpy.asarray(samples, dtype=float))
    if samples.ndim != 2:
        raise ValueError(
            "samples must be (n_samples, n_dimensions), got shape {}".format(
                samples.shape
            )
        )
    n_samples, n_dimensions = samples.shape

    for name, value, expected in (
        ("labels", labels, n_dimensions),
        ("truths", truths, n_dimensions),
        ("span", span, n_dimensions),
        ("colour_by", colour_by, n_samples),
        ("weights", weights, n_samples),
    ):
        if value is not None and len(value) != expected:
            raise ValueError(
                "{} has length {} but {} were expected.".format(
                    name, len(value), expected
                )
            )

    colours, colormap = palette(palette_name)
    if colour is None:
        colour = colours[0]
    if fill is None:
        fill = colour_by is None
    if labels is None:
        labels = ["$x_{{{}}}$".format(i) for i in range(n_dimensions)]
    if span is None:
        span = _default_span(samples)

    style = new_rcParams(width="page")
    if size is None:
        size = (3.0 * n_dimensions, 3.0 * n_dimensions)

    with pylab.matplotlib.rc_context(style):
        if fig is None:
            figure, axes = pylab.subplots(
                n_dimensions, n_dimensions, figsize=size, squeeze=False
            )
        else:
            figure = fig
            existing = numpy.array(figure.axes)
            # A colourbar adds an axes that is not part of the grid.
            grid = [ax for ax in existing if getattr(ax, "_gwforge_panel", False)]
            if len(grid) != n_dimensions**2:
                raise ValueError(
                    "The supplied figure has {} corner panels but {} are needed "
                    "for {} dimensions.".format(
                        len(grid), n_dimensions**2, n_dimensions
                    )
                )
            axes = numpy.array(grid).reshape(n_dimensions, n_dimensions)

        scatter = None
        for row in range(n_dimensions):
            for column in range(n_dimensions):
                axis = axes[row, column]
                axis._gwforge_panel = True
                if column > row:
                    axis.set_visible(False)
                    continue

                if row == column:
                    _draw_marginal(
                        axis,
                        samples[:, row],
                        weights,
                        span[row],
                        colour,
                        alpha,
                        linewidth,
                        fill,
                        bins,
                        smooth,
                        quantiles if show_titles else None,
                        title_fmt,
                        labels[row],
                    )
                else:
                    handle = _draw_panel(
                        axis,
                        samples[:, column],
                        samples[:, row],
                        weights,
                        span[column],
                        span[row],
                        levels,
                        colour,
                        alpha,
                        linewidth,
                        fill,
                        bins,
                        smooth,
                        colour_by,
                        colormap,
                        scatter_size,
                    )
                    scatter = scatter or handle

                if truths is not None:
                    axis.axvline(
                        truths[column], color=TRUTH_COLOUR, linewidth=1.0, zorder=5
                    )
                    if row != column:
                        axis.axhline(
                            truths[row], color=TRUTH_COLOUR, linewidth=1.0, zorder=5
                        )

                axis.set_xlim(span[column])
                if row != column:
                    axis.set_ylim(span[row])
                if row == n_dimensions - 1:
                    axis.set_xlabel(labels[column])
                else:
                    axis.set_xticklabels([])
                if column == 0 and row != 0:
                    axis.set_ylabel(labels[row])
                elif row != column:
                    axis.set_yticklabels([])
                # Prune the end ticks: adjacent panels share an edge, so an
                # unpruned first/last label collides with its neighbour's.
                axis.xaxis.set_major_locator(
                    MaxNLocator(nbins=4, prune="both")
                )
                if row != column:
                    axis.yaxis.set_major_locator(
                        MaxNLocator(nbins=4, prune="both")
                    )

        if scatter is not None and colour_by is not None:
            colourbar = figure.colorbar(
                scatter, ax=axes.ravel().tolist(), fraction=0.03, pad=0.02
            )
            if colour_label:
                colourbar.set_label(colour_label)

        _update_legend(figure, colour, label)
        if fig is None:
            figure.subplots_adjust(wspace=0.06, hspace=0.06)
    return figure


def _draw_marginal(
    axis,
    values,
    weights,
    span,
    colour,
    alpha,
    linewidth,
    fill,
    bins,
    smooth,
    quantiles,
    title_fmt,
    label,
):
    """One diagonal panel: the smoothed 1-D marginal, optionally titled."""
    density, edges = numpy.histogram(
        values, bins=bins, range=span, weights=weights, density=True
    )
    if smooth:
        density = gaussian_filter(density, smooth)
    centres = 0.5 * (edges[1:] + edges[:-1])
    axis.plot(centres, density, color=colour, linewidth=linewidth)
    if fill:
        axis.fill_between(centres, density, color=colour, alpha=alpha, linewidth=0)
    axis.set_ylim(bottom=0.0)
    axis.set_yticks([])

    if quantiles is not None:
        low, middle, high = _quantile(values, quantiles, weights)
        axis.set_title(
            "${}^{{+{}}}_{{-{}}}$".format(
                format(middle, title_fmt),
                format(high - middle, title_fmt),
                format(middle - low, title_fmt),
            ),
            fontsize="small",
        )


def _draw_panel(
    axis,
    x_values,
    y_values,
    weights,
    span_x,
    span_y,
    levels,
    colour,
    alpha,
    linewidth,
    fill,
    bins,
    smooth,
    colour_by,
    colormap,
    scatter_size,
):
    """One lower-triangle panel: optional scatter, then contours over it."""
    handle = None
    if colour_by is not None:
        handle = axis.scatter(
            x_values,
            y_values,
            c=colour_by,
            cmap=colormap,
            s=scatter_size,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )

    density, x_centres, y_centres = _panel_density(
        x_values, y_values, weights, span_x, span_y, bins, smooth
    )
    contours = enclosed_levels(density, levels)
    if fill:
        axis.contourf(
            x_centres,
            y_centres,
            density,
            levels=numpy.concatenate([contours, [density.max() + 1.0]]),
            colors=colour,
            alpha=alpha,
            zorder=2,
        )
    axis.contour(
        x_centres,
        y_centres,
        density,
        levels=contours,
        colors=colour,
        linewidths=linewidth,
        zorder=3,
    )
    return handle


def _update_legend(figure, colour, label):
    """Append one proxy handle to the figure legend, keeping earlier entries."""
    if not label:
        return
    from matplotlib.lines import Line2D

    handles = list(getattr(figure, "_gwforge_handles", []))
    handles.append(Line2D([], [], color=colour, linewidth=2.0, label=label))
    figure._gwforge_handles = handles
    if figure.legends:
        figure.legends[-1].remove()
    figure.legend(handles=handles, loc="upper right", frameon=False)
