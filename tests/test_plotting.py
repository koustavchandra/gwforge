"""Tests for :mod:`GWForge.plotting`.

A plotting test that only checks "it did not raise" is nearly worthless, so
these check the two things that can be quietly *wrong* rather than broken.

**The contour levels.** ``corner_plot`` advertises its levels as enclosed
probabilities: a 90% contour should contain 90% of the samples. That is a claim
about :func:`~GWForge.plotting.enclosed_levels`, and it is checked by counting
the samples inside the drawn path — not by asserting that a contour exists.

**The style is not global.** The figures apply their rcParams through
``matplotlib.rc_context``. The reference implementation this is modelled on does
a bare ``rcParams.update`` inside its corner function, which silently restyles
everything the caller draws afterwards. That is the bug these tests pin shut.
"""

import matplotlib

matplotlib.use("Agg")

import numpy
import pylab
import pytest
from matplotlib.path import Path

from GWForge.plotting import (
    LVK_COLOUR,
    OKABE_ITO,
    TEALROSE,
    XG_COLOUR,
    corner_plot,
    enclosed_levels,
    labels_for,
    new_rcParams,
    palette,
)
from GWForge.population_fisher import (
    BrokenPowerLawTwoPeakMass,
    DefaultSpin,
    MadauDickinsonRedshift,
)
from GWForge.plotting import POPULATION_LATEX_LABELS

COVARIANCE = numpy.array([[1.0, 0.7, 0.2], [0.7, 2.0, -0.3], [0.2, -0.3, 0.5]])
CENTRE = [0.0, 1.0, 2.0]


@pytest.fixture
def samples():
    generator = numpy.random.default_rng(20260822)
    return generator.multivariate_normal(CENTRE, COVARIANCE, 30000)


def panels(figure):
    """The corner grid, excluding any colourbar axes."""
    return [axis for axis in figure.axes if getattr(axis, "_gwforge_panel", False)]


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def test_corner_has_a_lower_triangular_grid(samples):
    figure = corner_plot(samples, labels=["$a$", "$b$", "$c$"])
    grid = numpy.array(panels(figure)).reshape(3, 3)
    assert len(panels(figure)) == 9
    for row in range(3):
        for column in range(3):
            assert grid[row, column].get_visible() == (column <= row)
    assert grid[2, 0].get_xlabel() == "$a$"
    assert grid[2, 2].get_xlabel() == "$c$"
    assert grid[1, 0].get_ylabel() == "$b$"
    pylab.close(figure)


def test_truths_are_drawn_where_asked(samples):
    figure = corner_plot(samples, truths=CENTRE)
    grid = numpy.array(panels(figure)).reshape(3, 3)
    vertical = [line.get_xdata()[0] for line in grid[1, 0].lines]
    horizontal = [line.get_ydata()[0] for line in grid[1, 0].lines]
    assert CENTRE[0] in vertical
    assert CENTRE[1] in horizontal
    pylab.close(figure)


def test_overlay_reuses_the_same_panels(samples):
    """A second dataset must draw *onto* the first, not beside it."""
    generator = numpy.random.default_rng(1)
    wider = generator.multivariate_normal(CENTRE, 4.0 * COVARIANCE, 30000)
    span = [(-8.0, 8.0), (-7.0, 9.0), (-1.0, 5.0)]
    figure = corner_plot(samples, span=span, colour=XG_COLOUR, label="XG")
    before = len(panels(figure))
    figure = corner_plot(
        wider, span=span, colour=LVK_COLOUR, label="LVK", show_titles=False, fig=figure
    )
    assert len(panels(figure)) == before
    grid = numpy.array(panels(figure)).reshape(3, 3)
    assert grid[1, 0].get_xlim() == pytest.approx(span[0])
    assert grid[1, 0].get_ylim() == pytest.approx(span[1])
    # Both datasets contributed contours to the same panel.
    contour_sets = [
        child
        for child in grid[1, 0].get_children()
        if isinstance(child, matplotlib.contour.ContourSet)
    ]
    assert len(contour_sets) >= 4
    pylab.close(figure)


def test_overlay_onto_a_mismatched_figure_is_refused(samples):
    figure = corner_plot(samples[:, :2])
    with pytest.raises(ValueError, match="corner panels"):
        corner_plot(samples, fig=figure)
    pylab.close(figure)


@pytest.mark.parametrize(
    "name,length",
    [("labels", 3), ("truths", 3), ("span", 3), ("colour_by", 30000)],
)
def test_mismatched_lengths_are_refused(samples, name, length):
    with pytest.raises(ValueError, match=name):
        corner_plot(samples, **{name: numpy.zeros(length - 1)})


# ---------------------------------------------------------------------------
# The contour levels mean what they say
# ---------------------------------------------------------------------------


def test_contours_enclose_the_stated_probability():
    """The claim the whole corner rests on, checked by counting samples.

    ``levels`` are enclosed probabilities, not density values: the level is
    found by sorting the smoothed histogram and walking down its cumulative
    sum. If that were wrong the contours would still look plausible, which is
    exactly why this is measured rather than eyeballed.
    """
    generator = numpy.random.default_rng(3)
    draws = generator.multivariate_normal(
        [0.0, 0.0], [[1.0, 0.5], [0.5, 1.0]], 200000
    )
    figure = corner_plot(draws, levels=(0.5, 0.9), bins=80)
    panel = numpy.array(panels(figure)).reshape(2, 2)[1, 0]

    enclosed = []
    for child in panel.get_children():
        if isinstance(child, matplotlib.contour.ContourSet) and not child.filled:
            for path in child.get_paths():
                if len(path.vertices) > 4:
                    enclosed.append(Path(path.vertices).contains_points(draws).mean())
    enclosed.sort()
    assert enclosed[0] == pytest.approx(0.5, abs=0.03)
    assert enclosed[1] == pytest.approx(0.9, abs=0.03)
    pylab.close(figure)


def test_enclosed_levels_are_ascending_and_unique():
    """``matplotlib.contour`` needs sorted, distinct levels or it drops one."""
    density = numpy.exp(-0.5 * numpy.linspace(-4, 4, 200) ** 2)
    values = enclosed_levels(numpy.outer(density, density), (0.5, 0.9, 0.99))
    assert numpy.all(numpy.diff(values) > 0)


# ---------------------------------------------------------------------------
# The colour-by-a-third-quantity mode
# ---------------------------------------------------------------------------


def test_colour_by_adds_a_scatter_and_a_colourbar(samples):
    """The pycbc idiom, and the reason this corner is written rather than configured."""
    third = samples[:, 0] + samples[:, 1]
    figure = corner_plot(
        samples, colour_by=third, colour_label=r"$\rho$", scatter_size=1.0
    )
    extra = [axis for axis in figure.axes if not getattr(axis, "_gwforge_panel", False)]
    assert len(extra) == 1, "expected exactly one colourbar axes"
    assert extra[0].get_ylabel() == r"$\rho$"

    panel = numpy.array(panels(figure)).reshape(3, 3)[1, 0]
    scatters = [c for c in panel.collections if c.get_array() is not None]
    assert scatters, "no colour-mapped scatter was drawn"
    assert scatters[0].get_cmap().name == "tealrose"
    pylab.close(figure)


def test_colour_by_turns_the_fill_off_by_default(samples):
    """A filled contour over a scatter hides the points the scatter exists to show."""
    third = samples[:, 2]
    panel_with = numpy.array(panels(corner_plot(samples, colour_by=third))).reshape(
        3, 3
    )[1, 0]
    panel_without = numpy.array(panels(corner_plot(samples))).reshape(3, 3)[1, 0]

    def filled(axis):
        return [
            child
            for child in axis.get_children()
            if isinstance(child, matplotlib.contour.ContourSet) and child.filled
        ]

    assert not filled(panel_with)
    assert filled(panel_without)
    pylab.close("all")


# ---------------------------------------------------------------------------
# Palettes and labels
# ---------------------------------------------------------------------------


def test_network_colours_come_from_the_palette():
    """So the discrete comparison colours and the continuous map cannot drift apart."""
    assert XG_COLOUR == TEALROSE[0] == "#009392"
    assert LVK_COLOUR == TEALROSE[-1] == "#d0587e"
    colours, _ = palette("tealrose")
    assert colours[:2] == (XG_COLOUR, LVK_COLOUR)


@pytest.mark.parametrize("name", ("tealrose", "okabe-ito", "Okabe_Ito", "TealRose"))
def test_palettes_are_valid_and_continuous(name):
    colours, colormap = palette(name)
    for colour in colours:
        assert matplotlib.colors.is_color_like(colour)
    values = colormap(numpy.linspace(0.0, 1.0, 64))
    assert values.shape == (64, 4)
    assert numpy.all(numpy.isfinite(values))


def test_unknown_palette_is_refused():
    with pytest.raises(ValueError, match="Unknown palette"):
        palette("viridis")


def test_okabe_ito_puts_black_last():
    """Black first would make the default two-dataset overlay black-on-orange."""
    assert OKABE_ITO[-1] == "#000000"
    assert palette("okabe-ito")[0][0] == "#e69f00"


def test_every_hyper_parameter_has_a_label():
    """A new parameter without a label prints a raw key with a mangled underscore."""
    models = [
        BrokenPowerLawTwoPeakMass(),
        MadauDickinsonRedshift(maximum_redshift=10.0),
        DefaultSpin(),
    ]
    missing = [
        name
        for model in models
        for name in model.parameter_names
        if name not in POPULATION_LATEX_LABELS
    ]
    assert not missing, "no LaTeX label for {}".format(missing)


def test_labels_for_escapes_the_fallback():
    assert labels_for(["alpha_1"]) == [r"$\alpha_{1}$"]
    assert labels_for(["not_a_parameter"]) == [r"not\_a\_parameter"]


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------


def test_rcparams_are_not_leaked_globally(samples):
    """The reference's ``plot_corner`` does a bare ``rcParams.update``; this must not.

    Drawing a figure has to leave the caller's matplotlib configuration exactly
    as it found it, or every plot made after a corner silently changes font.
    """
    before = dict(pylab.rcParams)
    figure = corner_plot(samples, labels=["$a$", "$b$", "$c$"])
    pylab.close(figure)
    after = dict(pylab.rcParams)
    changed = [key for key in before if before[key] != after.get(key)]
    assert not changed, "corner_plot leaked {}".format(changed)


def test_new_rcparams_widths():
    column = new_rcParams(width="column")
    page = new_rcParams(width="page")
    assert page["figure.figsize"][0] > column["figure.figsize"][0]
    # Golden ratio, and TeX points rather than PostScript points.
    width, height = column["figure.figsize"]
    assert width / height == pytest.approx(1.618033988749895)
    assert width == pytest.approx(2 * 246.0 / 72.27)
    assert column["mathtext.fontset"] == "stix"
    with pytest.raises(ValueError, match="must be 'column' or 'page'"):
        new_rcParams(width="banner")
