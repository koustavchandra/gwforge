#!/usr/bin/env python
"""Check a generated population against the GWTC-5.0 models it was drawn from.

Run after ``gwforge_population`` to confirm that what came out of the sampler is
what the model says. Every panel compares *samples* with the *density*, so a
disagreement means the two implementations have drifted apart -- which is
exactly how the mass-ratio, taper and spin bugs were found.

    python validation/gwtc5_population_check.py --population bbh.h5 \
        --config GWForge/population/population_configuration_files/bgp-gwtc5.ini

The corner plot is the headline, but the numeric checks below it are what
actually pass or fail. Two of them cannot be seen in any one-dimensional
histogram:

* ``p(q | m_1)`` is a *conditional*, so it is checked in bins of ``m_1``. Drawing
  ``q`` from the marginal instead -- as this package once did -- reproduces the
  marginal perfectly and the conditional not at all.
* the cos-tilt mixture is joint over the binary, so it is checked with the
  correlation between the two tilts. Factorising it leaves both marginals
  untouched and sets that correlation to zero.
"""
import argparse
import configparser
import json
import logging

import matplotlib

matplotlib.use("Agg")

import matplotlib.lines
import numpy
import pylab
from scipy import stats
from scipy.integrate import cumulative_trapezoid

from GWForge.population._smoothed_mass import (
    BrokenPowerLawTwoPeakSmoothedMassDistribution,
)
from GWForge.population.extrinsic import Extrinsic
from GWForge.population.mass import BGP_PARAMETERS, primary_mass_grid_nodes
from GWForge.population.redshift import Redshift
from GWForge.population.spin import (
    DEFAULT_BBH_SPIN_PARAMETERS,
    default_spin_magnitude_density,
    default_spin_tilt_density,
)
from GWForge.plotting import corner_plot, labels_for, new_rcParams, palette

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Parameters the corner shows, in order: everything ``gwforge_population``
# draws from a distribution, the two spin components kept separate so the panel
# shows the secondary as well as the primary.
CORNER_PARAMETERS = [
    "mass_1_source",
    "mass_2_source",
    "redshift",
    "a_1",
    "a_2",
    "cos_tilt_1",
    "cos_tilt_2",
    "ra",
    "dec",
    "theta_jn",
    "psi",
]

# The extrinsic parameters, and the range each is defined on. These are bilby
# priors rather than population models -- an isotropic sky and orientation --
# but they are drawn and they are in the figure, so they are checked too.
EXTRINSIC_SPANS = {
    "ra": (0.0, 2.0 * numpy.pi),
    "dec": (-0.5 * numpy.pi, 0.5 * numpy.pi),
    "theta_jn": (0.0, numpy.pi),
    "psi": (0.0, numpy.pi),
}

# Primary masses to integrate over when marginalising p(m_1) p(q | m_1) down to
# p(m_2). Coarser than the grid p(m_1) itself is drawn on, because the mesh is
# two-dimensional and the integrand is smooth in m_1.
SECONDARY_MARGINAL_NODES = 1200

# Bins of primary mass in which p(q | m_1) is checked separately.
PRIMARY_MASS_BINS = ((8.0, 12.0), (20.0, 30.0), (40.0, 60.0))


def read_config(path):
    """Mass, spin and redshift settings from a gwforge_population ini."""
    config = configparser.ConfigParser()
    if not config.read(path):
        raise FileNotFoundError(path)

    def collection(section, option, fallback):
        """``fallback`` when the ini omits the option, as the executable does.

        A configuration may name the model and leave the parameters out, in
        which case it gets the GWTC-5.0 medians the model modules carry.
        """
        if not config.has_option(section, option):
            return dict(fallback)
        return json.loads(config.get(section, option).replace("'", '"'))

    return {
        "mass": collection("Mass", "mass-parameters", BGP_PARAMETERS),
        "spin": collection("Spin", "spin-parameters", DEFAULT_BBH_SPIN_PARAMETERS),
        "redshift": collection(
            "Redshift",
            "redshift-parameters",
            {"gamma": 2.7, "kappa": 5.6, "z_peak": 1.9},
        ),
        "inclination_distribution": config.get(
            "Extrinsic", "inclination-distribution", fallback=None
        ),
        "maximum_redshift": config.getfloat("Redshift", "maximum-redshift"),
        "rate": config.getfloat("Redshift", "local-merger-rate-density"),
        "cosmology": config.get("Redshift", "cosmology", fallback="Planck18"),
        "time_delay_model": config.get(
            "Redshift", "time-delay-model", fallback="inverse"
        ),
    }


def mass_model(parameters):
    """The same ``_smoothed_mass`` object ``Mass`` builds, at the same resolution."""
    maximum_mass = parameters.get("maximum_mass", 200)
    nodes = primary_mass_grid_nodes(parameters, parameters["mmin"], maximum_mass)
    return BrokenPowerLawTwoPeakSmoothedMassDistribution(
        mmin=parameters["mmin"],
        mmax=maximum_mass,
        normalization_shape=(nodes, 1000),
    )


def primary_density(model, parameters, grid):
    """Analytic ``pi(m_1)`` on ``grid``."""
    keywords = {
        name: value
        for name, value in parameters.items()
        if name not in ("beta", "maximum_mass", "mmin_2", "delta_m_2")
    }
    return model.p_m1({"mass_1": grid}, **keywords)


def ratio_density_mesh(model, parameters, primary, ratio):
    """Analytic ``p(q | m_1)`` over matching arrays of ``m_1`` and ``q``."""
    return model.p_q(
        {"mass_1": primary, "mass_ratio": ratio},
        beta=parameters["beta"],
        mmin=parameters["mmin"],
        delta_m=parameters["delta_m"],
        mmin_2=parameters.get("mmin_2"),
        delta_m_2=parameters.get("delta_m_2"),
    )


def ratio_density(model, parameters, primary, grid):
    """Analytic ``p(q | m_1)`` at one ``m_1``."""
    return ratio_density_mesh(
        model, parameters, numpy.full_like(grid, primary), grid
    )


def secondary_density(model, parameters, primary_grid, grid):
    """Analytic ``pi(m_2)``, marginalised over the primary.

    The model is written as ``pi(m_1)`` and ``p(q | m_1)``, so the secondary's
    marginal has to be integrated out of the pair:

    .. math::
        \\pi(m_2) = \\int \\pi(m_1)\\, p(m_2 / m_1 \\mid m_1)\\, \\frac{dm_1}{m_1}

    over the primaries heavy enough to admit that secondary. ``p_q`` rebuilds its
    per-``m_1`` normalisation on every call, so this is one vectorised evaluation
    over the whole ``(m_1, m_2)`` mesh rather than a loop over primaries.
    """
    primary_mesh, secondary_mesh = numpy.meshgrid(primary_grid, grid, indexing="ij")
    ratio = secondary_mesh / primary_mesh
    allowed = ratio <= 1.0
    # p_q is only asked about ratios it is defined on; the rest are zeroed after,
    # so the clipped values it returns for them are never used.
    density = ratio_density_mesh(
        model, parameters, primary_mesh, numpy.clip(ratio, 1e-6, 1.0)
    )
    integrand = numpy.where(allowed, density, 0.0)
    integrand *= primary_density(model, parameters, primary_grid)[:, numpy.newaxis]
    integrand /= primary_mesh
    return numpy.trapezoid(integrand, primary_grid, axis=0)


def extrinsic_density(name, grid, inclination_distribution=None):
    """Analytic density for one of the extrinsic priors.

    The defaults in :class:`GWForge.population.extrinsic.Extrinsic` are an
    isotropic sky and orientation, written out here. ``theta_jn`` is the one that
    is configurable, so it follows ``inclination-distribution``.
    """
    if name == "ra":
        return numpy.full_like(grid, 1.0 / (2.0 * numpy.pi))
    if name == "psi":
        return numpy.full_like(grid, 1.0 / numpy.pi)
    if name == "dec":
        return 0.5 * numpy.cos(grid)
    if name != "theta_jn":
        raise ValueError("No extrinsic density for {!r}".format(name))
    if inclination_distribution is None:
        return 0.5 * numpy.sin(grid)
    if "schutz" not in inclination_distribution.lower():
        raise ValueError(
            "Unknown inclination distribution {!r}".format(inclination_distribution)
        )
    # The same sin(theta) Jacobian the sampler applies: schutz_inclination_prob
    # is a density in cos(iota), and this is a density in theta.
    density = Extrinsic(1).schutz_inclination_prob(grid) * numpy.sin(grid)
    return density / numpy.trapezoid(density, grid)


def magnitude_density(spin, grid):
    """Analytic ``p(chi)``, from the model module rather than restated here."""
    return default_spin_magnitude_density(
        grid, spin["mu_chi"], spin["sigma_chi"], amax=spin.get("amax", 1.0)
    )


def cosine_tilt_density(spin, grid):
    """Analytic marginal ``p(cos t)``: the mixture, marginalised over the partner."""
    return default_spin_tilt_density(
        grid,
        spin["mu_t"],
        spin["sigma_t"],
        spin["xi_spin"],
        t_min=spin.get("t_min", -1.0),
    )


def tilt_correlation(spin):
    """Predicted ``corr(cos t_1, cos t_2)`` for the joint mixture.

    Zero if the mixture is factorised per component, which is the whole point of
    checking it.
    """
    low = (-1.0 - spin["mu_t"]) / spin["sigma_t"]
    high = (1.0 - spin["mu_t"]) / spin["sigma_t"]
    truncated = stats.truncnorm(low, high, loc=spin["mu_t"], scale=spin["sigma_t"])
    mean, variance = truncated.mean(), truncated.var()
    xi = spin["xi_spin"]
    first = xi * mean
    second = xi * (variance + mean**2) + (1.0 - xi) / 3.0
    return (xi * mean**2 - first**2) / (second - first**2)


def redshift_density(settings, grid):
    """Analytic ``p(z)``: the time-delayed merger rate, normalised."""
    model = Redshift(
        redshift_model="MadauDickinson",
        local_merger_rate_density=settings["rate"],
        maximum_redshift=settings["maximum_redshift"],
        gps_start_time=0,
        cosmology=settings["cosmology"],
        parameters=settings["redshift"],
        time_delay_model=settings["time_delay_model"],
    )
    density = model.coalescence_rate()(grid)
    return density / numpy.trapezoid(density, grid)


# A histogram bin is compared only if the model predicts at least this many
# counts in it -- below that the Poisson approximation to its error is poor and
# the comparison says nothing.
MINIMUM_EXPECTED_COUNTS = 25

# Bins for the one-dimensional comparisons.
HISTOGRAM_BINS = 60

# Probability that a *correct* sampler fails one check by chance. The threshold
# on the worst residual is derived from this and the number of bins actually
# tested, rather than fixed: the worst of 60 bins is systematically larger than
# the worst of 20, so one flat cut cannot mean the same thing in both. It used
# to be flat at 4 sigma, which is a 0.4% false alarm at 60 bins -- tolerable
# while every check fell to zero at its edges and kept 20-40 bins, but the flat
# extrinsic priors keep all 60, and at thirteen checks a run would cry wolf
# about 5% of the time.
FALSE_ALARM_PROBABILITY = 1e-3


def compare(samples, density, grid, label, bins=HISTOGRAM_BINS):
    """Histogram against density, in units of the histogram's own scatter.

    A fixed percentage tolerance is the wrong test: with ~30,000 events over 60
    bins the counting noise alone is several percent in the bulk and far more in
    the tails, so a 5% threshold flags a perfectly good sampler and a 20% one
    would miss a real bias in the core. Each bin is therefore compared in units
    of its Poisson uncertainty, and the summary is the worst standardised
    residual together with a reduced chi-squared over all usable bins.

    Returns
    -------
    tuple
        ``(worst_standardised_residual, passed)``.
    """
    samples = numpy.asarray(samples)
    # Enough bins to resolve structure, few enough that each holds a testable
    # number of counts. A sparsely populated slice gets a coarser histogram
    # rather than no test at all.
    bins = int(numpy.clip(samples.size // (2 * MINIMUM_EXPECTED_COUNTS), 8, bins))
    edges = numpy.linspace(grid[0], grid[-1], bins + 1)
    counts, _ = numpy.histogram(samples, bins=edges)

    # Expected counts are the density *integrated over each bin*, not its value
    # at the bin centre times the width. The BGP low-mass peak is sigma = 0.66
    # Msun wide against bins several Msun across, so the midpoint rule is wrong
    # there by tens of percent and would report a perfectly good sampler as
    # failing by 36 sigma. Differencing the cumulative integral is exact to the
    # resolution of ``grid``.
    cumulative = cumulative_trapezoid(density, grid, initial=0.0)
    expected = samples.size * (
        numpy.interp(edges[1:], grid, cumulative)
        - numpy.interp(edges[:-1], grid, cumulative)
    )

    usable = expected >= MINIMUM_EXPECTED_COUNTS
    if usable.sum() < 5:
        logging.warning(
            "%-28s only %d usable bins; too few events to test", label, usable.sum()
        )
        return 0.0, True
    residual = (counts[usable] - expected[usable]) / numpy.sqrt(expected[usable])
    worst = numpy.max(numpy.abs(residual))
    reduced = numpy.sum(residual**2) / residual.size
    # Look elsewhere: the worst of N bins is being tested, not one chosen in
    # advance, so the per-bin cut is the global rate divided by the number of
    # bins. Measured against 40 independent draws of a known-correct sampler
    # this sits just above the largest residual seen, where a flat 4 sigma
    # fired once in 160.
    threshold = stats.norm.isf(0.5 * FALSE_ALARM_PROBABILITY / residual.size)
    passed = worst < threshold and reduced < 2.0
    logging.info(
        "%-28s worst residual %5.2f sigma  (of %4.2f)   chi2/dof %5.2f   "
        "(%d bins)  %s",
        label,
        worst,
        threshold,
        reduced,
        residual.size,
        "OK" if passed else "FAIL",
    )
    return worst, passed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", required=True, help="gwforge_population HDF5")
    parser.add_argument("--config", required=True, help="the ini it was generated from")
    parser.add_argument("--output-directory", default=".")
    parser.add_argument("--seed", type=int, default=250114)
    options = parser.parse_args(argv)

    import h5py
    import os

    settings = read_config(options.config)
    with h5py.File(options.population, "r") as handle:
        events = {
            key: handle[key][:]
            for key in (
                "mass_1_source",
                "mass_2_source",
                "redshift",
                "a_1",
                "a_2",
                "tilt_1",
                "tilt_2",
                "ra",
                "dec",
                "theta_jn",
                "psi",
            )
        }
    events["cos_tilt_1"] = numpy.cos(events["tilt_1"])
    events["cos_tilt_2"] = numpy.cos(events["tilt_2"])
    events["mass_ratio"] = events["mass_2_source"] / events["mass_1_source"]
    count = len(events["redshift"])
    logging.info("Loaded %d injections from %s", count, options.population)

    mass = settings["mass"]
    spin = settings["spin"]
    model = mass_model(mass)
    os.makedirs(options.output_directory, exist_ok=True)
    results = []

    # ---- one-dimensional marginals ------------------------------------
    upper = min(mass.get("maximum_mass", 200), 120.0)
    primary_grid = numpy.linspace(mass["mmin"], upper, 4000)
    results.append(
        compare(
            events["mass_1_source"],
            primary_density(model, mass, primary_grid),
            primary_grid,
            "p(m1)",
        )
    )
    # p(m_2) is not a model input -- it has to be marginalised out of
    # pi(m_1) p(q | m_1) -- so checking it tests the pair jointly, in a way
    # neither the primary marginal nor the conditional does on its own.
    secondary_grid = numpy.linspace(
        mass.get("mmin_2", mass["mmin"]), upper, 1000
    )
    marginal_grid = numpy.linspace(mass["mmin"], upper, SECONDARY_MARGINAL_NODES)
    secondary = secondary_density(model, mass, marginal_grid, secondary_grid)
    results.append(
        compare(events["mass_2_source"], secondary, secondary_grid, "p(m2)")
    )
    redshift_grid = numpy.linspace(1e-3, settings["maximum_redshift"], 2000)
    results.append(
        compare(
            events["redshift"],
            redshift_density(settings, redshift_grid),
            redshift_grid,
            "p(z)",
        )
    )
    magnitude_grid = numpy.linspace(0.0, 1.0, 2000)
    results.append(
        compare(
            numpy.concatenate([events["a_1"], events["a_2"]]),
            magnitude_density(spin, magnitude_grid),
            magnitude_grid,
            "p(chi)",
        )
    )
    cosine_grid = numpy.linspace(-1.0, 1.0, 2000)
    results.append(
        compare(
            numpy.concatenate([events["cos_tilt_1"], events["cos_tilt_2"]]),
            cosine_tilt_density(spin, cosine_grid),
            cosine_grid,
            "p(cos t) marginal",
        )
    )

    # ---- the extrinsic priors ------------------------------------------
    extrinsic_grids = {}
    for name, (low, high) in EXTRINSIC_SPANS.items():
        grid = numpy.linspace(low, high, 2000)
        density = extrinsic_density(
            name, grid, settings["inclination_distribution"]
        )
        extrinsic_grids[name] = (grid, density)
        results.append(compare(events[name], density, grid, "p({})".format(name)))

    # ---- the conditional, which a marginal check cannot see ------------
    ratio_grid = numpy.linspace(1e-3, 1.0, 1000)
    for low, high in PRIMARY_MASS_BINS:
        inside = (events["mass_1_source"] > low) & (events["mass_1_source"] < high)
        if inside.sum() < 500:
            logging.warning(
                "only %d events in m1 (%g, %g); skipping", inside.sum(), low, high
            )
            continue
        # p(q | m_1) varies across the bin, so the comparison is against the
        # conditional *averaged over the primary masses actually in it* -- not
        # against the conditional at the bin's median, which would show a bias
        # that is an artefact of binning a conditional.
        chosen = events["mass_1_source"][inside]
        if chosen.size > 200:
            chosen = numpy.random.default_rng(0).choice(chosen, 200, replace=False)
        # One vectorised call: p_q rebuilds its per-m1 normalisation grid on every
        # invocation, so evaluating it once over the whole (m1, q) mesh is orders
        # of magnitude faster than looping.
        primary_mesh, ratio_mesh = numpy.meshgrid(chosen, ratio_grid, indexing="ij")
        averaged = ratio_density_mesh(
            model, mass, primary_mesh, ratio_mesh
        ).mean(axis=0)
        results.append(
            compare(
                events["mass_ratio"][inside],
                averaged,
                ratio_grid,
                "p(q | m1 in {:.0f}-{:.0f})".format(low, high),
            )
        )

    # ---- the joint, which a marginal check cannot see ------------------
    predicted = tilt_correlation(spin)
    measured = numpy.corrcoef(events["cos_tilt_1"], events["cos_tilt_2"])[0, 1]
    uncertainty = 1.0 / numpy.sqrt(count)
    deviation = abs(measured - predicted) / uncertainty
    logging.info(
        "corr(cos t1, cos t2)         measured %.4f +- %.4f;  joint model %.4f, "
        "factorised model 0  ->  %.1f sigma from joint  %s",
        measured,
        uncertainty,
        predicted,
        deviation,
        "OK" if deviation < 4 else "FAIL",
    )
    # The two models are only distinguishable here if the joint prediction is
    # itself well above the sampling noise. It scales as the squared mean of the
    # truncated Gaussian, so for a mu_t far from +-1 it is tiny and *this*
    # catalogue cannot tell the two apart however many events it has. Say so
    # rather than let a vacuous pass look like a verified one.
    separation = abs(predicted) / uncertainty
    if separation < 3.0:
        logging.warning(
            "  the joint and factorised tilt models differ by only %.1f sigma at "
            "these parameters, so this catalogue cannot distinguish them. The "
            "sampler's joint structure is verified instead by "
            "tests/test_population_spin.py, at mu_t = 1 where the correlation "
            "is 0.44.",
            separation,
        )
    results.append((deviation, deviation < 4))

    # ---- the figure ----------------------------------------------------
    draws = numpy.column_stack([events[name] for name in CORNER_PARAMETERS])
    span = [
        (mass["mmin"], min(upper, numpy.percentile(events["mass_1_source"], 99.9))),
        (mass.get("mmin_2", mass["mmin"]), numpy.percentile(events["mass_2_source"], 99.9)),
        (0.0, settings["maximum_redshift"]),
        (0.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 1.0),
    ] + [EXTRINSIC_SPANS[name] for name in ("ra", "dec", "theta_jn", "psi")]
    figure = corner_plot(
        draws,
        labels=labels_for(CORNER_PARAMETERS),
        span=span,
        levels=(0.5, 0.9),
        label="generated ({} injections)".format(count),
        colour=palette()[0][0],
    )
    # Overlay the analytic marginal on each diagonal: the model the sampler was
    # supposed to draw from, on top of what it drew.
    panels = [axis for axis in figure.axes if getattr(axis, "_gwforge_panel", False)]
    grid = numpy.array(panels).reshape(len(CORNER_PARAMETERS), len(CORNER_PARAMETERS))
    magnitude = (magnitude_grid, magnitude_density(spin, magnitude_grid))
    cosine = (cosine_grid, cosine_tilt_density(spin, cosine_grid))
    analytic = {
        "mass_1_source": (primary_grid, primary_density(model, mass, primary_grid)),
        "mass_2_source": (secondary_grid, secondary),
        "redshift": (redshift_grid, redshift_density(settings, redshift_grid)),
        # The two components are identically distributed, so drawing the same
        # curve on both diagonals is itself the assertion.
        "a_1": magnitude,
        "a_2": magnitude,
        "cos_tilt_1": cosine,
        "cos_tilt_2": cosine,
    }
    analytic.update(extrinsic_grids)
    with matplotlib.rc_context(new_rcParams(width="page")):
        for index, name in enumerate(CORNER_PARAMETERS):
            if name not in analytic:
                continue
            axis = grid[index, index]
            x, y = analytic[name]
            axis.plot(x, y, color="k", linewidth=1.2, zorder=6)
    figure.legends[-1].remove()
    figure.legend(
        handles=[
            matplotlib.lines.Line2D(
                [], [], color=palette()[0][0], linewidth=2.0,
                label="generated ({} injections)".format(count)),
            matplotlib.lines.Line2D([], [], color="k", linewidth=1.2, label="model"),
        ],
        loc="upper right",
        frameon=False,
    )
    path = os.path.join(options.output_directory, "gwtc5_population_check.pdf")
    figure.savefig(path, bbox_inches="tight")
    pylab.close(figure)
    logging.info("Wrote %s", path)

    passed = sum(1 for _, ok in results if ok)
    logging.info(
        "%d/%d checks passed; worst standardised residual %.2f sigma",
        passed,
        len(results),
        max(value for value, _ in results[:-1]),
    )
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
