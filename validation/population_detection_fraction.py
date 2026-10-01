#!/usr/bin/env python
"""Detection fraction of a generated catalogue, overall and where it bites.

Run between ``gwforge_optimal_snr`` and ``gwforge_population_fisher``:

    python validation/population_detection_fraction.py \\
        --population bbh_population.h5 --snr bbh_snr.h5 \\
        --detectors CE40 CE20 ET --snr-thresholds 10 20 \\
        --output-directory detection_fraction

Where P_det sits in the forecast
--------------------------------

It is easy to read the Term-I estimator as having no selection function, because
nothing in it evaluates one. That is not right. The first term of Gair+2022
Eq. 21 is an expectation *under the detected density*,

.. math::

   \\Gamma_{1}^{ab} = -\\int \\frac{\\partial^2}{\\partial\\Lambda^a \\partial\\Lambda^b}
     \\ln\\frac{p(\\theta|\\Lambda)}{\\beta(\\Lambda)}\\; p_{\\rm det}(\\theta|\\Lambda)\\,
     d\\theta,
   \\qquad
   p_{\\rm det} = \\frac{P_{\\rm det}(\\theta)\\, p(\\theta|\\Lambda)}{\\beta(\\Lambda)},

so :math:`P_{\\rm det}` is in there twice: as the measure, and through the
detectable fraction :math:`\\beta`. It collapses out of the *estimator* for two
reasons, not one. :math:`P_{\\rm det}(\\theta)` carries no hyper-parameter, so it
drops from the score; and :math:`\\partial_a \\ln\\beta = E_{\\rm det}[s^a]`, which
is the mean subtraction. What remains is carried entirely by **which events are
in the sum** -- the detected catalogue is the Monte-Carlo draw from
:math:`p_{\\rm det}`.

So the SNR threshold is a first-class input, not a preprocessing detail: change
it and :math:`\\Gamma_1` changes. This script reports what that threshold does.
The binned fractions are the interesting part -- they show which part of the
population is being removed, which is what the centring has to undo.

:math:`P_{\\rm det}` returns explicitly at higher order (it appears as
:math:`\\partial_\\theta \\ln P_{\\rm det}` in Gair's :math:`\\Gamma_4`), but
GWForge implements Term I only.
"""
import argparse
import logging
import os

import matplotlib

matplotlib.use("Agg")

import numpy
import pylab

from GWForge.plotting import LVK_COLOUR, XG_COLOUR, new_rcParams
from GWForge.population_fisher.catalogue import (
    expand_detectors,
    network_snr_from,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Bins for the two conditional fractions.
REDSHIFT_BINS = numpy.linspace(0.0, 10.0, 21)
PRIMARY_MASS_BINS = numpy.geomspace(5.0, 120.0, 21)

# A bin with fewer than this many injections is not plotted -- the binomial
# error on a handful of events is wider than anything it could show.
MINIMUM_PER_BIN = 20


def read(path, keys=None):
    """One-dimensional datasets from an HDF5 file."""
    import h5py

    data = {}
    with h5py.File(path, "r") as handle:
        for key in handle.keys():
            dataset = handle[key]
            if keys is not None and key not in keys:
                continue
            if isinstance(dataset, h5py.Dataset) and dataset.ndim == 1:
                data[key] = dataset[:]
    return data


def binned_fraction(values, detected, edges):
    """Detected fraction per bin, with its binomial error and the bin centres."""
    total, _ = numpy.histogram(values, bins=edges)
    found, _ = numpy.histogram(values[detected], bins=edges)
    usable = total >= MINIMUM_PER_BIN
    centre = 0.5 * (edges[1:] + edges[:-1])
    with numpy.errstate(invalid="ignore", divide="ignore"):
        fraction = numpy.where(usable, found / numpy.maximum(total, 1), numpy.nan)
        error = numpy.sqrt(fraction * (1.0 - fraction) / numpy.maximum(total, 1))
    return centre[usable], fraction[usable], error[usable], total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", required=True)
    parser.add_argument("--snr", required=True)
    parser.add_argument("--detectors", nargs="+", default=["CE40", "CE20", "ET"])
    parser.add_argument("--snr-thresholds", type=float, nargs="+", default=[10.0, 20.0])
    parser.add_argument("--output-directory", default="detection_fraction")
    options = parser.parse_args(argv)

    population = read(options.population)
    snr_data = read(options.snr)
    n_total = len(population["mass_1"])

    rows = numpy.asarray(snr_data["index"], dtype=int)
    detectors = expand_detectors(options.detectors)
    network_snr, used = network_snr_from(snr_data, detectors)
    logging.info("Network: %s", ", ".join(used))
    logging.info(
        "%d injections, %d with SNRs computed", n_total, len(network_snr)
    )
    logging.info(
        "Network SNR: median %.1f, 5th pct %.1f, max %.1f",
        numpy.median(network_snr),
        numpy.percentile(network_snr, 5),
        network_snr.max(),
    )

    redshift = population["redshift"][rows]
    primary = population["mass_1_source"][rows]

    os.makedirs(options.output_directory, exist_ok=True)
    lines = ["Detection fraction", "=" * 60, ""]
    lines.append("network        : {}".format(", ".join(used)))
    lines.append("injections     : {}".format(n_total))
    lines.append("SNRs computed  : {}".format(len(network_snr)))
    lines.append("")

    colours = [XG_COLOUR, LVK_COLOUR, "0.35", "0.6"]
    with matplotlib.rc_context(new_rcParams(width="page")):
        figure, axes = pylab.subplots(1, 2, figsize=(11.0, 4.2))

        for index, threshold in enumerate(options.snr_thresholds):
            detected = network_snr >= threshold
            fraction = detected.sum() / len(network_snr)
            lines.append(
                "SNR >= {:g}:  P_det = {:.4f}   N_det = {}".format(
                    threshold, fraction, int(detected.sum())
                )
            )
            logging.info(
                "SNR >= %g:  P_det = %.4f   N_det = %d",
                threshold,
                fraction,
                detected.sum(),
            )

            colour = colours[index % len(colours)]
            label = r"$\rho_{{\rm net}} \geq {:g}$".format(threshold)
            for axis, values, edges, scale in (
                (axes[0], redshift, REDSHIFT_BINS, "linear"),
                (axes[1], primary, PRIMARY_MASS_BINS, "log"),
            ):
                centre, value, error, _ = binned_fraction(values, detected, edges)
                axis.errorbar(
                    centre,
                    value,
                    yerr=error,
                    color=colour,
                    marker="o",
                    markersize=3.0,
                    linewidth=1.4,
                    capsize=2.0,
                    label=label,
                )
                axis.set_xscale(scale)

        axes[0].set_xlabel(r"$z$")
        axes[1].set_xlabel(r"$m_{1}^{\mathrm{source}}\ [M_{\odot}]$")
        # matplotlib's default log minor labels collide badly over this narrow a
        # decade range; name the few ticks worth reading instead.
        axes[1].set_xticks([5, 10, 20, 40, 80])
        axes[1].set_xticklabels(["5", "10", "20", "40", "80"])
        axes[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        for axis in axes:
            axis.set_ylabel(r"$P_{\mathrm{det}}$")
            axis.set_ylim(-0.03, 1.03)
            axis.legend(frameon=False, loc="lower left")
        figure.tight_layout()
        path = os.path.join(options.output_directory, "detection_fraction.pdf")
        figure.savefig(path)
        pylab.close(figure)
    logging.info("Wrote %s", path)

    lines.append("")
    lines.append(
        "P_det is not evaluated by the Term-I estimator; it is carried by which\n"
        "events enter the sum. See this script's module docstring."
    )
    summary = os.path.join(options.output_directory, "detection_fraction.txt")
    with open(summary, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    logging.info("Wrote %s", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
