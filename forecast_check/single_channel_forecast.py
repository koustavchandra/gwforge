#!/usr/bin/env python
"""Delay-time posterior for a single-channel (all isolated) population: Procedure 1 of the notes.

Takes a ``gwforge_population_fisher`` result, whose ``fiducial`` already holds the
maximum-likelihood Madau-Dickinson parameters fitted to all injected redshifts
(Eq. 15) and whose covariance is the inverse Term-I Fisher matrix scaled to the
detected count (Eq. A1), marginalises the mass and spin hyper-parameters by
taking the (gamma, kappa, z_peak) block, and evaluates the Eq. 18 grid posterior
on (alpha, log10 tau_min) with the star-formation rate pinned at Madau-Dickinson.
Everything but the wiring is ``forecast.py``.

Usage:
    single_channel_forecast.py --result pf_xg/mrs_xg_snr10.npz --label CE40+CE20+ET
Needs ../time_delay_forecasting_for_bbhs (forecast.py, time_delay_model.py) next to this directory.
"""
import argparse
import sys
from pathlib import Path

import numpy
import pylab

# The grid-posterior step and the delay model live in the (untracked) study directory.
STUDY = Path(__file__).resolve().parent.parent / "time_delay_forecasting_for_bbhs"
sys.path.insert(0, str(STUDY))
import forecast  # noqa: E402
from GWForge.plotting import TEALROSE, corner_plot, new_rcParams  # noqa: E402
from GWForge.population_fisher import PopulationFisherResult  # noqa: E402
from time_delay_model import ConvolutionGrid, PowerLawDelay  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
parser.add_argument("--result", required=True, type=Path, help="gwforge_population_fisher npz with the mass, redshift and spin blocks")
parser.add_argument("--label", required=True, help="network name, also the posterior file suffix")
parser.add_argument("--output-directory", type=Path, default=Path("delay"))
opts = parser.parse_args()

result = PopulationFisherResult.from_npz(str(opts.result))
mean, covariance = forecast.lambda_block(result)
print("{}: N_det = {}, lambda_hat = {}, sigma = {}".format(opts.label, result.n_events, mean.round(4), numpy.sqrt(numpy.diag(covariance)).round(4)))

forecast.OUTPUT = opts.output_directory
opts.output_directory.mkdir(parents=True, exist_ok=True)
grid = ConvolutionGrid(z_max=forecast.MAXIMUM_REDSHIFT, n_z=200, n_tau=2000, tau_floor=1e-3)  # as the reference
case = "single_channel_{}".format(opts.label.replace("+", "_"))
out = forecast.infer(case, mean, covariance, grid, PowerLawDelay(), 241, 281, opts.label)  # the notes' grid

samples = out["samples"]
print("  sigma(alpha) = {:.4f}, sigma(log10 tau_min) = {:.4f}, correlation = {:.3f}".format(samples[:, 0].std(), samples[:, 1].std(), numpy.corrcoef(samples.T)[0, 1]))
with pylab.rc_context(new_rcParams("column", aspect_ratio=1.0)):
    figure = corner_plot(samples, [r"$\alpha$", r"$\log_{10}\tau_{\min}/\mathrm{Gyr}$"], colour=TEALROSE[0], truths=[1.0, numpy.log10(0.02)], label=opts.label)
    figure.savefig(opts.output_directory / "time_delay_corner_{}.pdf".format(case))
print("  wrote", opts.output_directory / "posterior_{}.npz".format(case))
