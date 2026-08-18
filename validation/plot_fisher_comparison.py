#!/usr/bin/env python
"""
Overlay GWForge and GWFast forecasts on one corner plot per source.

Reads the HDF5 written by ``fisher_vs_gwfast.py``, turns each covariance matrix
into the multivariate normal it describes, and plots three things together:

* GWForge order 1 -- the Fisher covariance, drawn as a multivariate normal
* GWForge order 2 -- the doublet-DALI posterior, sampled with emcee
* GWFast          -- its Fisher covariance, drawn the same way as GWForge's

The first and third should sit on top of each other; the second shows where the
non-Gaussian correction moves the contours.

Usage:
    python validation/plot_fisher_comparison.py comparison.h5 --outdir plots
"""
import argparse
import logging
import os
import sys

import h5py
import numpy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from GWForge.fisher.dali import DALILikelihood, sample_posterior  # noqa: E402
from GWForge.fisher.plot import overlay_corner, samples_from_covariance  # noqa: E402

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

#: How many sigmas either side of the fiducial point to sample the DALI
#: posterior over. Wide enough not to truncate the posterior -- a prior narrower
#: than the likelihood silently sets the width -- but not so wide that the
#: expansion is sampled far outside where it is trustworthy.
PRIOR_WIDTH_IN_SIGMAS = 8.0


class _Tensors:
    """Minimal stand-in exposing what :class:`DALILikelihood` reads."""

    def __init__(self, names, parameters, fisher, doublet, quartic):
        self.names = names
        self.parameters = parameters
        self.matrix = fisher
        self.doublet = doublet
        self.quartic = quartic


def dali_samples(group, names, parameters, outdir, label):
    """Sample the stored doublet-DALI posterior, or None if it was not stored."""
    if "dali_doublet" not in group:
        return None
    from bilby.core.prior import PriorDict, Uniform

    tensors = _Tensors(
        list(names),
        parameters,
        group["gwforge_fisher"][:],
        group["dali_doublet"][:],
        group["dali_quartic"][:],
    )
    likelihood = DALILikelihood(tensors, order=2)
    sigmas = numpy.sqrt(numpy.diag(group["gwforge_covariance"][:]))
    priors = PriorDict(
        {
            name: Uniform(
                parameters[name] - PRIOR_WIDTH_IN_SIGMAS * sigma,
                parameters[name] + PRIOR_WIDTH_IN_SIGMAS * sigma,
                name=name,
            )
            for name, sigma in zip(names, sigmas)
        }
    )
    result = sample_posterior(
        likelihood, priors, outdir=outdir, label=label, nwalkers=100, nsteps=3000, nburn=1000
    )
    return result.posterior[list(names)].to_numpy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison", help="HDF5 written by fisher_vs_gwfast.py")
    parser.add_argument("--outdir", default="fisher_comparison", help="Where to write plots")
    parser.add_argument(
        "--parameters",
        nargs="+",
        default=None,
        help="Subset to plot; the full eleven make a very large figure",
    )
    parser.add_argument(
        "--no-dali", action="store_true", help="Skip the order-2 overlay and its sampling"
    )
    parser.add_argument("--seed", type=int, default=42)
    opts = parser.parse_args()

    os.makedirs(opts.outdir, exist_ok=True)
    with h5py.File(opts.comparison, "r") as handle:
        names = [
            name.decode() if isinstance(name, bytes) else str(name)
            for name in handle.attrs["parameters"]
        ]
        chosen = opts.parameters or names
        missing = [name for name in chosen if name not in names]
        if missing:
            raise SystemExit(
                "unknown parameter(s) {}; the file has {}".format(
                    ", ".join(missing), ", ".join(names)
                )
            )
        columns = [names.index(name) for name in chosen]

        for key in sorted(handle.keys(), key=lambda name: int(name.split("_")[1])):
            group = handle[key]
            parameters = {
                name: float(group["injection"].attrs[name])
                for name in group["injection"].attrs
            }
            mean = numpy.array([parameters[name] for name in names])

            datasets = [
                (
                    "GWForge Fisher",
                    samples_from_covariance(
                        mean, group["gwforge_covariance"][:], seed=opts.seed
                    )[:, columns],
                ),
                (
                    "GWFast Fisher",
                    samples_from_covariance(
                        mean, group["gwfast_covariance"][:], seed=opts.seed + 1
                    )[:, columns],
                ),
            ]
            if not opts.no_dali:
                samples = dali_samples(
                    group, names, parameters, opts.outdir, "dali_{}".format(key)
                )
                if samples is not None:
                    datasets.insert(1, ("GWForge DALI (order 2)", samples[:, columns]))

            path = os.path.join(opts.outdir, "{}.png".format(key))
            overlay_corner(
                datasets,
                chosen,
                truths=[parameters[name] for name in chosen],
                save=path,
                title="{}: network SNR {:.0f} (GWForge) vs {:.0f} (GWFast)".format(
                    key,
                    float(group["gwforge_snr"][()]),
                    float(group["gwfast_snr"][()]),
                ),
            )
    logging.info("Done!")


if __name__ == "__main__":
    main()
