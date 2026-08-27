#!/usr/bin/env python
"""
Forecast the same binary with two different waveform models and overlay them.

An aligned-spin binary is described equally well by ``IMRPhenomXHM`` and
``SEOBNRv5HM``, so both should recover the same signal-to-noise ratio and,
because the covariance is set by how the signal changes with the parameters,
very nearly the same error ellipses. Where they differ is informative: it is a
lower bound on the waveform systematic in any Fisher forecast, and it is not
something a single-model validation against another Fisher code can reveal.

Usage:
    python validation/compare_waveforms.py --outdir waveform_comparison
"""
import argparse
import logging
import os
import sys
import warnings

import numpy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import bilby  # noqa: E402
from GWForge.fisher import FisherMatrix, covariance  # noqa: E402
from GWForge.fisher.parameters import DEFAULT_PARAMETERS  # noqa: E402
from GWForge.fisher.plot import overlay_corner, samples_from_covariance  # noqa: E402
from GWForge.ifo.detectors import Network  # noqa: E402

warnings.filterwarnings("ignore", message=".*gwsignal_binary_black_hole is deprecated.*")
warnings.filterwarnings("ignore", message=".*UNREVIEWED.*")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

GPS_TIME = 1187008882.4

# One aligned-spin binary, described identically to both models.
SOURCE = dict(
    chirp_mass=28.0,
    symmetric_mass_ratio=0.24,
    chi_1=0.3,
    chi_2=-0.2,
    luminosity_distance=1500.0,
    theta_jn=0.9,
    psi=0.9,
    phase=1.3,
    geocent_time=GPS_TIME,
    ra=1.375,
    dec=-1.2108,
)

# ``(label, approximant, source model name)``.
MODELS = (
    ("IMRPhenomXHM", "IMRPhenomXHM", "lal_binary_black_hole"),
    ("SEOBNRv5HM", "SEOBNRv5HM", "gwsignal_binary_black_hole"),
)


def build(approximant, source_model_name, detectors, duration, sampling_frequency,
          minimum_frequency):
    """Fisher matrix for one waveform model on a freshly built network."""
    from GWForge.fisher.config import SOURCE_MODELS

    source_model = SOURCE_MODELS[source_model_name]
    if source_model is None:
        raise SystemExit("{} is unavailable".format(source_model_name))

    interferometers = Network(ifos=list(detectors)).initialise_ifos()
    for interferometer in interferometers:
        interferometer.minimum_frequency = minimum_frequency
    interferometers.set_strain_data_from_zero_noise(
        sampling_frequency=sampling_frequency,
        duration=duration,
        start_time=GPS_TIME - duration + 2,
    )
    waveform_generator = bilby.gw.WaveformGenerator(
        duration=duration,
        sampling_frequency=sampling_frequency,
        frequency_domain_source_model=source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=dict(
            waveform_approximant=approximant,
            reference_frequency=minimum_frequency,
            minimum_frequency=minimum_frequency,
        ),
    )
    return FisherMatrix(
        interferometers,
        waveform_generator,
        SOURCE,
        fisher_parameters=list(DEFAULT_PARAMETERS),
        order=1,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="waveform_comparison")
    parser.add_argument("--detectors", nargs="+", default=["CE40", "CE20", "ET"])
    parser.add_argument("--duration", type=float, default=8.0)
    parser.add_argument("--sampling-frequency", type=float, default=2048.0)
    parser.add_argument("--minimum-frequency", type=float, default=20.0)
    parser.add_argument(
        "--parameters",
        nargs="+",
        default=["chirp_mass", "symmetric_mass_ratio", "luminosity_distance",
                 "theta_jn", "chi_1", "chi_2"],
        help="Subset to plot; all eleven make an unreadably large figure",
    )
    opts = parser.parse_args()
    os.makedirs(opts.outdir, exist_ok=True)

    results = {}
    for label, approximant, source_model_name in MODELS:
        logging.info("Building the %s Fisher matrix", label)
        results[label] = build(
            approximant,
            source_model_name,
            opts.detectors,
            opts.duration,
            opts.sampling_frequency,
            opts.minimum_frequency,
        )

    first, second = (label for label, _, _ in MODELS)
    names = results[first].names

    logging.info("=" * 78)
    logging.info(
        "network SNR: %s %.2f, %s %.2f  (%.2f%%)",
        first,
        results[first].optimal_snrs["network"],
        second,
        results[second].optimal_snrs["network"],
        100
        * (
            results[second].optimal_snrs["network"]
            / results[first].optimal_snrs["network"]
            - 1
        ),
    )
    logging.info(
        "%-22s %12s %12s %8s %10s", "parameter", first, second, "ratio", "noise"
    )
    noise = results[second].derivatives.noise_ratios
    for name in names:
        one, two = results[first].sigmas[name], results[second].sigmas[name]
        logging.info(
            "%-22s %12.4e %12.4e %8.3f %9.0f%%",
            name,
            one,
            two,
            two / one,
            100 * noise.get(name, 0.0),
        )
    logging.info("=" * 78)
    logging.info(
        "The 'noise' column is %s's own numerical noise as a fraction of the "
        "change one finite-difference step produces -- the floor on how well "
        "that sigma can agree with anything.",
        second,
    )

    columns = [names.index(name) for name in opts.parameters]
    mean = numpy.array([SOURCE[name] for name in names])
    datasets = []
    for index, label in enumerate((first, second)):
        matrix, _ = covariance(results[label].matrix, names=names)
        datasets.append(
            (
                label,
                samples_from_covariance(mean, matrix, seed=42 + index)[:, columns],
            )
        )
    path = os.path.join(opts.outdir, "waveform_comparison.png")
    overlay_corner(
        datasets,
        opts.parameters,
        truths=[SOURCE[name] for name in opts.parameters],
        save=path,
        title="Same binary, two waveform models (network SNR {:.0f})".format(
            results[first].optimal_snrs["network"]
        ),
    )
    logging.info("Done!")


if __name__ == "__main__":
    main()
