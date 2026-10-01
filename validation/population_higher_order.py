#!/usr/bin/env python
r"""How far the Term-I population Fisher can be trusted for a given catalogue.

Run after ``gwforge_population_fisher`` and a per-event ``gwforge_fisher`` pass:

    python validation/population_higher_order.py \\
        --config-file spectral_siren_aplus.ini \\
        --result spectral_siren_aplus_joint_snr10.npz \\
        --event-fisher 'chunks/fisher_*.pkl' \\
        --output-directory higher_order

What this measures, and what it does not
----------------------------------------

Gair, Antonelli & Barbieri (`arXiv:2205.07893 <https://arxiv.org/abs/2205.07893>`_)
expand the hyper-parameter Fisher into five terms. GWForge computes the first.
Term I is the leading order of an expansion in the size of the *per-event
measurement uncertainty* relative to the curvature of the population density, so
the question "is Term I enough?" has a dimensionless answer, and it is cheaper
and more robust to compute than the corrections themselves.

For each detected event, with :math:`C_k` the measurement covariance of
:math:`\theta = (m_1^{\rm det}, m_2^{\rm det}, d_L)` and

.. math::

   H_k = \left.\frac{\partial^2 \ln p(\theta \mid \Lambda)}
                    {\partial\theta\,\partial\theta}\right|_{\theta_k},

the expansion parameter is the spectral radius :math:`\rho(H_k C_k)`. Terms
II--V are all built from :math:`A_k = F_k - H_k` with :math:`F_k = C_k^{-1}`,
and they are small exactly when :math:`\rho(H_k C_k) \ll 1`. The per-direction
version :math:`|H_{aa}| C_{aa}` says *which* coordinate is responsible.

Two mistakes are easy here and both are avoided:

* :math:`C_k` must be the **marginal** covariance -- a slice of
  :math:`\Sigma = \Gamma^{-1}` -- not a slice of the single-event Fisher, which
  conditions on the nuisance directions and is smaller by orders of magnitude.
* the spectral radius, not :math:`\max_a |H_{aa}| C_{aa}`. The correction is a
  matrix expansion, and correlations between the mass and distance directions
  make the true parameter larger than any diagonal entry.

The corrections themselves are deliberately **not** computed. The seminumeric
reference implements :math:`\Gamma_{\rm II}`--:math:`\Gamma_{\rm V}` by nested
finite differences and then declines to use them for a spectral-siren forecast,
because for a mass peak this sharp the correction sits below the
finite-difference noise floor and the raw sum returns non-positive-definite
matrices. A diagnostic that says "trust it" or "do not" is the honest product of
that machinery.
"""

import argparse
import glob
import logging
import os
import pickle

import matplotlib

matplotlib.use("Agg")

import numpy
import pylab

from GWForge.plotting import LVK_COLOUR, XG_COLOUR, new_rcParams

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# The observables the population density is written over.
OBSERVABLES = ("mass_1", "mass_2", "luminosity_distance")

# Single-event Fisher basis the marginal covariance is taken from, in the order
# that maps onto :data:`OBSERVABLES` through :func:`jacobian`.
EVENT_BASIS = ("chirp_mass", "symmetric_mass_ratio", "luminosity_distance")

# Fractional step for the theta-Hessian. The reference uses the same value;
# log_prob is smooth in theta, so this is not delicate.
THETA_STEP = 1e-3

# Half-width of the window around a mass peak, in units of its own sigma.
PEAK_WINDOW = 2.0


def load_event_fishers(pattern):
    """Per-event marginal covariances in :data:`EVENT_BASIS`.

    Parameters
    ----------
    pattern : str
        Glob matching the pickles ``gwforge_fisher`` wrote.

    Returns
    -------
    tuple
        ``(index, covariance, snr, n_failed)``. ``covariance`` has shape
        ``(n, 3, 3)`` and ``index`` numbers rows of the file the Fisher pass was
        given.
    """
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise ValueError("No per-event Fisher pickles match {!r}".format(pattern))

    index, covariance, snr, n_failed = [], [], [], 0
    for path in paths:
        with open(path, "rb") as handle:
            records = pickle.load(handle)
        for record in records:
            if "error" in record or "covariance" not in record:
                n_failed += 1
                continue
            names = list(record["names"])
            try:
                columns = [names.index(name) for name in EVENT_BASIS]
            except ValueError:
                raise ValueError(
                    "Event {} was forecast over {}, which does not contain {}. "
                    "The population basis cannot be recovered from it.".format(
                        record["index"], names, list(EVENT_BASIS)
                    )
                )
            # Slicing the *covariance* marginalises; slicing the Fisher would
            # condition on the other nine parameters instead.
            block = numpy.asarray(record["covariance"])[numpy.ix_(columns, columns)]
            if not numpy.all(numpy.isfinite(block)):
                n_failed += 1
                continue
            index.append(int(record["index"]))
            covariance.append(block)
            snr.append(float(record["optimal_snrs"]["network"]))
    order = numpy.argsort(index)
    return (
        numpy.asarray(index)[order],
        numpy.asarray(covariance)[order],
        numpy.asarray(snr)[order],
        n_failed,
    )


def jacobian(chirp_mass, symmetric_mass_ratio, step=1e-6):
    """``d(m1, m2, dL) / d(chirp_mass, symmetric_mass_ratio, dL)`` per event.

    Central differences rather than the closed form: the map runs through
    ``sqrt(1 - 4 eta)``, which is easy to differentiate wrongly and impossible
    to get wrong numerically.

    Parameters
    ----------
    chirp_mass, symmetric_mass_ratio : numpy.ndarray
    step : float
        Fractional displacement.

    Returns
    -------
    numpy.ndarray
        ``(n, 3, 3)``.
    """

    def components(mc, eta):
        total = mc / eta ** 0.6
        root = numpy.sqrt(numpy.clip(1.0 - 4.0 * eta, 0.0, None))
        return 0.5 * total * (1.0 + root), 0.5 * total * (1.0 - root)

    n = len(chirp_mass)
    matrix = numpy.zeros((n, 3, 3))
    matrix[:, 2, 2] = 1.0
    for column, values in enumerate((chirp_mass, symmetric_mass_ratio)):
        delta = step * numpy.abs(values)
        up = [chirp_mass, symmetric_mass_ratio]
        down = [chirp_mass, symmetric_mass_ratio]
        up[column] = values + delta
        down[column] = values - delta
        m1_up, m2_up = components(*up)
        m1_down, m2_down = components(*down)
        matrix[:, 0, column] = (m1_up - m1_down) / (2.0 * delta)
        matrix[:, 1, column] = (m2_up - m2_down) / (2.0 * delta)
    return matrix


def rotate(covariance, chirp_mass, symmetric_mass_ratio):
    """Covariance carried from :data:`EVENT_BASIS` into :data:`OBSERVABLES`.

    A covariance transforms with the Jacobian of the *forward* map,
    ``Sigma_y = J Sigma_x J^T``. Using the inverse here is the classic error and
    shrinks the mass errors by a large factor without warning.
    """
    matrix = jacobian(chirp_mass, symmetric_mass_ratio)
    return numpy.einsum("nai,nij,nbj->nab", matrix, covariance, matrix)


def theta_hessian(model, events, parameters, step=THETA_STEP):
    """Hessian of ``log_prob`` in the observables, ``(n_events, 3, 3)``.

    Central differences, vectorised over events -- ``log_prob`` returns one
    value per event, so the whole catalogue is differenced at once and the cost
    is 19 density evaluations regardless of how many events there are.
    """
    keys = list(OBSERVABLES)
    base = {key: numpy.asarray(events[key], dtype=float) for key in events}
    scale = numpy.array([numpy.abs(base[key]) * step for key in keys])
    scale = numpy.where(scale > 0.0, scale, step)

    def shifted(offsets):
        moved = dict(base)
        for position, key in enumerate(keys):
            if offsets[position]:
                moved[key] = base[key] + offsets[position] * scale[position]
        return model.log_prob(moved, parameters)

    n = len(base[keys[0]])
    centre = shifted([0, 0, 0])
    hessian = numpy.zeros((n, 3, 3))
    for i in range(3):
        plus = shifted([1 if k == i else 0 for k in range(3)])
        minus = shifted([-1 if k == i else 0 for k in range(3)])
        hessian[:, i, i] = (plus + minus - 2.0 * centre) / scale[i] ** 2
    for i in range(3):
        for j in range(i + 1, 3):
            def offset(si, sj):
                out = [0, 0, 0]
                out[i], out[j] = si, sj
                return shifted(out)

            mixed = (
                offset(1, 1) - offset(1, -1) - offset(-1, 1) + offset(-1, -1)
            ) / (4.0 * scale[i] * scale[j])
            hessian[:, i, j] = mixed
            hessian[:, j, i] = mixed
    return hessian


def diagnostic(hessian, covariance):
    """Per-direction and matrix expansion parameters.

    Returns
    -------
    tuple
        ``(per_direction, spectral_radius)`` of shapes ``(n, 3)`` and ``(n,)``.
    """
    per_direction = numpy.abs(
        numpy.einsum("naa->na", hessian) * numpy.einsum("naa->na", covariance)
    )
    product = numpy.einsum("nij,njk->nik", hessian, covariance)
    finite = numpy.all(numpy.isfinite(product), axis=(1, 2))
    radius = numpy.full(len(product), numpy.nan)
    radius[finite] = numpy.max(
        numpy.abs(numpy.linalg.eigvals(product[finite])), axis=1
    )
    return per_direction, radius


def quantiles(values, label):
    """One formatted line of summary statistics.

    Non-finite entries are excluded and counted rather than silently dropped:
    they are events a one-part-in-a-thousand displacement pushes outside the
    model's support, and how many there are is worth knowing.
    """
    finite = numpy.isfinite(values)
    good = values[finite]
    if not len(good):
        return "  {:>22s}   (nothing finite)".format(label)
    excluded = int((~finite).sum())
    return (
        "  {:>22s}   median {:9.3e}   90th {:9.3e}   max {:9.3e}   "
        "above 1: {:5.1f}%{}".format(
            label,
            numpy.median(good),
            numpy.percentile(good, 90),
            good.max(),
            100.0 * (good > 1.0).mean(),
            "   [{} non-finite]".format(excluded) if excluded else "",
        )
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-file", required=True, help="The Fisher run's ini")
    parser.add_argument(
        "--result", required=True, help=".npz the Fisher run saved, for the fiducial"
    )
    parser.add_argument(
        "--event-fisher", required=True, help="Glob for the gwforge_fisher pickles"
    )
    parser.add_argument(
        "--detected",
        required=True,
        help="HDF5 the per-event Fisher pass was run on, in the same row order",
    )
    parser.add_argument("--label", default="forecast", help="Name used in the figures")
    parser.add_argument("--output-directory", default="higher_order")
    options = parser.parse_args(argv)

    import configparser

    import h5py

    from GWForge.population_fisher import PopulationFisherResult
    from GWForge.population_fisher.config import build_model

    config = configparser.ConfigParser()
    if not config.read(options.config_file):
        raise FileNotFoundError(options.config_file)
    model, _, _, _ = build_model(config)
    result = PopulationFisherResult.from_npz(options.result)

    # The model carries the config values; the forecast was evaluated at the
    # fitted ones, and the Hessian has to be taken at the same point.
    parameters = dict(model.fiducial)
    parameters.update(result.fiducial)
    logging.info(
        "Hessian at the fitted fiducial: %s",
        ", ".join(
            "{}={:.4g}".format(name, result.fiducial[name])
            for name in ("H0", "mpp_1", "kappa")
            if name in result.fiducial
        ),
    )

    index, covariance, snr, n_failed = load_event_fishers(options.event_fisher)
    logging.info(
        "%d per-event covariances (%d unusable)", len(index), n_failed
    )

    with h5py.File(options.detected, "r") as handle:
        rows = {key: handle[key][:] for key in handle.keys()}
    events = {key: rows[key][index] for key in OBSERVABLES}
    chirp_mass = rows["chirp_mass"][index]
    eta = rows["symmetric_mass_ratio"][index]

    covariance = rotate(covariance, chirp_mass, eta)
    hessian = theta_hessian(model, events, parameters)
    per_direction, radius = diagnostic(hessian, covariance)

    os.makedirs(options.output_directory, exist_ok=True)
    rule = "-" * 78
    lines = [
        rule,
        "Term-I validity for {}".format(options.label),
        rule,
        "events with a usable covariance : {}".format(len(index)),
        "events dropped                  : {}".format(n_failed),
        "network SNR                     : median {:.1f}, max {:.1f}".format(
            numpy.median(snr), snr.max()
        ),
        "",
        "Expansion parameter -- Term I is the leading order in this, so it is",
        "trustworthy where the number is small and not where it approaches 1.",
        "",
    ]
    for position, key in enumerate(OBSERVABLES):
        lines.append(quantiles(per_direction[:, position], "|H| C  " + key))
    lines.append(quantiles(radius, "rho(H C)  [the one]"))

    # The global median hides the events that carry the information: the narrow
    # low-mass peak is a standard scale, and it is what H0 is measured against.
    peak, width = parameters.get("mpp_1"), parameters.get("sigpp_1")
    if peak is not None and width is not None:
        source = rows["mass_1_source"][index]
        window = numpy.abs(source - peak) <= PEAK_WINDOW * width
        lines += [
            "",
            "Restricted to |m1_source - {:.2f}| <= {:g} sigma ({} events, {:.1f}% "
            "of the sample):".format(peak, PEAK_WINDOW, int(window.sum()),
                                     100.0 * window.mean()),
            "",
        ]
        if window.any():
            for position, key in enumerate(OBSERVABLES):
                lines.append(
                    quantiles(per_direction[window, position], "|H| C  " + key)
                )
            lines.append(quantiles(radius[window], "rho(H C)  [the one]"))
    lines.append(rule)

    report = "\n".join(lines)
    print(report)
    summary = os.path.join(
        options.output_directory, "{}_higher_order.txt".format(options.label)
    )
    with open(summary, "w") as handle:
        handle.write(report + "\n")
    logging.info("Wrote %s", summary)

    with matplotlib.rc_context(new_rcParams(width="page")):
        figure, axes = pylab.subplots(1, 2, figsize=(11.0, 4.2))
        good = numpy.isfinite(radius)
        bins = numpy.geomspace(
            max(numpy.nanpercentile(radius[good], 0.5), 1e-6),
            max(numpy.nanmax(radius[good]), 2.0),
            40,
        )
        axes[0].hist(radius[good], bins=bins, color=XG_COLOUR, alpha=0.5)
        if peak is not None and window.any():
            axes[0].hist(
                radius[good & window], bins=bins, color=LVK_COLOUR, alpha=0.5,
                label=r"$m_1^{\rm src}$ within $2\sigma$ of the low-mass peak",
            )
            axes[0].legend(frameon=False, loc="upper left", fontsize="small")
        axes[0].axvline(1.0, color="0.2", linestyle="--", linewidth=1.2)
        axes[0].set_xscale("log")
        axes[0].set_xlabel(r"$\rho(H_k C_k)$")
        axes[0].set_ylabel("events")

        axes[1].scatter(snr[good], radius[good], s=5, color=XG_COLOUR, alpha=0.4)
        axes[1].axhline(1.0, color="0.2", linestyle="--", linewidth=1.2)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"network $\rho$")
        axes[1].set_ylabel(r"$\rho(H_k C_k)$")
        figure.tight_layout()
        path = os.path.join(
            options.output_directory, "{}_higher_order.pdf".format(options.label)
        )
        figure.savefig(path)
        pylab.close(figure)
    logging.info("Wrote %s", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
