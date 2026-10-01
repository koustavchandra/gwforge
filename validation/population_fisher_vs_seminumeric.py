#!/usr/bin/env python
"""Reproduce the seminumeric population-Fisher forecasts with GWForge.

The `population_fisher_seminumeric` tree is an independent implementation of the
same estimator -- Gair, Antonelli & Barbieri (2022) (`arXiv:2205.07893 <https://arxiv.org/abs/2205.07893>`_)
Eq. 21 Term I -- written before GWForge had one, and its published forecasts are
checked in. This runs GWForge's models on *that code's own catalogue*, so the
detected events are identical and any difference in sigma is a difference in the
code rather than Monte-Carlo scatter.

    python validation/population_fisher_vs_seminumeric.py \
        --reference-directory /path/to/population_fisher_seminumeric

Two conventions have to be reconciled first, and getting either wrong produces a
plausible-looking disagreement:

* **The Madau-Dickinson exponent.** GWForge writes
  ``psi ∝ (1+z)^gamma / [1 + ((1+z)/(1+z_p))^kappa]``; the reference's
  denominator exponent is ``alpha + beta``. So ``kappa = alpha + beta`` is a
  change of variables, not a rename, and the Fisher must be pushed through the
  Jacobian below before sigma can be compared.
* **The low-mass taper.** The reference applies one taper to *both* component
  masses. GWForge's BGP has independent primary and secondary tapers
  (GWTC-5.0 Tab. 5), so it is a different density unless ``mmin_2 = mmin`` and
  ``delta_m_2 = delta_m``. This script collapses them; a comparison without that
  collapse would be measuring physics, not agreement.

What is *not* a convention difference, checked rather than assumed: both codes
normalise the three mixture components individually, so ``lam_0`` really is
``lambda_0``; and the reference's Planck18 ``dV_c/dz`` and GWForge's extra
``[1 + (1+z_p)^-kappa]`` factor are both event-independent, so they shift score
columns by a constant that the centring removes and cannot move sigma.

The spin sector has no counterpart: the reference defines no spin model at all,
only a sampler with no ``log_prob``. It keeps the finite-difference oracle in
``tests/test_population_fisher.py``.

What this found
---------------

The population sector reproduces the published forecasts essentially exactly:
score columns agreeing to 1e-11 per event, and sigma to a few times 1e-6. The
redshift block agrees to every digit printed once the Jacobian above is applied.

The cosmology sector does not, quite -- sigma(H0), sigma(Om) and sigma(w0) come
out ~1.3% smaller. That is the reference's finite differences, not a GWForge
error; see :data:`TOLERANCE` for the per-event evidence.
"""
import argparse
import logging
import os
import sys

import numpy

from GWForge.population_fisher import (
    BrokenPowerLawTwoPeakMass,
    JointPopulationModel,
    MadauDickinsonRedshift,
    SpectralSirenModel,
    population_fisher,
)

logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# Where the reference lives, overridable from the environment so the test that
# wraps this can skip cleanly when it is not checked out.
DEFAULT_REFERENCE_DIRECTORY = os.environ.get(
    "GWFORGE_SEMINUMERIC_DIRECTORY",
    "/Users/kchandra/projects/CE_STM_CBC/population_fisher_seminumeric",
)

# Its catalogue: CE40 + CE20 + ET at 5 Hz, the file every checked-in forecast
# was produced from.
REFERENCE_CATALOGUE = (
    "network_bbh_CE40km_1p5MW_Aplus_coat_5.0hz_CE20km_1p5MW_Aplus_coat_5.0hz_ETD_5.0hz.h5"
)

# Reference mass name -> GWForge mass name. A pure rename in both directions.
MASS_NAMES = {
    "alpha_1": "alpha_1",
    "alpha_2": "alpha_2",
    "m_break": "m_break",
    "lambda_0": "lam_0",
    "lambda_1": "lam_1",
    "mu_1": "mpp_1",
    "sigma_1": "sigpp_1",
    "mu_2": "mpp_2",
    "sigma_2": "sigpp_2",
    "beta_q": "beta",
    "delta_m": "delta_m",
}

# Reference redshift names, in the order the Jacobian below assumes.
REDSHIFT_NAMES = ["alpha", "beta", "z_peak"]

# d(gamma, kappa, z_peak) / d(alpha, beta, z_peak) for kappa = alpha + beta.
# The Fisher is a rank-2 covariant tensor, so Gamma_ref = J^T Gamma_gwforge J.
REDSHIFT_JACOBIAN = numpy.array(
    [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
)

# Reference cosmology name -> GWForge name.
COSMOLOGY_NAMES = {"H0": "H0", "Om": "Om0", "w0": "w0"}

# Relative tolerance on sigma for the pure-population comparison. The two codes'
# score columns agree to 1e-11 per event there; the residual is the five
# significant figures of ``m_break`` and ``delta_m`` in the summary file.
POPULATION_TOLERANCE = 1e-4

# Relative tolerance for the spectral-siren comparison, where **every** marginal
# is looser -- not just the cosmology ones. The reference computes its ``Om``
# and ``w0`` score columns by embedded finite differences; GWForge computes them
# analytically. Per event the two agree to a median 8e-6 there, but ~0.05% of
# events differ by up to 5%, the sparse spiky signature of finite-difference
# noise rather than an analytic error. Marginalisation then spreads those few
# events across the whole matrix: H0 is 0.95 correlated with Om and 0.90 with
# mu_1, so sigma(mu_1) moves even though the mass score columns are identical to
# 1e-11. The script prints the measured column agreement under each table so
# this is checkable rather than asserted.
SPECTRAL_SIREN_TOLERANCE = 2e-2

# Where the mass+redshift forecast lives, by SNR threshold.
MASS_REDSHIFT_NPZ = "forecast_mass_redshift_output/forecast_joint_snr{:.0f}.npz"

# A full corner is 12 or 15 parameters and legible in none of them; these are
# the ones worth looking at. Mass shape, then the redshift block.
MASS_REDSHIFT_PLOT = ["alpha_1", "alpha_2", "beta_q", "alpha", "beta", "z_peak"]

# For spectral sirens, the cosmology plus the two population parameters it is
# most correlated with -- which is where the 1.3% comes from.
SPECTRAL_SIREN_PLOT = ["H0", "Om", "w0", "mu_1", "alpha", "z_peak"]


def reference_modules(directory):
    """Import the reference package from ``directory``, or explain why not."""
    if not os.path.isdir(directory):
        raise SystemExit(
            "No seminumeric reference at {}. Pass --reference-directory or set "
            "GWFORGE_SEMINUMERIC_DIRECTORY.".format(directory)
        )
    if directory not in sys.path:
        sys.path.insert(0, directory)
    import matplotlib

    matplotlib.use("Agg")
    import forecast_mass_redshift
    import load
    import spectral_sirens_forecast

    return load, forecast_mass_redshift, spectral_sirens_forecast


def reference_fiducial(directory, forecast_mass_redshift, spectral_sirens_forecast):
    """The reference's *published* fiducial, not a re-derivation of it.

    The published forecasts were produced at a particular maximum-likelihood
    optimum, and re-running that fit does not land on it again: the BGP
    likelihood is flat enough here that the optimiser stops 1e-5 to 6e-4 away,
    parameter by parameter, and sigma inherits exactly that. So the fiducial is
    reassembled from what the reference recorded:

    * the nine free parameters from each forecast's own ``.npz``, exactly;
    * ``m_min`` and ``m_max`` recomputed from the injected mass range, which is
      how the reference sets them -- not fitted, so this reproduces them to the
      last digit;
    * ``m_break`` and ``delta_m`` parsed from ``forecast_summary.txt``, the only
      place the run recorded them, at the five significant figures it printed.

    Comparing at a re-fitted fiducial would be comparing two different models
    and calling the difference a code discrepancy.
    """
    path = os.path.join(directory, REFERENCE_CATALOGUE)
    mass_1, mass_2 = forecast_mass_redshift.load_injected_masses(path)
    fiducial = {
        "m_min": 0.999 * float(min(mass_1.min(), mass_2.min())),
        "m_max": 1.001 * float(mass_1.max()),
    }
    fiducial.update(
        _parse_summary_fiducial(
            os.path.join(
                directory, "forecast_mass_redshift_output", "forecast_summary.txt"
            )
        )
    )
    fiducial.update(forecast_mass_redshift.MD_FIDUCIAL)
    fiducial.update(spectral_sirens_forecast.FIDUCIAL_COSMOLOGY)
    return fiducial


def _parse_summary_fiducial(path):
    """The ``name = value`` block the reference prints above its tables."""
    values = {}
    with open(path) as handle:
        for line in handle:
            if line.strip().startswith("====="):
                break
            if "=" in line and len(line.split("=")) == 2:
                name, value = (part.strip() for part in line.split("="))
                try:
                    values[name] = float(value)
                except ValueError:
                    continue
    return values


def reference_sigma(directory, relative_path):
    """``sigma``, and the exact fiducial, from one checked-in ``.npz`` forecast.

    Returns ``(sigma_by_name, fiducial_by_name)``; the second is the published
    optimum for the free parameters, which is what the comparison must be run
    at. ``None`` when that forecast is not checked in.
    """
    path = os.path.join(directory, relative_path)
    if not os.path.exists(path):
        return None, None
    archive = numpy.load(path, allow_pickle=True)
    names = [str(name) for name in archive["parameter_names"]]
    sigma = numpy.sqrt(numpy.diag(archive["covariance"]))
    fiducial = dict(zip(names, [float(value) for value in archive["fiducial"]]))
    return dict(zip(names, sigma)), fiducial


def reference_covariance(directory, relative_path):
    """The published covariance matrix, in the reference's own parameter order."""
    archive = numpy.load(os.path.join(directory, relative_path), allow_pickle=True)
    return numpy.asarray(archive["covariance"], dtype=float)


def gwforge_mass_model(fiducial):
    """A GWForge BGP model collapsed onto the reference's single-taper density.

    ``mmin_2``/``delta_m_2`` are set equal to their primary counterparts, which
    is the only configuration in which the two codes describe the same
    ``p(m_2 | m_1)``.
    """
    return BrokenPowerLawTwoPeakMass(
        mmin=fiducial["m_min"],
        m_high=fiducial["m_max"],
        maximum_mass=fiducial["m_max"],
        mmin_2=fiducial["m_min"],
        delta_m_2=fiducial["delta_m"],
        fiducial={
            gwforge: fiducial[reference]
            for reference, gwforge in MASS_NAMES.items()
        },
    )


def gwforge_redshift_model(fiducial, maximum_redshift, cosmology=None):
    """A GWForge Madau-Dickinson model at the reference's fiducial.

    ``kappa = alpha + beta``; see the module docstring.
    """
    return MadauDickinsonRedshift(
        maximum_redshift=maximum_redshift,
        cosmology=cosmology,
        fiducial={
            "gamma": fiducial["alpha"],
            "kappa": fiducial["alpha"] + fiducial["beta"],
            "z_peak": fiducial["z_peak"],
        },
    )


def to_reference_basis(result, names):
    """``sigma`` in the reference's parameterisation, for the names it uses.

    Mass and cosmology are renames, so their sigmas transfer untouched. The
    redshift block is a change of variables, so its Fisher sub-block is pushed
    through :data:`REDSHIFT_JACOBIAN` and re-inverted. Doing it on the sub-block
    is exact here only because the transformation is block-diagonal -- it mixes
    the three redshift parameters with each other and with nothing else -- so
    the full-matrix transform reduces to this.
    """
    gwforge_names = list(result.parameter_names)
    fisher = numpy.asarray(result.fisher, dtype=float)

    order = []
    for name in names:
        if name in MASS_NAMES:
            order.append(gwforge_names.index(MASS_NAMES[name]))
        elif name in COSMOLOGY_NAMES:
            order.append(gwforge_names.index(COSMOLOGY_NAMES[name]))
        elif name == "alpha":
            order.append(gwforge_names.index("gamma"))
        elif name == "beta":
            order.append(gwforge_names.index("kappa"))
        elif name == "z_peak":
            order.append(gwforge_names.index("z_peak"))
        else:
            raise KeyError("No GWForge counterpart for {!r}".format(name))

    permuted = fisher[numpy.ix_(order, order)]

    redshift_block = [index for index, name in enumerate(names) if name in REDSHIFT_NAMES]
    if len(redshift_block) == 3:
        jacobian = numpy.eye(len(names))
        jacobian[numpy.ix_(redshift_block, redshift_block)] = REDSHIFT_JACOBIAN
        permuted = jacobian.T @ permuted @ jacobian

    return numpy.linalg.inv(permuted)


def gwforge_covariance(result, names):
    """GWForge's own covariance, its columns ordered like ``names``.

    Untransformed -- this is already in GWForge's parameterisation; only the
    ordering is borrowed from the reference so the two matrices line up.
    """
    order = [
        list(result.parameter_names).index(gwforge_name(name)) for name in names
    ]
    return numpy.asarray(result.covariance, dtype=float)[numpy.ix_(order, order)]


def sigma_of(covariance, names):
    """``{name: sigma}`` from a covariance ordered like ``names``."""
    return dict(zip(names, numpy.sqrt(numpy.diag(covariance))))


def gwforge_name(name):
    """The GWForge spelling of one reference parameter name."""
    if name in MASS_NAMES:
        return MASS_NAMES[name]
    if name in COSMOLOGY_NAMES:
        return COSMOLOGY_NAMES[name]
    return {"alpha": "gamma", "beta": "kappa", "z_peak": "z_peak"}[name]


def overlay(reference_covariance, covariance, names, fiducial, title, path, show):
    """Overlay the two forecasts' contours, in **GWForge's** parameterisation.

    The tables above compare sigma in the reference's variables, because that is
    what its ``.npz`` files hold. The figure goes the other way and pushes the
    reference's covariance *forward* into GWForge's, for two reasons: the
    parameters then carry their proper labels, and ``alpha``/``beta`` mean
    different things in the two codes, so a figure labelled in reference
    variables invites exactly the confusion the Jacobian exists to avoid.

    Covariance is contravariant where the Fisher is covariant, so where
    ``Gamma_ref = J^T Gamma_gwforge J`` this is ``Sigma_gwforge = J Sigma_ref J^T``.
    """
    from GWForge.population_fisher.plot import corner_plot
    from GWForge.population_fisher.population_fisher_term_I import (
        PopulationFisherResult,
    )

    jacobian = numpy.eye(len(names))
    block = [index for index, name in enumerate(names) if name in REDSHIFT_NAMES]
    if len(block) == 3:
        jacobian[numpy.ix_(block, block)] = REDSHIFT_JACOBIAN
    pushed = jacobian @ reference_covariance @ jacobian.T

    gwforge_names = [gwforge_name(name) for name in names]
    values = {}
    for name, gwforge in zip(names, gwforge_names):
        values[gwforge] = fiducial[name]
    values["kappa"] = fiducial["alpha"] + fiducial["beta"]
    values["gamma"] = fiducial["alpha"]

    def wrap(matrix):
        return PopulationFisherResult(
            parameter_names=list(gwforge_names),
            fixed_parameters={},
            fiducial={name: values[name] for name in gwforge_names},
            fisher=numpy.linalg.inv(matrix),
            covariance=matrix,
            condition_number=float(numpy.linalg.cond(matrix)),
            n_events=0,
            n_total=0,
            scores=None,
            mean_score=None,
            degenerate=[],
        )

    corner_plot(
        [("seminumeric", wrap(pushed)), ("GWForge", wrap(covariance))],
        save=path,
        title=title,
        seed=250114,
        parameters=[gwforge_name(name) for name in show if name in names],
    )
    print("  wrote {}".format(path))


def report(title, reference, measured, tolerance):
    """Print one comparison table; return True if every parameter agrees."""
    print("\n" + "=" * 72)
    print("  " + title)
    print("=" * 72)
    print(
        "  {:<12s}{:>14s}{:>14s}{:>12s}   {}".format(
            "parameter", "seminumeric", "GWForge", "ratio", "")
    )
    print("  " + "-" * 66)
    passed = True
    for name, expected in reference.items():
        if name not in measured:
            continue
        got = measured[name]
        ratio = got / expected
        ok = abs(ratio - 1.0) <= tolerance
        passed = passed and ok
        print(
            "  {:<12s}{:>14.6g}{:>14.6g}{:>12.7f}   {}".format(
                name, expected, got, ratio, "OK" if ok else "FAIL"
            )
        )
    return passed


def compare_mass_and_redshift(directory, load, fiducial, snr_threshold, output_directory=None):
    """The `forecast_mass_redshift.py` forecasts: mass, redshift, and the joint."""
    reference, published = reference_sigma(
        directory,
        "forecast_mass_redshift_output/forecast_joint_snr{:.0f}.npz".format(
            snr_threshold
        ),
    )
    if reference is None:
        print("  no checked-in joint forecast at SNR >= {:.0f}".format(snr_threshold))
        return True
    fiducial = dict(fiducial, **published)
    catalogue = load.load_population_catalogue(
        os.path.join(directory, REFERENCE_CATALOGUE),
        snr_threshold=snr_threshold,
        with_fisher=False,
    )
    events = {
        "mass_1_source": catalogue.events["mass_1_source"],
        "mass_2_source": catalogue.events["mass_2_source"],
        "redshift": catalogue.events["redshift"],
    }
    print(
        "\nSNR >= {:.0f}: {} detected events".format(
            snr_threshold, len(events["redshift"])
        )
    )

    mass = gwforge_mass_model(fiducial)
    redshift = gwforge_redshift_model(
        fiducial, maximum_redshift=float(events["redshift"].max()) * 1.001
    )
    model = JointPopulationModel([mass, redshift])

    free = [MASS_NAMES[name] for name in reference if name in MASS_NAMES]
    free += ["gamma", "kappa", "z_peak"]
    result = population_fisher(model, events, free_parameters=free)

    names = list(reference)
    covariance = to_reference_basis(result, names)
    passed = report(
        "mass + redshift, SNR >= {:.0f}".format(snr_threshold),
        reference,
        sigma_of(covariance, names),
        POPULATION_TOLERANCE,
    )
    gwforge = gwforge_covariance(result, names)
    if output_directory is not None:
        overlay(
            reference_covariance(directory, MASS_REDSHIFT_NPZ.format(snr_threshold)),
            gwforge,
            names,
            fiducial,
            "mass + redshift, SNR >= {:.0f}".format(snr_threshold),
            os.path.join(
                output_directory,
                "mass_redshift_snr{:.0f}.pdf".format(snr_threshold),
            ),
            MASS_REDSHIFT_PLOT,
        )
    return passed


def compare_spectral_sirens(directory, load, fiducial, snr_threshold, output_directory=None):
    """The `spectral_sirens_forecast.py` joint (H0, Om, w0) forecast."""
    relative = "spectral_sirens_output/spectral_joint_snr{:.0f}.npz".format(
        snr_threshold
    )
    reference, published = reference_sigma(directory, relative)
    if reference is None:
        print("  no checked-in spectral forecast at SNR >= {:.0f}".format(snr_threshold))
        return True
    from GWForge.cosmology import FlatwCDM

    fiducial = dict(fiducial, **published)
    cosmology = FlatwCDM(
        H0=fiducial["H0"], Om0=fiducial["Om"], w0=fiducial["w0"]
    )

    catalogue = load.load_population_catalogue(
        os.path.join(directory, REFERENCE_CATALOGUE),
        snr_threshold=snr_threshold,
        with_fisher=False,
    )
    redshift = catalogue.events["redshift"]
    one_plus_z = 1.0 + redshift
    events = {
        "mass_1": catalogue.events["mass_1_source"] * one_plus_z,
        "mass_2": catalogue.events["mass_2_source"] * one_plus_z,
        "luminosity_distance": cosmology.luminosity_distance(redshift),
    }

    mass = gwforge_mass_model(fiducial)
    redshift_model = gwforge_redshift_model(
        fiducial,
        maximum_redshift=float(redshift.max()) * 1.001,
        cosmology=cosmology,
    )
    model = SpectralSirenModel(
        mass_model=mass, redshift_model=redshift_model, cosmology=cosmology
    )

    # Free exactly what the published run freed -- the .npz names it. Leaving
    # m_break and delta_m free would make the matrix singular: at this fiducial
    # the taper is 0.0039 Msun wide against a lightest secondary of 2.06, so no
    # event sits inside it and the delta_m column is identically zero.
    free = []
    for name in reference:
        if name in MASS_NAMES:
            free.append(MASS_NAMES[name])
        elif name in COSMOLOGY_NAMES:
            free.append(COSMOLOGY_NAMES[name])
        else:
            free.append({"alpha": "gamma", "beta": "kappa", "z_peak": "z_peak"}[name])
    result = population_fisher(model, events, free_parameters=free)
    names = list(reference)
    covariance = to_reference_basis(result, names)
    passed = report(
        "spectral sirens, SNR >= {:.0f}".format(snr_threshold),
        reference,
        sigma_of(covariance, names),
        SPECTRAL_SIREN_TOLERANCE,
    )
    gwforge = gwforge_covariance(result, names)
    if output_directory is not None:
        overlay(
            reference_covariance(directory, relative),
            gwforge,
            names,
            fiducial,
            "spectral sirens, SNR >= {:.0f}".format(snr_threshold),
            os.path.join(
                output_directory,
                "spectral_sirens_snr{:.0f}.pdf".format(snr_threshold),
            ),
            SPECTRAL_SIREN_PLOT,
        )
    score_agreement(
        directory, fiducial, catalogue, model, events, snr_threshold
    )
    return passed


def score_agreement(directory, fiducial, catalogue, model, events, snr_threshold):
    """Print how well the two codes' *score columns* agree, per event.

    The evidence for the loosened tolerance above. sigma is a marginal, so one
    noisy column contaminates all of them; the columns themselves are where the
    two implementations can be compared without that mixing.
    """
    import pop_models
    import spectral_sirens_forecast

    cosmology = {name: fiducial[name] for name in ("H0", "Om", "w0")}
    reference_mass = pop_models.BrokenPowerLawPlusTwoPeaksMass(
        fiducial={
            name: fiducial[name]
            for name in pop_models.BrokenPowerLawPlusTwoPeaksMass().parameter_names
        }
    )
    reference_redshift = pop_models.MadauDickinsonRedshift(
        fiducial={name: fiducial[name] for name in REDSHIFT_NAMES}
    )
    reference_model = spectral_sirens_forecast.SpectralSirenModel(
        reference_mass, reference_redshift, cosmology
    )
    reference_events = spectral_sirens_forecast.detector_frame_events(
        catalogue, cosmology
    )
    reference_score = reference_model.analytic_score(
        reference_events, {n: fiducial[n] for n in reference_model.parameter_names}
    )
    score = model.analytic_score(events, model.fiducial)

    print("\n  Per-event score agreement (median |dS| / column sd):")
    for name in list(MASS_NAMES) + REDSHIFT_NAMES + list(COSMOLOGY_NAMES):
        if name == "delta_m":
            continue
        reference_column = reference_score[
            :, reference_model.parameter_names.index(name)
        ]
        if name == "alpha":
            column = (
                score[:, model.parameter_names.index("gamma")]
                + score[:, model.parameter_names.index("kappa")]
            )
        elif name == "beta":
            column = score[:, model.parameter_names.index("kappa")]
        else:
            gwforge = MASS_NAMES.get(name) or COSMOLOGY_NAMES.get(name) or name
            column = score[:, model.parameter_names.index(gwforge)]
        first = reference_column - reference_column.mean()
        second = column - column.mean()
        spread = max(first.std(), 1e-300)
        print(
            "    {:<10s} {:9.2e}   (99th pct {:8.2e}, worst {:8.2e})".format(
                name,
                numpy.median(numpy.abs(first - second)) / spread,
                numpy.percentile(numpy.abs(first - second), 99) / spread,
                numpy.abs(first - second).max() / spread,
            )
        )
    print(
        "\n  The mass and redshift columns are identical to roundoff; only Om and\n"
        "  w0 -- the reference's finite-differenced ones -- disagree, and only on\n"
        "  a sparse minority of events. That is what widens every sigma above."
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-directory", default=DEFAULT_REFERENCE_DIRECTORY)
    parser.add_argument("--snr-thresholds", type=float, nargs="+", default=[10.0, 20.0])
    parser.add_argument(
        "--output-directory",
        default="population_fisher_comparison",
        help="Where the overlay figures go. Pass '' to skip plotting.",
    )
    parser.add_argument(
        "--skip-spectral-sirens",
        action="store_true",
        help="Mass and redshift only; the spectral-siren comparison is the slow one.",
    )
    options = parser.parse_args(argv)

    if options.output_directory:
        os.makedirs(options.output_directory, exist_ok=True)
    else:
        options.output_directory = None

    load, forecast_mass_redshift, spectral_sirens_forecast = reference_modules(
        options.reference_directory
    )
    print("Re-deriving the reference's mass fiducial from its injected masses ...")
    fiducial = reference_fiducial(
        options.reference_directory, forecast_mass_redshift, spectral_sirens_forecast
    )

    outcomes = []
    for threshold in options.snr_thresholds:
        outcomes.append(
            compare_mass_and_redshift(
                options.reference_directory,
                load,
                fiducial,
                threshold,
                options.output_directory,
            )
        )
        if not options.skip_spectral_sirens:
            outcomes.append(
                compare_spectral_sirens(
                    options.reference_directory,
                    load,
                    fiducial,
                    threshold,
                    options.output_directory,
                )
            )

    print("\n" + "=" * 72)
    print(
        "  {}/{} comparisons agree".format(sum(1 for ok in outcomes if ok), len(outcomes))
    )
    print(
        "  spin is not compared: the reference defines no spin model, only a "
        "sampler.\n  It keeps the finite-difference oracle in "
        "tests/test_population_fisher.py."
    )
    print("=" * 72)
    return 0 if all(outcomes) else 1


if __name__ == "__main__":
    raise SystemExit(main())
