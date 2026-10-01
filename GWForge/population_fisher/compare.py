"""Fit a population block to the injections and report what came back.

Two different jobs wear the same clothes here, and the difference matters when
reading the output.

For the **redshift** block the fit is load-bearing.
:class:`GWForge.population.redshift.Redshift` convolves the star-formation rate
with a formation-to-merger time-delay distribution, so a generated catalogue is
not distributed as the direct Madau-Dickinson form. The fitted
:math:`(\\gamma, \\kappa, z_p)` are the *effective* parameters of the model being
forecast, and they will not equal the generation values -- that is expected, and
the ratio column says how far the time delay moved them.

For the **mass** and **spin** blocks the fit is a check.
``gwforge_population`` draws those straight from the same densities this package
evaluates, so the fit must return the generation values to within Monte-Carlo
error. A ratio that is not close to one there is a bug, not a result.
"""

import logging

import numpy

from .fit import compare_to, fit_model

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

# How each block's fit should be read.
BLOCK_ROLE = {
    "redshift": (
        "fiducial (the generator convolves psi_MD with a time delay, so the "
        "generation values are not this model's fiducials)"
    ),
    "mass": "consistency check (should return the generation values)",
    "spin": "consistency check (should return the generation values)",
}


def compare_block(model, block, injections, free_parameters=None):
    """Fit one block to the injections and format the comparison.

    Parameters
    ----------
    model : GWForge.population_fisher.base.PopulationModel
        Fitted **in place**: on return ``model.fiducial`` holds the fitted
        values, so the Fisher that follows is evaluated at them.
    block : str
        ``"mass"``, ``"redshift"`` or ``"spin"``.
    injections : dict
        Every injection, detected or not.
    free_parameters : sequence of str or None
        Vary only these, holding the rest at their config values. Defaults to
        every parameter the model has.

        The driver passes the *Fisher's* free list, so the fit varies exactly
        what the forecast will vary. That matters when a parameter is pinned
        because the data cannot constrain it: ``delta_m`` for a taper narrower
        than the sampling grid has no maximum to find, so an optimiser given it
        rails against whichever bound it was handed and drags the correlated
        parameters with it. Since the catalogue was *generated* here, such
        parameters are known exactly and there is nothing to fit.

    Returns
    -------
    list of str
        Lines for the summary file.
    """
    generated = dict(model.fiducial)
    count = len(numpy.atleast_1d(injections[model.event_keys[0]]))
    if free_parameters is not None:
        free_parameters = [
            name for name in free_parameters if name in model.parameter_names
        ]
        if not free_parameters:
            free_parameters = None
    logging.info(
        "Fitting the %s block to %d injections over %s",
        block,
        count,
        free_parameters if free_parameters is not None else "all parameters",
    )
    fitted = fit_model(model, injections, free_parameters=free_parameters)
    for name in model.parameter_names:
        model.fiducial[name] = float(fitted[name])
    varied = free_parameters or model.parameter_names
    held = [name for name in model.parameter_names if name not in varied]
    lines = [
        "{} block fit to {} injections -- {}".format(
            block, count, BLOCK_ROLE.get(block, "")
        ),
        compare_to(fitted, generated, varied),
        "  mean log p = {:.6f} over {} samples{}".format(
            fitted["_mean_log_prob"],
            fitted["_n_samples"],
            "" if fitted["_success"] else "   (optimiser: {})".format(fitted["_message"]),
        ),
    ]
    if held:
        lines.append(
            "  held at their config values: "
            + ", ".join("{}={:g}".format(name, generated[name]) for name in held)
        )
    if fitted["_n_dropped"]:
        lines.append(
            "  WARNING: {} sample(s) were outside the model's support and were "
            "dropped from the fit.".format(fitted["_n_dropped"])
        )
    return lines
