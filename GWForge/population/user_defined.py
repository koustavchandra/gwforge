"""User-defined population distributions from tabulated (xx, yy) priors.

This lets people feed a population produced by an external population-synthesis
code straight into GWForge, without shoehorning it into a parameterised model.
The distribution of any parameter is supplied as a support array ``xx`` and a
(possibly unnormalised) probability-density array ``yy`` in a JSON file with the
structure::

    {
      "mass_1_source": {"xx": [...], "yy": [...]},
      "mass_ratio":    {"xx": [...], "yy": [...]},
      "a_1":           {"xx": [...], "yy": [...]}
    }

Each ``(xx, yy)`` pair is turned into a :class:`bilby.core.prior.Interped`
prior, which normalises ``yy`` and samples via inverse-CDF.
"""

import json
import numpy
import bilby


def get_xx_yy_from_population_priors(file, parameter):
    """Extract the ``xx`` (support) and ``yy`` (density) arrays for a parameter.

    Parameters
    ----------
    file : str
        Path to the JSON file of population priors, structured as
        ``{parameter: {"xx": [...], "yy": [...]}}``.
    parameter : str
        The parameter name to extract (e.g. ``'mass_1_source'``, ``'mass_ratio'``,
        ``'a_1'``).

    Returns
    -------
    xx : numpy.ndarray
        The support (x-values) of the prior distribution.
    yy : numpy.ndarray
        The probability density (y-values) of the prior distribution.
    """
    with open(file, "r") as f:
        priors = json.load(f)
    if parameter not in priors:
        raise KeyError(
            "Parameter '{}' not found in {}. Available: {}".format(
                parameter, file, sorted(priors)
            )
        )
    xx = numpy.asarray(priors[parameter]["xx"], dtype=float)
    yy = numpy.asarray(priors[parameter]["yy"], dtype=float)
    return xx, yy


def interped_prior_from_file(file, parameter, name=None):
    """Build a :class:`bilby.core.prior.Interped` prior from a tabulated prior.

    Parameters
    ----------
    file : str
        Path to the JSON population-priors file.
    parameter : str
        Key to read from the file.
    name : str, optional
        Name to give the resulting prior. Defaults to ``parameter``.

    Returns
    -------
    bilby.core.prior.Interped
    """
    xx, yy = get_xx_yy_from_population_priors(file, parameter)
    if numpy.any(yy < 0):
        raise ValueError(
            "Probability density 'yy' for '{}' contains negative values".format(
                parameter
            )
        )
    return bilby.core.prior.Interped(
        xx=xx,
        yy=yy,
        minimum=numpy.min(xx),
        maximum=numpy.max(xx),
        name=name or parameter,
    )
