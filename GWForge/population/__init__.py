import json
import numpy


def get_xx_yy_from_population_priors(file, parameter):
    """
    Extract xx (support) and yy (probability density) arrays for a given parameter
    from a population priors JSON file.

    Parameters
    ----------
    file : str
        Path to the JSON file containing population priors.
        The file should have structure: {parameter: {"xx": [...], "yy": [...]}}.
    parameter : str
        The parameter name to extract (e.g. 'mass_1', 'mass_ratio', 'a_1').

    Returns
    -------
    xx : numpy.ndarray
        The support (x-values) of the prior distribution.
    yy : numpy.ndarray
        The probability density (y-values) of the prior distribution.
    """
    with open(file, "r") as f:
        priors = json.load(f)
    xx = numpy.asarray(priors[parameter]["xx"])
    yy = numpy.asarray(priors[parameter]["yy"])
    return xx, yy
