"""Base classes to implement models with pairing functions."""

import inspect


def _primary_secondary_general(dataset, p_m1, p_m2):
    return p_m1 * p_m2 * (dataset["mass_1"] >= dataset["mass_2"]) * 2


class _PairingMassDistribution(object):
    r"""
    Generic mass distribution with a pairing function base class.

    Implements mass distributions of the form:

    .. math::
        p(m_1, m_2) = p_1(m_1) * p_2(m_2) * f_p(q) : m_1 \geq m_2

    """

    def __init__(self, mmin=0.5, mmax=350.0):
        self.mmin = mmin
        self.mmax = mmax
        self.qmin = mmin / mmax

    def __call__(self, dataset, **kwargs):
        raise NotImplementedError

    def p1_m1(self, *args, **kwargs):
        raise NotImplementedError

    def p2_m2(self, *args, **kwargs):
        raise NotImplementedError

    def pairing(self, *args, **kwargs):
        raise NotImplementedError

    def p_m1_m2(self, dataset, **kwargs):
        # parse arguments for use in pm(m) vs fp(q)
        from gwpopulation.utils import xp

        pm1_args = [k for k, v in inspect.signature(self.p1_m1).parameters.items() if k not in ["self", "mass"]]
        pm1_dict = {k: kwargs[k] for k in dict(kwargs) if k in pm1_args}

        pm2_args = [k for k, v in inspect.signature(self.p2_m2).parameters.items() if k not in ["self", "mass"]]
        pm2_dict = {k: kwargs[k] for k in dict(kwargs) if k in pm2_args}

        fp_args = [k for k, v in inspect.signature(self.pairing).parameters.items() if k not in ["self", "dataset"]]
        fp_dict = {k: kwargs[k] for k in dict(kwargs) if k in fp_args}

        p_m1 = xp.where((dataset["mass_1"] >= self.mmin) * (dataset["mass_1"] <= self.mmax), self.p1_m1(dataset["mass_1"], **pm1_dict), 0.0)
        p_m2 = xp.where((dataset["mass_2"] >= self.mmin) * (dataset["mass_2"] <= self.mmax), self.p2_m2(dataset["mass_2"], **pm2_dict), 0.0)
        fp = self.pairing(dataset, **fp_dict)

        return _primary_secondary_general(dataset, p_m1, p_m2) * fp


class _IdenticalPairingMassDistribution(_PairingMassDistribution):
    r"""
    Base class for mass distribution with a mass ratio-dependent pairing function,
    where p(m1) and p(m2) are identical.

    Implements mass distributions of the form:

    .. math::
        p(m_1, m_2) &= p_m(m_1) * p_m(m_2) * f_p(q) : m_1 \geq m_2

        q &= m_2/m_1
    """

    def __init__(self, mmin=0.5, mmax=350.0):
        self.mmin = mmin
        self.mmax = mmax
        self.qmin = mmin / mmax

    def __call__(self, dataset, **kwargs):
        raise NotImplementedError

    def p_m(self, *args, **kwargs):
        raise NotImplementedError

    def pairing(self, *args, **kwargs):
        raise NotImplementedError

    def p_m1_m2(self, dataset, **kwargs):
        from gwpopulation.utils import xp

        # parse arguments for use in pm(m) vs fp(q)
        pm_args = [k for k, v in inspect.signature(self.p_m).parameters.items() if k not in ["self", "mass"]]
        pm_dict = {k: kwargs[k] for k in dict(kwargs) if k in pm_args}

        fp_args = [k for k, v in inspect.signature(self.pairing).parameters.items() if k not in ["self", "dataset"]]
        fp_dict = {k: kwargs[k] for k in dict(kwargs) if k in fp_args}

        # evaluate probabilities
        p_m1 = xp.where((dataset["mass_1"] >= self.mmin) * (dataset["mass_1"] <= self.mmax), self.p_m(dataset["mass_1"], **pm_dict), 0.0)
        p_m2 = xp.where((dataset["mass_2"] >= self.mmin) * (dataset["mass_2"] <= self.mmax), self.p_m(dataset["mass_2"], **pm_dict), 0.0)
        fp = self.pairing(dataset, **fp_dict)

        return _primary_secondary_general(dataset, p_m1, p_m2) * fp
