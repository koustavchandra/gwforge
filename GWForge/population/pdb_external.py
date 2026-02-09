""" Power Law + Dip + Break models. None are reviewed."""
from gwpopulation.utils import powerlaw, truncnorm
from pairing import _IdenticalPairingMassDistribution

def power_law_dip_break_1d(mass, 
                           A, A2, NSmin, NSmax, BHmin, BHmax, UPPERmin, UPPERmax,
                           n0, n1, n2, n3, n4, n5, 
                           alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
                           absolute_mmin, absolute_mmax
):
    r"""
    The one-dimensional mass distribution considered in Fishbach, Essick, Holz. Does
    Matter Matter? ApJ Lett 899, 1 (2020) : arXiv:2006.13178

    .. math::
        p(m|\lambda) = n(m|\gamma_{\text{low}}, \gamma_{\text{high}}, A) \times
            l(m|m_{\text{max}}, \eta) \\
                 \times \begin{cases}
                         & m^{\alpha_1} \text{ if } m < \gamma_{\text{low}} \\
                         & m^{\alpha_2} \text{ if } m > \gamma_{\text{low}} \\
                         & 0 \text{ otherwise }
                 \end{cases}.
    
    where $l(m|m_{\text{max}}, \eta)$ is the low pass filter with powerlaw $\eta$
    applied at mass $m_{\text{max}}$,
    $n(m|\gamma_{\text{low}}, \gamma_{\text{high}}, A)$ is the notch
    filter with depth $A$ applied between $\gamma_{\text{low}}$ and 
    $\gamma_{\text{high}}$, and
    $\lambda$ is the subset of hyperparameters $\{ \gamma_{\text{low}},
    \gamma_{\text{high}}, A, \alpha_1, \alpha_2, m_{\text{min}}, m_{\
    text{max}}\}$.
    
    Parameters
    ----------
    mass: array-like
        Mass to evaluate probability at (:math:`m`).
    alpha_1: float
        Powerlaw exponent for compact objects below NSmax (:math:`\alpha_1`).
    alpha_2: float
        Powerlaw exponent for compact objects above BHmin (:math:`\alpha_2`).
    alpha_dip: float
        Powerlaw exponent for compact objects between NSmax and BHmin (:math: `\alpha_d`).
    NSmin: float
        Minimum compact object mass (:math:`m_\min`).
    NSmax: float
        Mass at which the notch filter starts for the lower mass gap (:math:`\gamma_{low1}`).
    BHmin: float
        Mass at which the notch filter ends for the lower mass gap (:math:`\gamma_{high1}`).
    BHmax: float
        Maximum mass in the powerlaw distributed component (:math:`m_\max`).
    A: float
        depth of the dip between NSmax and BHmin (A).
    UPPERmin: float
        Mass at which the notch filter starts for the upper mass gap (:math:`\gamma_{low2}`).
    UPPERmax: float
        Mass at which the notch filter ends for the upper mass gap (:math:`\gamma_{high2}`).
    mu1: float
        Location of the upper peak where an overdensity of merging compact objects is observed (:math:`\mu_{peak1}`).. 
    sig1: float
        Width of the upper peak where an overdensity of merging compact objects is observed (:math:`\sigma_{peak1}`)..
    mix1: float
        Mixing fraction of the first gaussian peak with the powerlaw + notches (:math: `c_1`)
    mu2: float
        Location of the lower peak where an overdensity of merging compact objects is observed (:math:`\mu_{peak2}`)..
    sig2: float
        Width of the lower peak where an overdensity of merging compact objects is observed (:math:`\mu_{peak2}`)..
    mix2: float
        Mixing fraction of the second gaussian peak with the powerlaw + notches (:math: `c_2`)
    absolute_mmin: float
        The minimum limit for the truncated normal distribution. 
    absolute_mmax: float
        The maximum limit for the truncated normal distribution. 
    n{0,5}:float
        Exponents to set the sharpness of the low mass cutoff and high mass cutoff, respectively (:math:`\eta_i`). 
    n{1,2}: float
        Exponents to set the sharpness of the lower edge and upper edge of the lower mass gap, respectively (:math:`\eta_i`). 
    n{3,4}: float
        Exponents to set the sharpness of the lower edge and upper edge of the upper mass gap, respectively (:math:`\eta_i`). 
        
    """
    from gwpopulation.utils import xp
    
    gaussian_peak1 = truncnorm(mass, mu1, sig1,low=absolute_mmin,high=absolute_mmax)
    gaussian_peak2 = truncnorm(mass, mu2, sig2,low=absolute_mmin,high=absolute_mmax)
    
    condlist = [mass<NSmax,(mass>=NSmax) & (mass<BHmin),mass>=BHmin]
    choicelist = [mass**alpha_1,
                  (mass**alpha_dip)*(NSmax**(alpha_1-alpha_dip)),
                  (mass**alpha_2)*(NSmax**(alpha_1-alpha_dip)) * (BHmin**(alpha_dip-alpha_2))
                  ]
    plaw = xp.select(condlist,choicelist,default=0.0)

    highpass_lower = (1 + (NSmin / mass) ** n0)
    notch_lower = 1.0 - A / ((1 + (NSmax / mass) ** n1) * (1 + (mass / BHmin) ** n2))
    notch_upper = 1.0 - A2 / ((1 + (UPPERmin / mass) ** n3) * (1 + (mass / UPPERmax) ** n4))
    lowpass_upper = 1 + (mass / BHmax) ** n5

    return (1 + mix1 * gaussian_peak1 + mix2 * gaussian_peak2) * plaw * notch_lower * notch_upper / highpass_lower / lowpass_upper

class _NotchFilterPairingMassDistribution(_IdenticalPairingMassDistribution):
    """
    Generic pairing function mass distribution with "matter matters" 1D mass
    model base class.
    """
    def p_m(self, mass, A, A2, NSmin, NSmax, BHmin, BHmax, 
            UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
            alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2):
        return power_law_dip_break_1d(mass, A, A2, NSmin, NSmax, BHmin, BHmax, 
                                      UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
                                      alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, 
                                      mu2, sig2, mix2,absolute_mmin=self.mmin,absolute_mmax=self.mmax)

class NotchFilterPowerLawPairingMassDistribution(_NotchFilterPairingMassDistribution):
    """
    Two-dimensional "Power Law + Dip + Break" Model.

    Parameters
    ----------
    dataset: dict
        Dictionary of numpy arrays for 'mass_1' (:math:`m_1`) and
        'mass_ratio' q (:math:`m_2=m_1*q`).
    alpha_1: float
        Powerlaw exponent for compact object below break (:math:`\alpha_1`).
    alpha_2: float
        Powerlaw exponent for compact object above break (:math:`\alpha_2`).
    NSmin: float
        Minimum compact object mass (:math:`m_\min`).
    NSmax: float
        Mass at which the notch filter starts (:math:`\gamma_{low}`)
    BHmin: float
        Mass at which the notch filter ends (:math:`\gamma_{high}`).
        Also, the mass value at which the power law exponent switches from alpha_1 to alpha_2.
    BHmax: float
        Maximum mass in the powerlaw distributed component (:math:`m_\max`).
    n{0,1,2,3}: float
        Exponents to set the sharpness of the low mass cutoff, low edge of dip,
        high edge of dip, and high mass cutoff, respectively (:math:`\eta_i`).
    A: float
        depth of the dip between NSmax and BHmin (A).
    """
    def pairing(self, dataset, beta_q):
        mass_ratio = dataset["mass_2"]/dataset["mass_1"]
        return powerlaw(mass_ratio, beta_q, 1, self.qmin)

    def __call__(self, dataset, A, A2, NSmin, NSmax, BHmin, BHmax, 
            UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
            alpha_1, alpha_2, mu1, sig1, mix1, mu2, sig2, mix2, beta_q
            ):
        # get arguments in a dict
        kwargs = locals()
        kwargs.pop('self')
        return self.p_m1_m2(**kwargs)
    
class NotchFilterBinnedPairingMassDistribution(_NotchFilterPairingMassDistribution):
    """
    2D "Power Law + Dip + Break" Model with a pairing function that depends on m2.
    """
    def pairing(self, dataset, beta_pair_1, beta_pair_2, mbreak):
        from gwpopulation.utils import xp
        mass_ratio = dataset["mass_2"]/dataset["mass_1"]

        beta_pair = xp.where(dataset["mass_2"] < mbreak, beta_pair_1, beta_pair_2)
        lower_q_bound = xp.where(dataset["mass_2"] < mbreak,
            self.qmin, mbreak/self.mmax)

        return powerlaw(mass_ratio, beta_pair, 1, lower_q_bound)

    def __call__(self, dataset, A, A2, NSmin, NSmax, BHmin, BHmax, 
            UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
            alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
            beta_pair_1, beta_pair_2, mbreak,
            ):
        # get arguments in a dict
        kwargs = locals()
        kwargs.pop('self')
        return self.p_m1_m2(**kwargs)

class NotchFilterBinnedPairing2MassDistribution(_NotchFilterPairingMassDistribution):
    """
    2D "Power Law + Dip + Break" Model with a pairing function that depends on m2. - bin edge is Nsmax
    """
    def pairing(self, dataset, beta_pair_1, beta_pair_2, NSmax):
        from gwpopulation.utils import xp
        mass_ratio = dataset["mass_2"]/dataset["mass_1"]

        beta_pair = xp.where(dataset["mass_2"] < NSmax, beta_pair_1, beta_pair_2)
        lower_q_bound = xp.where(dataset["mass_2"] < NSmax,
            self.qmin, NSmax/self.mmax)
        return powerlaw(mass_ratio, beta_pair, 1, lower_q_bound)

    def __call__(self, dataset, A, A2, NSmin, NSmax, BHmin, BHmax, 
            UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
            alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
            beta_pair_1, beta_pair_2
            ):
        # get arguments in a dict
        kwargs = locals()
        kwargs.pop('self')
        return self.p_m1_m2(**kwargs)
    
class NotchFilterBinnedPairing3MassDistribution(_NotchFilterPairingMassDistribution):
    """
    2D "Power Law + Dip + Break" Model with a pairing function that depends on m2. - bin edge is BHmin
    """
    def pairing(self, dataset, beta_pair_1, beta_pair_2, BHmin):
        from gwpopulation.utils import xp
        mass_ratio = dataset["mass_2"]/dataset["mass_1"]

        beta_pair = xp.where(dataset["mass_2"] < BHmin, beta_pair_1, beta_pair_2)
        lower_q_bound = xp.where(dataset["mass_2"] < BHmin,
            self.qmin, BHmin/self.mmax)
        return powerlaw(mass_ratio, beta_pair, 1, lower_q_bound)

    def __call__(self, dataset, A, A2, NSmin, NSmax, BHmin, BHmax, 
            UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5, 
            alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
            beta_pair_1, beta_pair_2
            ):
        # get arguments in a dict
        kwargs = locals()
        kwargs.pop('self')
        return self.p_m1_m2(**kwargs)
