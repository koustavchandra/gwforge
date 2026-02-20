"""
Sampling function for NotchFilterBinnedPairingMassDistribution model.

This module provides a function to sample binary masses according to the 
NotchFilterBinnedPairingMassDistribution population model.
"""

import sys
import numpy
from gwpopulation.utils import powerlaw
from tqdm import tqdm
import configparser
from pathlib import Path
import bilby

# Import from GWForge
from GWForge.population.pdb_external import NotchFilterBinnedPairingMassDistribution


def rejection_sampling_uniform_grid(
    n_samples,
    A, A2, NSmin, NSmax, BHmin, BHmax,
    UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5,
    alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
    beta_pair_1, beta_pair_2, mbreak,
    mmin=0.5, mmax=350.0,
    max_iterations=10000,
    verbose=False
):
    """
    Sample binary masses from the NotchFilterBinnedPairingMassDistribution model
    using rejection sampling.
    
    Parameters
    ----------
    n_samples : int
        Number of samples to draw
    A : float
        Depth of the lower mass gap dip
    A2 : float
        Depth of the upper mass gap dip
    NSmin : float
        Minimum mass (lower cutoff)
    NSmax : float
        Start of lower mass gap
    BHmin : float
        End of lower mass gap
    BHmax : float
        Maximum mass in power-law component
    UPPERmin : float
        Start of upper mass gap
    UPPERmax : float
        End of upper mass gap
    n0 : float
        Sharpness of low mass cutoff
    n1 : float
        Sharpness of lower edge of lower mass gap
    n2 : float
        Sharpness of upper edge of lower mass gap
    n3 : float
        Sharpness of lower edge of upper mass gap
    n4 : float
        Sharpness of upper edge of upper mass gap
    n5 : float
        Sharpness of high mass cutoff
    alpha_1 : float
        Power-law exponent for m < NSmax
    alpha_2 : float
        Power-law exponent for m > BHmin
    alpha_dip : float
        Power-law exponent between NSmax and BHmin
    mu1 : float
        Mean of first Gaussian peak
    sig1 : float
        Width of first Gaussian peak
    mix1 : float
        Mixing fraction of first Gaussian peak
    mu2 : float
        Mean of second Gaussian peak
    sig2 : float
        Width of second Gaussian peak
    mix2 : float
        Mixing fraction of second Gaussian peak
    beta_pair_1 : float
        Mass ratio exponent for m2 < mbreak
    beta_pair_2 : float
        Mass ratio exponent for m2 >= mbreak
    mbreak : float
        Mass breakpoint for pairing function transition
    mmin : float, optional
        Minimum mass (default: 0.5)
    mmax : float, optional
        Maximum mass (default: 350.0)
    max_iterations : int, optional
        Maximum number of rejection sampling iterations (default: 1000)
    verbose : bool, optional
        Print progress information (default: False)
        
    Returns
    -------
    m1_samples : ndarray
        Sampled primary masses
    m2_samples : ndarray
        Sampled secondary masses
    acceptance_rate : float
        Overall acceptance rate of the rejection sampling
    """
    
    # Initialize the model class
    model = NotchFilterBinnedPairingMassDistribution(mmin=mmin, mmax=mmax)
    
    # Hyperparameters dict for the model
    hyperparams = {
        'A': A, 'A2': A2, 'NSmin': NSmin, 'NSmax': NSmax, 'BHmin': BHmin, 'BHmax': BHmax,
        'UPPERmin': UPPERmin, 'UPPERmax': UPPERmax,
        'n0': n0, 'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4, 'n5': n5,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_dip': alpha_dip,
        'mu1': mu1, 'sig1': sig1, 'mix1': mix1, 'mu2': mu2, 'sig2': sig2, 'mix2': mix2,
        'beta_pair_1': beta_pair_1, 'beta_pair_2': beta_pair_2, 'mbreak': mbreak,
    }
    
    # Initialize storage for accepted samples
    m1_accepted = []
    m2_accepted = []
    
    # Estimate maximum of the probability density for rejection sampling
    # Sample a large grid to find approximate maximum
    test_m1 = numpy.linspace(mmin, mmax, 100)
    test_m2 = numpy.linspace(mmin, mmax, 100)
    test_m1_grid, test_m2_grid = numpy.meshgrid(test_m1, test_m2)
    
    # Enforce m1 >= m2 constraint
    valid = test_m1_grid >= test_m2_grid
    test_m1_flat = test_m1_grid[valid]
    test_m2_flat = test_m2_grid[valid]
    
    # Create dataset dict for model evaluation
    test_dataset = {
        'mass_1': test_m1_flat,
        'mass_2': test_m2_flat,
    }
    
    # Compute joint probability using the model
    joint_prob = model(test_dataset, **hyperparams)
    
    p_max = numpy.max(joint_prob)
    if p_max <= 0:
        raise ValueError("Maximum probability is non-positive. Check hyperparameters.")
    
    if verbose:
        print(f"Estimated maximum probability: {p_max}")
        print(f"Starting rejection sampling for {n_samples} samples...")
    
    total_proposals = 0
    iteration = 0
    
    pbar = tqdm(total=n_samples, desc="Sampling", disable=not verbose, unit="samples")
    
    while len(m1_accepted) < n_samples and iteration < max_iterations:
        # Number of proposals to draw in this iteration
        n_proposals = max(100, int(2 * (n_samples - len(m1_accepted))))
        
        # Sample m1 uniformly from [mmin, mmax]
        m1_proposal = numpy.random.uniform(mmin, mmax, n_proposals)
        
        # Sample m2 uniformly from [mmin, m1] to ensure m1 >= m2
        m2_proposal = numpy.random.uniform(mmin, m1_proposal)
        
        # Create dataset dict for model evaluation
        dataset = {
            'mass_1': m1_proposal,
            'mass_2': m2_proposal,
        }
        
        # Evaluate joint probability using the model
        joint_prob = model(dataset, **hyperparams)
        
        # Acceptance test with uniform random numbers
        u = numpy.random.uniform(0, 1, len(m1_proposal))
        acceptance_threshold = joint_prob / p_max
        accepted = u < acceptance_threshold
        
        # Store accepted samples
        m1_accepted.extend(m1_proposal[accepted])
        m2_accepted.extend(m2_proposal[accepted])
        
        n_accepted = numpy.sum(accepted)
        total_proposals += len(m1_proposal)
        iteration += 1
        
        pbar.update(min(n_accepted, n_samples - pbar.n))
    
    pbar.close()
    
    if len(m1_accepted) < n_samples and verbose:
        print(f"Warning: Could only generate {len(m1_accepted)} out of {n_samples} requested samples "
              f"after {max_iterations} iterations.")
    
    # Convert to arrays - return however many samples were obtained
    m1_samples = numpy.array(m1_accepted)
    m2_samples = numpy.array(m2_accepted)
    
    actual_samples = len(m1_samples)
    acceptance_rate = actual_samples / total_proposals if total_proposals > 0 else 0
    
    if verbose:
        print(f"Sampling complete!")
        print(f"Generated {actual_samples} samples (requested {n_samples})")
        print(f"Total proposals: {total_proposals}")
        print(f"Acceptance rate: {acceptance_rate:.4f}")
    
    return m1_samples, m2_samples, acceptance_rate

def importance_sampling_m1_m2_prop(
    n_samples,
    A, A2, NSmin, NSmax, BHmin, BHmax,
    UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5,
    alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
    beta_pair_1, beta_pair_2, mbreak,
    mmin=0.5, mmax=350.0,
    oversample_factor=5,
    verbose=False
):
    """
    Sample binary masses using importance sampling instead of rejection sampling.
    """

    model = NotchFilterBinnedPairingMassDistribution(mmin=mmin, mmax=mmax)

    hyperparams = {
        'A': A, 'A2': A2, 'NSmin': NSmin, 'NSmax': NSmax, 'BHmin': BHmin, 'BHmax': BHmax,
        'UPPERmin': UPPERmin, 'UPPERmax': UPPERmax,
        'n0': n0, 'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4, 'n5': n5,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_dip': alpha_dip,
        'mu1': mu1, 'sig1': sig1, 'mix1': mix1, 'mu2': mu2, 'sig2': sig2, 'mix2': mix2,
        'beta_pair_1': beta_pair_1, 'beta_pair_2': beta_pair_2, 'mbreak': mbreak,
    }

    # ---- Build 1D proposal p(m) ----
    hyperpars_subset = {k: hyperparams[k] for k in [
        'alpha_1', 'alpha_dip', 'alpha_2',
        'NSmin', 'NSmax', 'BHmin', 'BHmax',
        'UPPERmin', 'UPPERmax',
        'n0', 'n1', 'n2', 'n3', 'n4', 'n5',
        'mix1', 'mu1', 'sig1', 'mix2', 'mu2', 'sig2',
        'A', 'A2'
    ]}

    m_input = numpy.linspace(mmin, mmax, int(1e4))
    p_m = model.p_m(m_input, **hyperpars_subset)

    prob_m = bilby.core.prior.Interped(
        m_input, p_m,
        minimum=mmin,
        maximum=mmax
    )

    # ---- Draw proposals ----
    N_prop = oversample_factor * n_samples

    if verbose:
        print(f"Drawing {N_prop} proposal samples...")

    m1_prop = prob_m.sample(N_prop)
    m2_prop = prob_m.sample(N_prop)

    # enforce m1 >= m2
    valid = m1_prop >= m2_prop
    m1_prop = m1_prop[valid]
    m2_prop = m2_prop[valid]

    if len(m1_prop) < n_samples:
        raise RuntimeError("Not enough valid proposals — increase oversample_factor")

    dataset = {
        'mass_1': m1_prop,
        'mass_2': m2_prop,
    }

    # ---- Evaluate true joint PDF ----
    joint_prob = model(dataset, **hyperparams)

    # ---- Evaluate proposal density q(m1,m2) ----
    q_m1 = prob_m.prob(m1_prop)
    q_m2 = prob_m.prob(m2_prop)
    proposal_prob = q_m1 * q_m2

    # Avoid division by zero
    mask = proposal_prob > 0
    joint_prob = joint_prob[mask]
    proposal_prob = proposal_prob[mask]
    m1_prop = m1_prop[mask]
    m2_prop = m2_prop[mask]

    weights = joint_prob / proposal_prob
    weights /= numpy.sum(weights)

    # ---- Resample according to weights ----
    rng = numpy.random.default_rng()
    indices = rng.choice(len(weights), size=n_samples, replace=True, p=weights)

    m1_samples = m1_prop[indices]
    m2_samples = m2_prop[indices]

    # ---- Effective sample size ----
    ess = 1.0 / numpy.sum(weights**2)

    if verbose:
        print(f"Importance sampling complete.")
        print(f"Effective sample size ≈ {ess:.1f} / {len(weights)}")

    return m1_samples, m2_samples, ess

def importance_sampling_m1_q_prop(
    n_samples,
    A, A2, NSmin, NSmax, BHmin, BHmax,
    UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5,
    alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
    beta_pair_1, beta_pair_2, mbreak,
    mmin=0.5, mmax=350.0,
    oversample_factor=3,
    verbose=False
):
    """
    Importance sampling using near-optimal proposal:
        m1 ~ p_m(m)
        q  ~ pairing power-law
        m2 = q * m1
    """

    model = NotchFilterBinnedPairingMassDistribution(mmin=mmin, mmax=mmax)

    hyperparams = {
        'A': A, 'A2': A2, 'NSmin': NSmin, 'NSmax': NSmax,
        'BHmin': BHmin, 'BHmax': BHmax,
        'UPPERmin': UPPERmin, 'UPPERmax': UPPERmax,
        'n0': n0, 'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4, 'n5': n5,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_dip': alpha_dip,
        'mu1': mu1, 'sig1': sig1, 'mix1': mix1,
        'mu2': mu2, 'sig2': sig2, 'mix2': mix2,
        'beta_pair_1': beta_pair_1,
        'beta_pair_2': beta_pair_2,
        'mbreak': mbreak,
    }

    # -------------------------
    # 1D mass proposal p(m)
    # -------------------------
    hyperpars_subset = {k: hyperparams[k] for k in [
        'alpha_1', 'alpha_dip', 'alpha_2',
        'NSmin', 'NSmax', 'BHmin', 'BHmax',
        'UPPERmin', 'UPPERmax',
        'n0', 'n1', 'n2', 'n3', 'n4', 'n5',
        'mix1', 'mu1', 'sig1', 'mix2', 'mu2', 'sig2',
        'A', 'A2'
    ]}

    m_grid = numpy.linspace(mmin, mmax, 10000)
    p_m_grid = model.p_m(m_grid, **hyperpars_subset)

    prob_m = bilby.core.prior.Interped(
        m_grid, p_m_grid,
        minimum=mmin,
        maximum=mmax
    )

    # -------------------------
    # Proposal sampling
    # -------------------------
    N_prop = oversample_factor * n_samples
    rng = numpy.random.default_rng()

    m1_prop = prob_m.sample(N_prop)

    # ---- Sample mass ratio q ----
    q_prop = numpy.empty(N_prop)

    for i in range(N_prop):

        m1 = m1_prop[i]

        # q_min = mmin / m1
        q_min = mmin / m1
        q_max = 1.0

        # choose correct beta
        beta = beta_pair_1 if m1 < mbreak else beta_pair_2

        # sample power-law in q
        if abs(beta + 1) > 1e-8:
            u = rng.uniform()
            q_prop[i] = (
                (u * (q_max**(beta+1) - q_min**(beta+1)) + q_min**(beta+1))
                ** (1/(beta+1))
            )
        else:
            # beta = -1 special case
            u = rng.uniform()
            q_prop[i] = q_min * (q_max / q_min)**u

    m2_prop = q_prop * m1_prop

    # -------------------------
    # Evaluate target density
    # -------------------------
    dataset = {'mass_1': m1_prop, 'mass_2': m2_prop}
    target = model(dataset, **hyperparams)

    # -------------------------
    # Evaluate proposal density
    # -------------------------

    # p(m1)
    q_m1 = prob_m.prob(m1_prop)

    # p(q | m1)
    q_pdf = numpy.zeros_like(q_prop)

    for i in tqdm(range(N_prop)):
        m1 = m1_prop[i]
        beta = beta_pair_1 if m1 < mbreak else beta_pair_2

        q_min = mmin / m1
        q_max = 1.0

        if abs(beta + 1) > 1e-8:
            norm = (q_max**(beta+1) - q_min**(beta+1)) / (beta+1)
            q_pdf[i] = q_prop[i]**beta / norm
        else:
            norm = numpy.log(q_max / q_min)
            q_pdf[i] = 1.0 / (q_prop[i] * norm)

    # Jacobian: m2 = q*m1  →  dm2 = m1 dq
    proposal = q_m1 * q_pdf / m1_prop

    mask = proposal > 0
    target = target[mask]
    proposal = proposal[mask]
    m1_prop = m1_prop[mask]
    m2_prop = m2_prop[mask]

    weights = target / proposal
    weights /= numpy.sum(weights)

    # -------------------------
    # Resample
    # -------------------------
    indices = rng.choice(len(weights), size=n_samples, replace=True, p=weights)

    m1_samples = m1_prop[indices]
    m2_samples = m2_prop[indices]

    ess = 1.0 / numpy.sum(weights**2)

    if verbose:
        print(f"Effective sample size ≈ {ess:.1f} / {len(weights)}")

    return m1_samples, m2_samples, ess

def lintsampling(
    n_samples,
    A, A2, NSmin, NSmax, BHmin, BHmax,
    UPPERmin, UPPERmax, n0, n1, n2, n3, n4, n5,
    alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2,
    beta_pair_1, beta_pair_2, mbreak,
    mmin=0.5, mmax=350.0,
    grid_size=256,
    verbose=False
):
    """
    Sample binary masses using LintSampler.
    """
    try:
        from lintsampler import LintSampler
    except ImportError:
        raise ImportError("Please install lintsampler using pip")

    model = NotchFilterBinnedPairingMassDistribution(mmin=mmin, mmax=mmax)

    hyperparams = {
        'A': A, 'A2': A2, 'NSmin': NSmin, 'NSmax': NSmax,
        'BHmin': BHmin, 'BHmax': BHmax,
        'UPPERmin': UPPERmin, 'UPPERmax': UPPERmax,
        'n0': n0, 'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4, 'n5': n5,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_dip': alpha_dip,
        'mu1': mu1, 'sig1': sig1, 'mix1': mix1,
        'mu2': mu2, 'sig2': sig2, 'mix2': mix2,
        'beta_pair_1': beta_pair_1,
        'beta_pair_2': beta_pair_2,
        'mbreak': mbreak,
    }

    # -------------------------------------------------
    # Define vectorized joint PDF for LintSampler
    # -------------------------------------------------
    def joint_pdf(x):
        """
        x : array of shape (N,2)
        """
        m1 = x[:, 0]
        m2 = x[:, 1]

        # Enforce physical constraint
        valid = m1 >= m2

        pdf = numpy.zeros(len(x))

        if numpy.any(valid):
            dataset = {
                'mass_1': m1[valid],
                'mass_2': m2[valid],
            }
            pdf[valid] = model(dataset, **hyperparams)

        return pdf

    # -------------------------------------------------
    # Build 2D domain grid
    # -------------------------------------------------
    m1_axis = numpy.geomspace(mmin, mmax, grid_size)
    m2_axis = numpy.geomspace(mmin, mmax, grid_size)

    domain = (m1_axis, m2_axis)

    if verbose:
        print(f"Building LintSampler grid: {grid_size} x {grid_size}")

    sampler = LintSampler(
        domain,
        joint_pdf,
        seed=numpy.random.default_rng(),
        vectorizedpdf=True
    )

    samples = sampler.sample(n_samples)

    m1_samples = samples[:, 0]
    m2_samples = samples[:, 1]

    return m1_samples, m2_samples
