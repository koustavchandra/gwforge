"""
Sampling function for NotchFilterBinnedPairingMassDistribution model.

This module provides a function to sample binary masses according to the 
NotchFilterBinnedPairingMassDistribution population model.
"""

import sys
import numpy as np
from gwpopulation.utils import powerlaw
from tqdm import tqdm
import configparser
from pathlib import Path

# Import from GWForge
from GWForge.population.pdb_external import NotchFilterBinnedPairingMassDistribution


def sample_notch_filter_binned_pairing_masses(
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
    test_m1 = np.linspace(mmin, mmax, 100)
    test_m2 = np.linspace(mmin, mmax, 100)
    test_m1_grid, test_m2_grid = np.meshgrid(test_m1, test_m2)
    
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
    
    p_max = np.max(joint_prob)
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
        m1_proposal = np.random.uniform(mmin, mmax, n_proposals)
        
        # Sample m2 uniformly from [mmin, m1] to ensure m1 >= m2
        m2_proposal = np.random.uniform(mmin, m1_proposal)
        
        # Create dataset dict for model evaluation
        dataset = {
            'mass_1': m1_proposal,
            'mass_2': m2_proposal,
        }
        
        # Evaluate joint probability using the model
        joint_prob = model(dataset, **hyperparams)
        
        # Acceptance test with uniform random numbers
        u = np.random.uniform(0, 1, len(m1_proposal))
        acceptance_threshold = joint_prob / p_max
        accepted = u < acceptance_threshold
        
        # Store accepted samples
        m1_accepted.extend(m1_proposal[accepted])
        m2_accepted.extend(m2_proposal[accepted])
        
        n_accepted = np.sum(accepted)
        total_proposals += len(m1_proposal)
        iteration += 1
        
        pbar.update(min(n_accepted, n_samples - pbar.n))
    
    pbar.close()
    
    if len(m1_accepted) < n_samples and verbose:
        print(f"Warning: Could only generate {len(m1_accepted)} out of {n_samples} requested samples "
              f"after {max_iterations} iterations.")
    
    # Convert to arrays - return however many samples were obtained
    m1_samples = np.array(m1_accepted)
    m2_samples = np.array(m2_accepted)
    
    actual_samples = len(m1_samples)
    acceptance_rate = actual_samples / total_proposals if total_proposals > 0 else 0
    
    if verbose:
        print(f"Sampling complete!")
        print(f"Generated {actual_samples} samples (requested {n_samples})")
        print(f"Total proposals: {total_proposals}")
        print(f"Acceptance rate: {acceptance_rate:.4f}")
    
    return m1_samples, m2_samples, acceptance_rate


# Example usage
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Load hyperparameters from INI file
    cfg = configparser.ConfigParser()
    # preserve key case from INI (case-sensitive keys)
    cfg.optionxform = str
    ini_path = Path(__file__).with_name('pdb_hyperpars.ini')
    if not ini_path.exists():
        raise FileNotFoundError(f"INI file not found: {ini_path}")
    cfg.read(ini_path)
    sec = cfg['hyperparameters']

    # Load entire INI section as float-valued hyperparameters dictionary
    hyperparams = {k: float(v) for k, v in sec.items()}
    print("Loaded hyperparameters from INI file:", hyperparams)

    # Control parameters (use defaults if not provided in INI)
    n_samples = 1000
    max_iterations = 10000
    mmin = float(sec.get('mmin', 0.5))
    mmax = float(sec.get('mmax', 350.0))
    verbose = sec.get('verbose', 'True').lower() in ('1', 'true', 'yes', 'on')

    print(f"Sampling {n_samples} binary masses from NotchFilterBinnedPairingMassDistribution...")
    m1, m2, acc_rate = sample_notch_filter_binned_pairing_masses(
        n_samples=n_samples,
        max_iterations=max_iterations,
        verbose=verbose,
        **hyperparams,
    )

    # Plot the samples if any
    if len(m1) > 0:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))

        # m1 distribution
        axes[0].hist(m1, bins=30, alpha=0.7, edgecolor='black')
        axes[0].set_xlabel('Primary Mass (M_sun)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Primary Mass Distribution')

        # m2 distribution
        axes[1].hist(m2, bins=30, alpha=0.7, edgecolor='black')
        axes[1].set_xlabel('Secondary Mass (M_sun)')
        axes[1].set_ylabel('Count')
        axes[1].set_title('Secondary Mass Distribution')

        # m1 vs m2 scatter
        axes[2].scatter(m1, m2, alpha=0.5, s=10)
        axes[2].plot([hyperparams['NSmin'], mmax],
                     [hyperparams['NSmin'], mmax], 'r--', label='m1 = m2')
        axes[2].set_xlabel('Primary Mass (M_sun)')
        axes[2].set_ylabel('Secondary Mass (M_sun)')
        axes[2].set_title('Mass Ratio Distribution')
        axes[2].legend()
        axes[2].set_aspect('equal')

        plt.tight_layout()
        plt.savefig('mass_samples_visualization.png', dpi=150, bbox_inches='tight')
        print("Saved visualization to 'mass_samples_visualization.png'")
    else:
        print("No samples were generated.")
