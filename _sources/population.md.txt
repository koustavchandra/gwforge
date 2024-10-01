# Simulating the source population

GWForge can simulate a wide range of **binary** source populations (at the moment). This page tells you how you can use them.

The first step in generating a source population is to determine the distribution of sources in luminosity distance $\(D_L\)$ (or equivalently redshift $\(z\))$ and the expected number of signals in the data. So, we will start by setting up the `[Redshift]` section.

## Redshift
 For this, you need to specify:

- Redshift distribution model
- Local merger rate density in $\(Gpc^{-3}yr^{-1}\)$
- Maximum redshift of the source
- Cosmological parameters such as $H_0,~O_{m0}, O_{de0}$ and $T_{cmb0}$, assuming [LambdaCDM cosmology](http://hyperphysics.phy-astr.gsu.edu/hbase/Astro/lambda.html)
- A reference start time when you switch on your detector.

The last is optional. If not provided, GWForge assumes [Planck18](https://docs.astropy.org/en/latest/api/astropy.cosmology.realizations.Planck18.html) cosmology. 
$$
H_0 = 67.66~km/s/Mpc,~O_{m0} = 0.30966,~T_{cmb0} = 2.7255 K
$$

Alternatively, you can select any of the [cosmological realisations implemented in astropy](https://docs.astropy.org/en/stable/cosmology/realizations.html).

If you choose to define a *custom* Universe, you can do so as follows:
```ini
[Redshift]
redshift-model = MadauDickinson
redshift-parameters = {'gamma': 2.7, 'kappa': 5.6, 'z_peak': 1.9}
local-merger-rate-density = 22
maximum-redshift = 30
; custom cosmology
cosmology = custom 
H0 = 70
Om0 = 0.3
Ode0 = 0.7
Tcmb0 = 2.735
; analysis start time
gps-start-time = 1893024018
```


For reference, the Madau-Dickinson distribution function is:
$$
p(z | \gamma, \kappa, z_\mathrm{peak}) \propto \frac{1}{1+z} \frac{dV_c}{dz} \psi(z | \gamma, \kappa, z_\mathrm{peak})
$$
where
$$
\psi(z | \gamma, \kappa, z_\mathrm{peak}) \propto \frac{(1+z)^\gamma}{1+\big(\frac{1+z}{1+z_\mathrm{peak}}\big)^\kappa} M_\odot/\mathrm{year}/\mathrm{Mpc}^3
$$
and
| parameter | descsription| 
| --- | ---|
| $\gamma$ | Slope of the distribution at low redshift |
| $\kappa$ | Slope of the distribution at high redshift |
| $z_\textrm{peak}$ | Redshift at which the distribution peaks. |
| $z_\textrm{max}$ | The maximum redshift allowed. |

It is important to note that this describes the progenitor formation rate distribution. To obtain the compact binary merger rate distribution, GWForge convolves it with a time-delay distribution.

 The reader can refer to the following for further details:
* [A Mock Data Challenge for the Einstein Gravitational-Wave Telescope](https://inspirehep.net/literature/1084847)
* [Mock data study for next-generation ground-based detectors: The performance loss of matched filtering due to correlated confusion noise](https://inspirehep.net/literature/2148213)


```{note}
The above implementation assumes that all compact binary systems are formed by [isolated binary evolution via the common-envelope phase](https://www.frontiersin.org/articles/10.3389/fspas.2020.00038/full). It also adopts a flat-in-log distribution for a time delay between binary formation and merger.
```

## Mass
The [Mass] section helps define the mass distribution of the binary population. Similar to the [Redshift] section, a model name and a dictionary of parameters must be provided. For example:

```ini
[Mass]
mass-model = PowerLaw+Peak
mass-parameters = {'alpha':3.37, 'beta': 0.76, 'delta_m':5.23,  'mmin':4.89, 'mmax':88.81, 'lam':0.04, 'mpp': 33.60, 'sigpp':4.59}
```


The currently available mass distribution models and their parameters are:
<details>
  <summary>List of mass distribution models</summary>

  | Model Name | Parameters | Description | 
  | ---|---| ---| 
  |[`PowerLaw+Peak`](https://colmtalbot.github.io/gwpopulation/_autosummary/gwpopulation.models.mass.SinglePeakSmoothedMassDistribution.html#gwpopulation.models.mass.SinglePeakSmoothedMassDistribution)| `alpha, beta, mmin, mmax, lam, mpp, sigpp, delta_m` | Powerlaw + peak model for two-dimensional mass distribution with low mass smoothing.
  |[`MultiPeak`](https://colmtalbot.github.io/gwpopulation/_autosummary/gwpopulation.models.mass.MultiPeakSmoothedMassDistribution.html#gwpopulation.models.mass.MultiPeakSmoothedMassDistribution)| `alpha, beta, mmin, mmax, lam, lam_1, mpp_1, mpp_2, sigpp_1, sigp_2, delta_m` | Powerlaw + two peak model for two-dimensional mass distribution with low mass smoothing.
  |[`BrokenPowerLaw`](https://colmtalbot.github.io/gwpopulation/_autosummary/gwpopulation.models.mass.BrokenPowerLawSmoothedMassDistribution.html#gwpopulation.models.mass.BrokenPowerLawSmoothedMassDistribution)| `alpha_1, alpha_2, beta, break_fraction, mmin, mmax, delta_m` | Broken power law for two-dimensional mass distribution with low mass smoothing. |
  |`UniformSecondary`| `alpha, beta, delta_m, mmin, mmax, 88.81, lam, mpp, sigpp, minimum_secondary_mass, maximum_secondary_mass` | PowerLaw + Peak for primary mass and uniform for secondary |
  |`DoubleGaussian`| `mu_1, sigma_1, mu_2, sigma_2, breaking_fraction, mmin, mmax` | Truncated Gaussian distribution for primary and secondary
  |`LogNormal`| `mu, sigma` | Log-normal distribution with mean mu and width sigma for primary and secondary | 
  |`PowerLawDipBreak`|`mmin, mmax, alpha_1, alpha_2, gamma_low, gamma_high, eta_low, eta_high, A, n` | Extension of power law break model |
  |`PowerLaw`| `alpha_1, mmin, mmax`  | Power law with bounds and alpha, spectral index for primary and secondary |      

The parameter names are heavily dependent on gwpopulation and bilby. Thus, it is essential to keep track of definition changes.

For more details, refer to the following publications:
* [Binary Black Hole Population Properties Inferred from the First and Second Observing Runs of Advanced LIGO and Advanced Virgo](https://inspirehep.net/literature/1706043)
* [Population Properties of Compact Objects from the Second LIGO-Virgo Gravitational-Wave Transient Catalog](https://inspirehep.net/literature/1826636)
* [Population of Merging Compact Binaries Inferred Using Gravitational Waves through GWTC-3](https://inspirehep.net/literature/1961598)

</details>

```{note}
GWForge overlooks special characters and converts everything to lower cases. So `PowerLaw+Peak` is equivalent to `powerlawpeak`.
```

## Spin
The `[Spin]` section determines the spin distribution of the population. For example:
```ini
[Spin]
spin-model = Beta-Aligned
spin-parameters = {'minimum_primary_spin' : 0, 'maximum_primary_spin':  0.99, 'minimum_secondary_spin' : 0, 'maximum_secondary_spin' : 0.5, 'mu_chi' : 0.26, 'sigma_squared_chi' : 0.02}
```
defines a quasi-circular (non-precessing) binary population whose spin magnitude is sampled from a beta distribution.

Here is the list of currently available spin distribution

<details>
  <summary>List of spin distribution models</summary>

  |Model | Parameters | Description|
  |---|---|---|
  |`Non-spinning`| `None` | Non-spinning 
  |`Aligned`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | [Aligned spin distribution Bilby-style](https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.prior.AlignedSpin.html#bilby.gw.prior.AlignedSpin)| 
  |`Aligned-Bilby`|`minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | [Aligned spin distribution Bilby-style](https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.prior.AlignedSpin.html#bilby.gw.prior.AlignedSpin)|
  |`Aligned-Uniform`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | Aligned component of spins are sampled from uniform distribution|
  |`Beta-Aligned`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi`| Bilby style aligned spin distribution Bilby-style with spin magnitudes obeying Beta distribution |
  |`Aligned-Gaussian-Uniform`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin,mu_chi_1, sigma_chi_1` | Aligned component of primary is sampled from Truncated Gaussian and secondary from uniform |
  |`Isotropic`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | Spin Magnitudes sampled from Uniform distribution + Isotropic distribution of spin angles |
  |`Isotropic-Beta`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi` | Spin Magnitudes sampled from Beta distribution. Isotropic distribution of spin angles
  |`Isotropic-Beta_Gaussian`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi, sigma_t` | Spin magnitudes sampled from Beta distribution. Truncated Gaussian distribution for cosine tilt angles.
  |`Isotropic-Beta_Gaussian_Uniform`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi, sigma_t, xi_spin` | Spin magnitudes sampled from Beta distribution. A fraction of the binaries have cosine tilt angles from Truncated Gaussian distribution and the rest from a uniform distribution between (-1,1)|
  | `Default` | `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi, sigma_t, xi_spin` | Same as `Isotropic-Beta_Gaussian_Uniform` |

For more details, refer to the following publications:
* [Binary Black Hole Population Properties Inferred from the First and Second Observing Runs of Advanced LIGO and Advanced Virgo](https://inspirehep.net/literature/1706043)
* [Population Properties of Compact Objects from the Second LIGO-Virgo Gravitational-Wave Transient Catalog](https://inspirehep.net/literature/1826636)
* [Population of Merging Compact Binaries Inferred Using Gravitational Waves through GWTC-3](https://inspirehep.net/literature/1961598)

</details>

## Extrinsic
The [Extrinsic] section is designed to handle sky location and binary orientation parameters. You can specify a bilby prior file as input. By default, it assumes an isotropic distribution for sky location and orientation parameters and a uniform distribution for the polarization angle.

For example:
```ini
[Extrinsic]
```
will use the second.

## EOS
The `[EOS]` section accepts an `eos-file` from which it reads the mass and tidal parameters. By default it uses SLy EOS but you can specify as follows:
```ini
[EOS]
eos-file = /ligo/home/ligo.org/koustav.chandra/projects/Cosmic-Explorer-MDC/gwforge/GWForge/inject/eos_tables/TOVSeq_SLy.dat
```
provided it is consistent with how Rahul likes to define them.

### Generating the population.

To generate the binary parameters for the population, execute the following:
```bash
gwforge_population --config-file bbh.ini --output-file bbh.h5
```
It should take at most a minute to generate the output file. By default `gwforge_population` assumes your source type is BBH. For other options, please check `gwforge_population --help`. Please note that the waveform approximant that you use for your waveform generation supports tidal parameters if the source-type is bns or nsbh.

```note
By default a year gwforge_population generates a year worth of population. If you want some other value, please add the `duration` flag and add a value in seconds.
Example: `duration=4096`. Please note that the population generated should be greater than the number of signals injected.
```

A few more example configuration file exist here: `~/.conda/envs/gwforge-venv/lib/python3.9/site-packages/GWForge/population/`. Feel free to modify and see what you get.

### Naive way to check the population
You can check the binary parameters of the population by doing the following:
```python
from GWForge.utils import cornerplot
cornerplot(file='bbh.h5', parameters=['mass_1_source', 'mass_2_source', 'spin_1z','spin_2z',  'redshift'], save='pop.png')
```
This will create a plot called `pop.png` in the current working directory with the parameters. The list of parameters can be found by doing `h5ls -r bbh.h5`. It list all the keys of an HDF5 file.