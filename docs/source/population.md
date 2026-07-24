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

It is important to note that $\psi(z)$ describes the progenitor *formation* rate distribution. To obtain the compact binary *merger* rate distribution, GWForge convolves it with a time-delay distribution $p(\tau)$ between the formation of the binary and its eventual merger:
$$
\mathcal{R}(z_m) \propto \int_{z_m}^{\infty} \psi(z_f)\, p\big(t(z_f) - t(z_m)\big)\, \frac{dt}{dz_f}\, dz_f,
$$
where $t(z)$ is the lookback time and the lower limit $z_f \geq z_m$ enforces that a binary must form before it merges.

By default GWForge uses the `inverse` time-delay distribution, $p(\tau) \propto 1/\tau$ (i.e. flat-in-log), between a minimum delay of $20~\mathrm{Myr}$ and the age of the Universe. This is the canonical choice for compact binaries formed through [isolated binary evolution via the common-envelope phase](https://www.frontiersin.org/articles/10.3389/fspas.2020.00038/full). If you want a different slope, use the `powerlaw` model, $p(\tau) \propto \tau^{-\mathrm{slope}}$, and set `time-delay-slope`:

```ini
[Redshift]
redshift-model = MadauDickinson
redshift-parameters = {'gamma': 2.7, 'kappa': 5.6, 'z_peak': 1.9}
local-merger-rate-density = 22
maximum-redshift = 10
gps-start-time = 1893024018
; power-law time delay instead of the default 1/tau
time-delay-model = powerlaw
time-delay-slope = 0.8
```

| parameter | description |
| --- | --- |
| `time-delay-model` | `inverse` (default, $\propto 1/\tau$) or `powerlaw` ($\propto \tau^{-\mathrm{slope}}$) |
| `time-delay-slope` | Slope of the power-law time delay. Ignored (fixed to 1) for `inverse` |

If, instead of Madau-Dickinson, you would rather scale the merger rate as a simple power law of redshift, $\psi(z) = (1+z)^{\lambda}$, use the `PowerLaw` model:
```ini
[Redshift]
redshift-model = PowerLaw
redshift-parameters = {'lamb': 2.7}
local-merger-rate-density = 22
maximum-redshift = 10
gps-start-time = 1893024018
```

 The reader can refer to the following for further details:
* [A Mock Data Challenge for the Einstein Gravitational-Wave Telescope](https://inspirehep.net/literature/1084847)
* [Mock data study for next-generation ground-based detectors: The performance loss of matched filtering due to correlated confusion noise](https://inspirehep.net/literature/2148213)


```{note}
The redshift, time-delay, mass and spin models are all implemented in-house on top of `bilby` — GWForge no longer depends on `gwpopulation` or `pycbc` to build the population. The numbers are validated against those original implementations, so nothing changes for you except a lighter set of dependencies.
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
  |`BGP`| `alpha_1, alpha_2, m_break, mmin, m_high, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, delta_m, beta` | Broken power law + two Gaussian peaks — the fiducial BBH mass model from GWTC-4.0/5.0. See [below](#the-bgp-model). |
  |`UniformSecondary`| `alpha, beta, delta_m, mmin, mmax, 88.81, lam, mpp, sigpp, minimum_secondary_mass, maximum_secondary_mass` | PowerLaw + Peak for primary mass and uniform for secondary |
  |`DoubleGaussian`| `mu_1, sigma_1, mu_2, sigma_2, breaking_fraction, mmin, mmax` | Truncated Gaussian distribution for primary and secondary
  |`LogNormal`| `mu, sigma` | Log-normal distribution with mean mu and width sigma for primary and secondary | 
  |`PowerLawDipBreak`|`mmin, mmax, alpha_1, alpha_2, gamma_low, gamma_high, eta_low, eta_high, A, n` | Extension of power law break model |
  |`PowerLaw`| `alpha, mmin, mmax`  | Power law with bounds and alpha, spectral index for primary and secondary |
  |`Uniform_components`| `mmin, mmax` | Both component masses drawn uniformly in `[mmin, mmax]` and ordered so that $m_1 \geq m_2$ |
  |`Uniform_M_q`| `minimum_total_mass, maximum_total_mass, minimum_mass_ratio, maximum_mass_ratio` | Uniform in total mass and mass ratio |
  |`FullPop_GWTC4`| `A, A2, NSmin, NSmax, BHmin, BHmax, UPPERmin, UPPERmax, n0..n5, alpha_1, alpha_2, alpha_dip, mu1, sig1, mix1, mu2, sig2, mix2, beta_pair_1, beta_pair_2, mbreak, mmin, mmax` | Full "Power Law + Dip + Break" compact-binary mass function with a pairing function (GWTC-4.0). See [below](#the-fullpop_gwtc4-model). |
  |`UserDefined`| `file, primary_parameter, (mass_ratio_parameter \| secondary_parameter)` | Tabulated distributions supplied by you, e.g. from a population-synthesis run. See [below](#user-defined-populations). |

The parameter definitions follow `gwpopulation` and the LVK population papers, even though the models are now implemented in-house. It is therefore still worth keeping an eye on how those references define things.

For more details, refer to the following publications:
* [Binary Black Hole Population Properties Inferred from the First and Second Observing Runs of Advanced LIGO and Advanced Virgo](https://inspirehep.net/literature/1706043)
* [Population Properties of Compact Objects from the Second LIGO-Virgo Gravitational-Wave Transient Catalog](https://inspirehep.net/literature/1826636)
* [Population of Merging Compact Binaries Inferred Using Gravitational Waves through GWTC-3](https://inspirehep.net/literature/1961598)

</details>

```{note}
GWForge overlooks special characters and converts everything to lower cases. So `PowerLaw+Peak` is equivalent to `powerlawpeak`.
```

### The BGP model
`BGP` (Broken power law + Gaussian Peaks) is the fiducial BBH mass model used in the GWTC-4.0 and GWTC-5.0 population analyses. The primary mass is a mixture of a broken power law and two left-truncated Gaussian peaks, with a low-mass Planck taper $S$ applied to the whole distribution:
$$
\pi(m_1) \propto \Big[\lambda_0\, p_\mathrm{BP}(m_1) + \lambda_1\, N_\mathrm{lt}(m_1 | \mu_1, \sigma_1) + (1 - \lambda_0 - \lambda_1)\, N_\mathrm{lt}(m_1 | \mu_2, \sigma_2)\Big]\, S(m_1 | m_\mathrm{min}, \delta_m),
$$
where the broken power law switches from slope $-\alpha_1$ to $-\alpha_2$ at $m_\mathrm{break}$ over $[m_\mathrm{min}, m_\mathrm{high}]$, and $N_\mathrm{lt}$ is a normal distribution truncated below at $m_\mathrm{min}$. The mass ratio is a power law $\propto q^{\beta}$ with the same low-mass taper applied to the secondary. For example:

```ini
[Mass]
mass-model = BGP
mass-parameters = {'alpha_1': 1.6, 'alpha_2': 5.0, 'm_break': 38.0, 'mmin': 5.0, 'm_high': 100.0, 'lam_0': 0.9, 'lam_1': 0.05, 'mpp_1': 33.0, 'sigpp_1': 4.0, 'mpp_2': 10.0, 'sigpp_2': 1.5, 'delta_m': 4.8, 'beta': 1.1}
```

| parameter | description |
| --- | --- |
| `alpha_1, alpha_2` | Power-law slopes below and above the break |
| `m_break` | Mass at which the power law breaks |
| `mmin, m_high` | Lower and upper edges of the power-law component |
| `lam_0, lam_1` | Mixing fractions of the power law and the first peak (the second peak gets $1 - \lambda_0 - \lambda_1$) |
| `mpp_1, sigpp_1` | Location and width of the first (high-mass) Gaussian peak |
| `mpp_2, sigpp_2` | Location and width of the second (low-mass) Gaussian peak |
| `delta_m` | Width of the low-mass smoothing |
| `beta` | Power-law slope of the mass ratio |
| `maximum_mass` | *(optional)* upper bound of the evaluation grid, default `200` $M_\odot$ |

The exact definitions are Eqs. (B10)–(B14) of the [GWTC-5.0 population paper](https://arxiv.org/abs/2605.27226).

### User-defined populations
Sometimes you already have a population — say, the output of a population-synthesis code — and you simply want GWForge to draw from it. The `UserDefined` model lets you do exactly that: you hand it a JSON file that tabulates the support `xx` and the (un-normalised) probability density `yy` of each parameter, and GWForge turns each into a prior it can sample from. The file looks like this:

```json
{
  "mass_1_source": {"xx": [...], "yy": [...]},
  "mass_ratio":    {"xx": [...], "yy": [...]}
}
```

You then point the `[Mass]` section at the file and tell it which key is the primary mass, and either a mass-ratio key or a secondary-mass key:

```ini
[Mass]
mass-model = UserDefined
mass-parameters = {'file': '/path/to/priors.json', 'primary_parameter': 'mass_1_source', 'mass_ratio_parameter': 'mass_ratio'}
```

| parameter | description |
| --- | --- |
| `file` | Path to the JSON population-priors file |
| `primary_parameter` | Key to use for the primary mass (default `mass_1_source`) |
| `mass_ratio_parameter` | Key to use for the mass ratio |
| `secondary_parameter` | Key to use for the secondary mass (use this *instead of* `mass_ratio_parameter`) |

```{note}
If you provide the primary and secondary masses as two independent distributions, GWForge samples them independently and then orders each pair so that $m_1 \geq m_2$; the pairing of the original population is therefore not preserved. If you care about the pairing, provide the mass ratio instead.
```

### The FullPop_GWTC4 model
`FullPop_GWTC4` is the strongly-parameterised "Power Law + Dip + Break" model that describes the *whole* compact-binary mass spectrum — from neutron stars, through the lower mass gap, to black holes — together with a pairing function that sets how strongly binaries favour equal masses. It carries a large number of hyperparameters (the notch edges, filter sharpnesses `n0..n5`, peak locations and pairing slopes listed in the table above). Because there is no closed-form inverse, the masses are drawn with an importance sampler; the default is `importance_m1_m2`. The model is described in App. B.2 of the [GWTC-5.0 population paper](https://arxiv.org/abs/2605.27226) (following Farah et al. 2022 and Mali & Essick 2025).

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
  |`Gaussian-Non_spinning`| `mu_chi_1, sigma_chi_1, minimum_primary_spin, maximum_primary_spin` | Primary aligned spin drawn from a Truncated Gaussian; secondary non-spinning |
  |`Aligned`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | [Aligned spin distribution Bilby-style](https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.prior.AlignedSpin.html#bilby.gw.prior.AlignedSpin)| 
  |`Aligned-Bilby`|`minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | [Aligned spin distribution Bilby-style](https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.prior.AlignedSpin.html#bilby.gw.prior.AlignedSpin)|
  |`Aligned-Uniform`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | Aligned component of spins are sampled from uniform distribution|
  |`Beta-Aligned`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin, mu_chi, sigma_squared_chi`| Bilby style aligned spin distribution Bilby-style with spin magnitudes obeying Beta distribution |
  |`Aligned-Gaussian-Uniform`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin,mu_chi_1, sigma_chi_1` | Aligned component of primary is sampled from Truncated Gaussian and secondary from uniform |
  |`Isotropic-Bilby`| `minimum_primary_spin, minimum_secondary_spin, maximum_primary_spin, maximum_secondary_spin` | Spin Magnitudes sampled from Uniform distribution + Isotropic distribution of spin angles |
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
The `[Extrinsic]` section handles the sky location and binary orientation. By default GWForge assumes an isotropic sky (right ascension uniform, declination cosine-distributed), isotropic orientation (inclination sine-distributed), and a uniform polarization angle. An empty section is enough to get these defaults:
```ini
[Extrinsic]
```

If you want something other than the defaults, you can point to your own [bilby prior file](https://lscsoft.docs.ligo.org/bilby/prior.html):
```ini
[Extrinsic]
extrinsic-prior-file = /path/to/my_extrinsic.prior
```

You can also override just the inclination without writing a full prior file. Setting `inclination-distribution = schutz` draws inclinations from the [Schutz (2011)](https://arxiv.org/abs/1102.5421) distribution, which weights orientations by their detectability — handy if you want a detection-like sample rather than a strictly isotropic one:
```ini
[Extrinsic]
inclination-distribution = schutz
```

| parameter | description |
| --- | --- |
| `extrinsic-prior-file` | *(optional)* bilby prior file overriding the default isotropic priors |
| `inclination-distribution` | *(optional)* set to `schutz` for the detectability-weighted inclination distribution |

## EOS (Equation of State)

The `[EOS]` section allows specifying an `eos-file` that provides the mass and tidal parameters for a neutron star equation of state. By default, the SLy EOS (Skyrme-Lyon) is used, but you can override this by specifying a different `eos-file` in the following way:

```ini
[EOS]
eos-file = /ligo/home/ligo.org/koustav.chandra/projects/Cosmic-Explorer-MDC/gwforge/GWForge/inject/eos_tables/TOVSeq_SLy.dat
```
provided it is consistent with how Rahul likes to define them. You can find examples of `eos-tables` in the [GWForge repository](https://github.com/koustavchandra/gwforge/tree/main/GWForge/inject/eos_tables).


**Structure of eos-file**
The EOS tables are structured in columns, where each column corresponds to different physical parameters. Below is a breakdown of few of these columns:

`C`: Compactness of neutron star.

`Mb`: Baryonic mass in solar masses

`M`: Mass in solar mass units

`R`: Radius in solar mass units

`kl` : Second love number


### Generating the population.

To generate the binary parameters for the population, execute the following:
```bash
gwforge_population --config-file bbh.ini --output-file bbh.h5
```
It should take at most a minute to generate the output file. By default `gwforge_population` assumes your source type is BBH. For other options, please check `gwforge_population --help`. Please note that the waveform approximant that you use for your waveform generation supports tidal parameters if the source-type is bns or nsbh.

If you want the population to be reproducible, pass a `--seed`:
```bash
gwforge_population --config-file bbh.ini --output-file bbh.h5 --seed 42
```
Running with the same seed and the same configuration gives you byte-for-byte the same population every time. Without a seed, each run gives you a fresh realisation.

```{note}
By default `gwforge_population` generates a year's worth of population. If you want some other value, add the `duration` option (in seconds) to the `[Redshift]` section — for example `duration = 4096`. Please note that the population generated should be greater than the number of signals injected.
```

A few more example configuration files ship with the package under `GWForge/population/population_configuration_files/` (inside your environment's `site-packages`, or in the [source tree](https://github.com/koustavchandra/gwforge/tree/main/GWForge/population/population_configuration_files)). Feel free to modify them and see what you get.

### Naive way to check the population
You can check the binary parameters of the population by doing the following:
```python
from GWForge.utils import cornerplot
cornerplot(file='bbh.h5', parameters=['mass_1_source', 'mass_2_source', 'spin_1z','spin_2z',  'redshift'], save='pop.png')
```
This will create a plot called `pop.png` in the current working directory with the parameters. The list of parameters can be found by doing `h5ls -r bbh.h5`. It list all the keys of an HDF5 file.
