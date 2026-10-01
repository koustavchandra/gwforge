# Population Fisher Forecasts

`gwforge_fisher` asks how well a network measures **one source**.
`gwforge_population_fisher` asks how well it measures **the population**: the
shape of the mass function, the history of the merger rate, the spin
distribution, and — through spectral sirens — the cosmology.

The estimator is the first term of the hyper-parameter Fisher expansion of
[Gair, Antonelli & Barbieri (2022)](https://arxiv.org/abs/2205.07893), Eq. 21. For
$N_{\rm det}$ detected events $\theta_k$ and hyper-parameters $\Lambda$,

$$s_k^i = \frac{\partial \ln p(\theta_k \mid \Lambda)}{\partial \Lambda^i},
\qquad
\langle s\rangle^i = \frac{1}{N_{\rm det}}\sum_k s_k^i,
\qquad
\Gamma^{ij} = \sum_k \big(s_k^i - \langle s\rangle^i\big)
                     \big(s_k^j - \langle s\rangle^j\big),$$

and $\Sigma = \Gamma^{-1}$.

## Where the selection function went

Not away — it is never *evaluated*, which is a different thing. Term I is an
expectation under the **detected** density,

$$\Gamma_1^{ab} = -\int \frac{\partial^2}{\partial\Lambda^a \partial\Lambda^b}
  \ln\frac{p(\theta|\Lambda)}{\alpha(\Lambda)}\;
  p_{\rm det}(\theta|\Lambda)\,\mathrm{d}\theta,$$

so the selection appears twice: as the measure, and through the detectable
fraction $\alpha(\Lambda)$. Two facts collapse both into the centred sum.
First, $p_{\rm det}(\theta)$ carries no hyper-parameter, so it drops from the
score. Second, the detected catalogue is itself a fair Monte-Carlo draw from
$p_{\rm det}(\theta)\,p(\theta|\Lambda)/\alpha(\Lambda)$, so

$$\frac{\partial \ln \alpha}{\partial \Lambda^i}
  = \int p_{\rm det}(\theta)\frac{p(\theta|\Lambda)}{\alpha(\Lambda)}
    \frac{\partial \ln p}{\partial \Lambda^i}\,\mathrm{d}\theta
  = \mathbb{E}_{\rm det}\!\left[s^i\right],$$

which the sample mean estimates. So there is no separate detection integral to
evaluate, no injection reweighting, and no rate normalisation — but the
selection is still doing work, carried entirely by **which events are in the
sum**. The SNR threshold is an input to the forecast, not a preprocessing step:
change it and $\Gamma$ changes. The scale of the forecast is set by how many
events were detected, so every $\sigma$ scales as $1/\sqrt{T_{\rm obs}}$.

`validation/population_detection_fraction.py` reports what a given threshold
removes, overall and binned in redshift and primary mass:

```bash
python validation/population_detection_fraction.py \
    --population bbh_population.h5 --snr bbh_snr.h5 \
    --detectors CE40 CE20 ET --snr-thresholds 10 20
```

$p_{\rm det}$ reappears *explicitly* at higher order — as
$\partial_\theta \ln p_{\rm det}$ in Gair's $\Gamma_4$ — which is one reason
this package stops at Term I.

```{note}
The identity is exact only if $p_{\rm det}$ carries no hyper-parameter. That is
why the spectral-siren model builds its density in the **detector frame** rather
than pushing events to the source frame at a trial cosmology, as described under
Spectral sirens below.
```

## Quick start

Three commands: generate a population, compute its SNRs, forecast.

```bash
gwforge_population --config-file GWForge/population/population_configuration_files/bgp-gwtc5.ini \
                   --output-file bbh_population.h5 --source-type bbh \
                   --reference-frequency 5 --seed 250114 --save-config

gwforge_optimal_snr --injection-file bbh_population.h5 --output-file bbh_snr.h5 \
                    --ifos CE40 CE20 ET --waveform-approximant IMRPhenomXPHM \
                    --minimum-frequency 5 --sampling-frequency 2048 --cores 10

gwforge_population_fisher --config-file mass_redshift_spin.ini
```

```ini
[Catalogue]
population-file = bbh_population.h5
snr-file = bbh_snr.h5
detectors = ['CE40', 'CE20', 'ET']
snr-thresholds = [10, 20]

[Model]
blocks = ['mass', 'redshift', 'spin']
; Must match the generation config: a different support is a different model.
maximum-redshift = 10
cosmology = Planck18
; Omit these five and they fall back to the GWTC-5.0 medians in
; GWForge.population.mass.BGP_PARAMETERS, which is also what
; bgp-gwtc5.ini generates.
mmin = 4.489078237
m-high = 300.0
maximum-mass = 300.0
mmin-2 = 3.487749145
delta-m-2 = 5.60056759
mass-parameters = {'alpha_1': 1.456442737, 'alpha_2': 5.100400428, 'm_break': 37.9806382, 'lam_0': 0.4206178692, 'lam_1': 0.5228817702, 'mpp_1': 9.989384122, 'sigpp_1': 0.6601839103, 'mpp_2': 33.2656028, 'sigpp_2': 4.582221363, 'delta_m': 3.123416382, 'beta': 0.8049660633}
redshift-parameters = {'gamma': 2.7, 'kappa': 5.6, 'z_peak': 1.9}
spin-parameters = {'mu_chi': 0.0751318233, 'sigma_chi': 0.3667469243, 'mu_t': 0.2788662234, 'sigma_t': 0.9661688072, 'xi_spin': 0.6650168221}

[Fit]
fit-redshift = True
fit-mass = True
fit-spin = True

[Fisher]
; delta_m is held at its generation value; see "Which mass parameters are free".
free-parameters = ['alpha_1', 'alpha_2', 'm_break', 'lam_0', 'lam_1', 'mpp_1', 'sigpp_1', 'mpp_2', 'sigpp_2', 'beta', 'gamma', 'kappa', 'z_peak', 'mu_chi', 'sigma_chi', 'mu_t', 'sigma_t', 'xi_spin']
score-method = analytic

[Output]
output-directory = population_fisher_output
label = mass_redshift_spin
plot-corner = True
plot-parameters = ['beta', 'gamma', 'kappa', 'z_peak', 'mu_chi']
```

Ready-made configs ship in `GWForge/population_fisher/configuration_files/`:
`redshift.ini`, `mass.ini`, `mass_redshift.ini`, `mass_redshift_spin.ini`,
`spectral_siren.ini`, and two that point the same analyses at 2G networks,
`mass_lvk.ini` (H1 + L1 + V1) and `spectral_siren_aplus.ini` (H1 + L1 + A1 at
A+). The catalogue they all read comes from `bgp-gwtc5.ini`, which lives with the
other population configs in `GWForge/population/population_configuration_files/`.

## The models

The densities are **GWForge's own, called directly rather than reimplemented**,
so a forecast and the catalogue it runs on describe the same population by
construction.

| block | hyper-parameters | density |
| --- | --- | --- |
| `mass` | `alpha_1, alpha_2, m_break, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, delta_m, beta` | {mod}`GWForge.population._smoothed_mass` — the BGP model, GWTC-5.0 Eqs. B10–B14 |
| `redshift` | `gamma, kappa, z_peak` | {func}`GWForge.population.redshift.madau_dickinson_psi_of_z` |
| `spin` | `mu_chi, sigma_chi, mu_t, sigma_t, xi_spin` | {func}`GWForge.population.spin.default_spin_magnitude_density` and {func}`GWForge.population.spin.default_spin_tilt_density`, both built on {func}`GWForge.population.spin.truncated_normal` |
| `cosmology` | `H0, Om0, w0` | {class}`GWForge.cosmology.FlatwCDM` |

The names are the ones the generation config uses, so the two files read the
same. Every fiducial above may be omitted, in which case it falls back to the
GWTC-5.0 median in {data}`GWForge.population.mass.BGP_PARAMETERS` or
{data}`GWForge.population.spin.DEFAULT_BBH_SPIN_PARAMETERS` — the single place
those numbers are written down.

```{warning}
`sigma_chi` is a **standard deviation**. The Beta spin models take a variance
under the name `sigma_squared_chi`, and a config carrying the old name is
rejected rather than silently reinterpreted.
```

### Which mass parameters are free, and why the rest are not

`mmin`, `m_high` and `maximum_mass` are **constructor arguments, not
hyper-parameters**, for two independent reasons.

*Statistically*, they are hard cutoffs. Inside the power-law interior their
scores are the same number for every event, so the centred score
$s_k - \langle s\rangle$ vanishes identically and the Fisher acquires a zero
row. Within Term I such a parameter is simply unmeasurable.

*Numerically*, they are the endpoints of the normalisation integrals, so their
derivatives carry Leibniz boundary terms and cannot be checked against a finite
difference on a fixed grid — as `m_high` crosses a node, the trapezoid jumps by
the integrand times the node spacing, so the finite difference is a staircase
rather than a derivative.

To study sensitivity to them, build a second model with different values. A
different support is a different model, not a different point in one.

`delta_m` is smooth — the Planck taper drives the integrand continuously to
zero at $m_{\min}$, so there is no edge for a boundary term to sit on — and it
*can* be freed. Every shipped config nevertheless holds it at its generation
value: it is the width of a sharp feature at a fixed source-frame mass, so it is
a standard scale, and it is the mass parameter that feeds $H_0$ most directly.
Freeing it lets the taper absorb the cosmology. Scan it by running the forecast
at several values instead.

## Derivatives

Two tiers, in strict order of preference.

| tier | used for | notes |
| --- | --- | --- |
| **analytic** | every model GWForge ships, cosmology included | the production path |
| **finite difference** | the test oracle, and a model added later whose score nobody has written out | never a production path |

There is deliberately no autodiff tier. Every model here has closed-form scores,
so JAX would have been a dependency the production path never executed.

Tier one covers everything in scope, **including the cosmology columns**. The
comoving distance is a quadrature,

$$D_C(z) = D_H \int_0^z \frac{\mathrm{d}z'}{E(z')},$$

so differentiating the *result* means finite-differencing an integral. But
differentiating *under* the integral leaves a closed-form integrand on the very
same nodes,

$$\frac{\partial}{\partial\Omega_{m,0}}\frac{1}{E}
   = -\frac{(1+z)^3 - (1+z)^{3(1+w_0)}}{2E^3},
\qquad
\frac{\partial}{\partial w_0}\frac{1}{E}
   = -\frac{3(1-\Omega_{m,0})(1+z)^{3(1+w_0)}\ln(1+z)}{2E^3},$$

which costs one extra weighted sum and is exact to the quadrature's accuracy.

```{note}
This is not a stylistic preference. Finite-differencing a quantity that is
itself a quadrature has an error floor set by the quadrature, not by the step,
so the error can *grow* as the step shrinks. Every derivative here is checked
against a finite difference in `tests/test_cosmology.py` and
`tests/test_population_fisher.py`, and the finite difference converges onto the
analytic value at the expected $O(h^2)$ — it is the approximation, not the
reference.
```

## Fiducials: the one place a fit is needed

`gwforge_population` samples masses and spins **directly from the same densities
this module evaluates**, so for those the fiducials are exactly the config
values.

Redshift is different. {class}`GWForge.population.redshift.Redshift` always
convolves $\psi_{\rm MD}$ with a formation-to-merger time-delay distribution, so
a generated catalogue is *not* distributed as
$\psi_{\rm MD}(z)\,\mathrm{d}V_c/\mathrm{d}z/(1+z)$. Running a forecast for the
direct Madau–Dickinson form at the generation $(\gamma, \kappa, z_p)$ would be a
forecast about the wrong point. `fit-redshift = True` therefore fits the direct
form to the **injected** redshifts — all of them, not just the detected ones,
because the population model describes the astrophysical population and fitting
to the detected subset would absorb the selection function into the
hyper-parameters.

With the `inverse` time delay the fit lands well away from the generation
values — the delay pushes the effective low-redshift slope $\gamma$ down by
about a third — and the summary prints both columns so you can see by how much.

`fit-mass` and `fit-spin` are **checks, not corrections**: they must return the
generation values to within Monte-Carlo error, and a ratio far from one there is
a bug rather than a result. The summary table prints both columns so the
distinction is visible.

## Spectral sirens

Adding `cosmology` to `blocks` switches the model from a source-frame product to
{class}`GWForge.population_fisher.model.SpectralSirenModel`, whose
density is over the **detector-frame observables**
$(m_1^{\rm det}, m_2^{\rm det}, d_L)$:

$$\ln p = \ln p_{\rm mass}(m_1^{\rm src}, m_2^{\rm src}) + \ln p(z)
  - 2\ln(1+z) - \ln\frac{\mathrm{d}d_L}{\mathrm{d}z}
  \;[+\ \ln p_{\rm spin}],$$

with $z = z(d_L; H_0, \Omega_{m,0}, w_0)$ and
$m^{\rm src} = m^{\rm det}/(1+z)$. The detector measures redshifted masses, so a
feature at a fixed *source-frame* mass appears at a detector-frame mass that
grows with $z$; combined with the measured $d_L$, that pins the
distance–redshift relation.

Working in the detector frame is the choice that has to be right. $P_{\rm det}$
is a threshold on an SNR computed from detector-frame quantities alone, so it
carries no hyper-parameter and the centred-score identity above holds exactly.

`configurations` runs a degeneracy ladder from one config file:

```ini
configurations = {'joint': {}, 'lcdm': {'w0': -1.0}, 'h0_only': {'Om0': 0.3111, 'w0': -1.0}}
```

What separates two published $\sigma(H_0)$ values is usually which rung they
quoted, not the physics, so it is worth reporting rather than choosing one
silently.

## Validation

Three checks establish that the forecast is right rather than merely
self-consistent. The first two are in `tests/test_population_fisher.py`; the
third compares against a different code entirely.

**Against the actual spread of maximum-likelihood estimates.** With no selection
applied, Term I is the plain asymptotic MLE covariance. The test draws many
independent catalogues from the Madau–Dickinson model, fits each, and checks
that the scatter of the estimates matches the forecast $\sigma$ to the
precision that many realisations can resolve.

**Against the sampler, end to end.** Independent one-year catalogues are drawn
through {class}`GWForge.population.mass.Mass` and fitted with this package's own
log density; every BGP and Default-spin parameter must come back unbiased to
within the realisation scatter, and that scatter must track the Fisher $\sigma$.
Individual parameters still fluctuate — one in eleven at $2\sigma$ is what one
in eleven at $2\sigma$ looks like.

**Against an independent implementation.** `population_fisher_seminumeric` is a
separate code for the same estimator, written before this one, with its
forecasts checked in. `validation/population_fisher_vs_seminumeric.py` runs
GWForge's models on *that code's own catalogue*, so the detected events are
identical and any difference in $\sigma$ is a difference in the code rather than
Monte-Carlo scatter. Two conventions have to be reconciled first, and the script
does both: the reference's Madau–Dickinson denominator exponent is
$\alpha + \beta$ where GWForge's is $\kappa$, so $\kappa = \alpha + \beta$ is a
change of variables needing a Jacobian rather than a rename; and the reference
tapers both component masses with one $(m_{\min}, \delta_m)$, so GWForge's
independent secondary taper has to be collapsed onto the primary's or the two
are simply different densities. The spin sector has no counterpart — the
reference defines no spin model, only a sampler — so it keeps the
finite-difference oracle above. The script writes an overlay corner per
comparison, drawn in **GWForge's** parameterisation so the parameters carry
their proper labels and `alpha`/`beta` — which mean different things in the two
codes — cannot be confused.

```bash
python validation/population_fisher_vs_seminumeric.py \
    --output-directory population_fisher_comparison
```

## Figures

The figures come from {mod}`GWForge.plotting`, an in-house corner rather than
{mod}`corner`. The reason is one requirement `corner` cannot meet: the
off-diagonal panels must be able to carry a **scatter coloured by a third
quantity** — events coloured by redshift or network SNR, the way pycbc's
posterior pages draw them. `corner`'s `plot_datapoints` colours every point the
same. The layout follows
[zeus's `cornerplot`](https://github.com/minaskar/zeus/blob/main/zeus/plotting.py):
marginals on the diagonal, contours in the lower triangle, upper triangle off.

```python
from GWForge.plotting import corner_plot, labels_for, XG_COLOUR, LVK_COLOUR

figure = corner_plot(samples, labels=labels_for(names), levels=(0.67, 0.90),
                     truths=fiducial, colour=XG_COLOUR, label="CE40+CE20+ET")
figure = corner_plot(other, labels=labels_for(names), colour=LVK_COLOUR,
                     label="H1+L1+V1", show_titles=False, fig=figure)
```

Two things differ from zeus, deliberately. Densities come from a **smoothed 2-D
histogram** rather than `seaborn.kdeplot` — fast enough to call in a loop, and
it lets the contour levels be *enclosed probabilities* computed from the sorted
histogram, so a "90% contour" contains 90% of the samples. The test suite checks
that by counting samples inside the drawn path, not by asserting a contour
exists. And the style is applied through `matplotlib.rc_context`, so drawing a
corner does not restyle everything the caller draws afterwards.

### Palettes

Two, selected with `--palette` or `palette_name=`:

| palette | kind | first two colours | notes |
| --- | --- | --- | --- |
| `tealrose` *(default)* | CARTOColors, diverging | `#009392`, `#d0587e` | the cycle walks inward from the ends, so the first two datasets get maximum separation |
| `okabe-ito` | Okabe & Ito, qualitative | `#e69f00`, `#56b4e9` | black moved last; its continuous map is a ramp built from palette members, not a published sequential scheme |

`XG_COLOUR` and `LVK_COLOUR` are the two ends of TealRose, so the discrete
comparison colours and the continuous scatter map are the same family by
construction. The hex values are stored in the module — `palettable` would be a
runtime dependency for fifteen strings.

### Colouring by a third quantity

```python
corner_plot(numpy.column_stack([m1, m2, z]),
            labels=labels_for(["mass_1_source", "mass_2_source", "redshift"]),
            colour_by=numpy.log10(catalogue.network_snr),
            colour_label=r"$\log_{10}\rho_{\rm net}$")
```

Each lower-triangle panel scatters the points on the palette's continuous map
with one shared colourbar, and the contour fill switches off by default — a
filled contour drawn over a scatter hides exactly the points the scatter exists
to show.

## Comparing networks

A saved forecast can be reloaded and overlaid without recomputing anything:

```bash
gwforge_population_fisher \
  --compare "CE40+CE20+ET=population_fisher_output/mass_snr10.npz" \
            "H1+L1+V1=population_fisher_output/mass_lvk_snr10.npz" \
  --parameters "['alpha_1', 'alpha_2', 'm_break', 'lam_0', 'beta']" \
  --output-file mass_xg_vs_lvk.pdf
```

It prints a table of $\sigma$ per parameter and network, with the ratio of the
last to the first, and draws the overlay corner. The forecasts must share their
free parameters; `mass_lvk.ini` ships alongside `mass.ini` and frees the same
ten against H1 + L1 + V1, with `delta_m` fixed in both. bilby supplies those
detectors with `minimum_frequency = 20` already set, so the SNR pass is the
same command with `--ifos H1 L1 V1 --minimum-frequency 20`. Because the
fiducials are fitted to **all injections** rather than to detections, they come
out identical for both networks and the comparison is like-for-like.

## Caveats worth quoting with any result

**Term I is the small-measurement-error limit.** It is valid when per-event
uncertainties are much smaller than the width of the population distribution.
For CE and ET that holds for most sources, but it fails where the mass function
has a sharp feature, and there Term I is *optimistic*. Terms II–V are not
implemented.

That validity is measurable rather than assumed, and
`validation/population_higher_order.py` measures it. Terms II–V are all built
from $A_k = F_k - H_k$, where $F_k = C_k^{-1}$ is the per-event measurement
Fisher and $H_k$ the Hessian of $\ln p(\theta \mid \Lambda)$ in the
*observables*, so they are small exactly when the spectral radius
$\rho(H_k C_k)$ is small. That single dimensionless number is the expansion
parameter, and it is far cheaper and far more robust than the corrections
themselves — which is also why the seminumeric reference computes
$\Gamma_{\rm II}$–$\Gamma_{\rm V}$ and then declines to use them for a
spectral-siren forecast: by nested finite differences the correction sits below
the noise floor and the raw sum returns non-positive-definite matrices.

The script needs per-event covariances, which means a `gwforge_fisher` pass over
the detected catalogue. Two things it is careful about, both easy to get wrong:
$C_k$ must be a slice of the **covariance** (marginal), not of the single-event
Fisher (conditional, smaller by orders of magnitude); and the diagnostic is the
spectral radius, not $\max_a |H_{aa}| C_{aa}$, because correlations between the
mass and distance directions make the true parameter larger than any diagonal
entry.

**`delta_m` is a standard scale, and it dominates $\sigma(H_0)$.** The taper at
$m_{\min}$ is a sharp feature at a fixed source-frame mass, so it feeds $H_0$
directly. The previous generation of this code measured $\sigma(H_0)$ swinging
by a factor of 77 across $\delta_m \in [0, 4]\,M_\odot$. Because GWForge *generates* the catalogue at a
known `delta_m`, that number is under control here rather than being set by
whatever a maximum-likelihood fit does against a hard injection cutoff — but the
sensitivity is real, and a spectral-siren forecast should scan it rather than
assume it.

**The condition number is part of the answer.** An 18-parameter forecast has
strongly correlated directions; a marginal $\sigma$ quoted without looking at
the correlation matrix can be dominated by a degeneracy rather than by the data.
`PopulationFisherResult.correlation_matrix()` is there for that.

## Output

Per SNR threshold and configuration, `gwforge_population_fisher` writes

* `{label}[_{configuration}]_snr{threshold}.npz` — the Fisher matrix, its
  covariance, the parameter names, the fiducial values, the SNR threshold and
  the event counts, which is everything `--compare` needs to reload it;
* `{label}[_{configuration}]_corner.pdf` — thresholds overlaid, restricted to
  `plot-parameters` if given (an 18-parameter corner is 324 panels and legible in
  none of them). Subsets are drawn from the **marginal** of the full covariance;
* `{label}_mass.pdf` and `{label}_redshift.pdf` — the reconstructed
  $\pi(m_1)$ and $p(z)$ with 90% credible bands, over a histogram of the
  injections. These are the quickest way to see what a table of $\sigma$s means,
  and to catch a fiducial that was fitted badly;
* `{label}_summary.txt` — the fit comparison, the catalogue cuts and the
  forecast tables, also echoed to stdout.
