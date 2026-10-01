# Fisher and DALI Forecasts

`gwforge_fisher` forecasts how well a network will measure the parameters of a
source. It works with **any** waveform `bilby.gw.WaveformGenerator` can produce.

At `order = 1` this is the Fisher matrix of
[Vallisneri (2008)](https://arxiv.org/abs/gr-qc/0703086),

$$\Gamma_{ab} = \langle \partial_a h \mid \partial_b h\rangle, \qquad
\Sigma = \Gamma^{-1},$$

with the noise-weighted inner product
$\langle a \mid b\rangle = 4\,\mathrm{Re}\int \tilde a^* \tilde b / S_n(f)\,\mathrm{d}f$
summed over the detectors of the network.

At `order = 2` it adds the doublet-DALI correction of
[Sellentin, Quartin & Amendola (2014)](https://arxiv.org/abs/1401.6892), as
applied to gravitational waves by
[GWDALI](https://arxiv.org/abs/2307.10154):

$$\log L = -\tfrac{1}{2} F_{ab}\Delta^a\Delta^b
           - \tfrac{1}{2} G_{abc}\Delta^a\Delta^b\Delta^c
           - \tfrac{1}{8} H_{abcd}\Delta^a\Delta^b\Delta^c\Delta^d,$$

with $F_{ab} = \langle h_{,a}|h_{,b}\rangle$,
$G_{abc} = \langle h_{,ab}|h_{,c}\rangle$ and
$H_{abcd} = \langle h_{,ab}|h_{,cd}\rangle$.

```{note}
`order = 3` (triplet-DALI) is **not** implemented and raises an error. Grouping
the expansion by highest derivative, the terms above already exhaust everything
expressible with first and second derivatives; the next group consists entirely
of *third* derivatives of the waveform.
```

## Quick start

```bash
gwforge_fisher --config-file fisher.ini
```

```ini
[Injection_Parameters]
; either write the source out here ...
chirp_mass = 28.0
symmetric_mass_ratio = 0.24
chi_1 = 0.3
chi_2 = -0.2
luminosity_distance = 1500.0
theta_jn = 0.7
psi = 0.9
phase = 1.3
geocent_time = 1893024018.0
ra = 1.375
dec = -1.2108

; ... or point at a gwforge_population file instead
; injection-file = bbh.h5
; start-index = 0
; end-index = 100

[Waveform_Generator]
waveform-approximant = SEOBNRv5HM_ROM
waveform-minimum-frequency = 20
frequency-domain-source-model = lal_binary_black_hole

[IFOS]
detectors = ['CE20', 'CE40', 'ET']
sampling-frequency = 4096

[Fisher]
order = 1

[Output]
output-file = fisher.pkl
plot-corner = True
```

```{important}
**Point `injection-file` at a catalogue and the redundant parametrisations are
stripped for you.** A `gwforge_population` file carries `mass_1`, `mass_2`,
`chirp_mass`, `symmetric_mass_ratio`, `mass_ratio`, `total_mass` *and* the
`_source` variants, all at once. bilby's
`convert_to_lal_binary_black_hole_parameters` prefers the component masses and
never looks at `chirp_mass`, so a Fisher over `(chirp_mass,
symmetric_mass_ratio)` would displace keys the waveform ignores: the mass
derivatives came back **exactly zero**, the Fisher was singular in the mass
block, and nothing raised — $\sigma(\mathcal{M})$ was wrong by six orders of
magnitude. {func}`GWForge.fisher.parameters.strip_shadowed_parameters` now
removes whatever duplicates a varied parameter, and logs what it dropped. The
same shadowing applies to spins, where `a_1` wins over `chi_1`.

If you build a {class}`~GWForge.fisher.matrix.FisherMatrix` yourself, this is
handled in its constructor; if you assemble derivatives directly, call the
helper first.
```

`frequency-domain-source-model` accepts `lal_binary_black_hole`,
`lal_binary_neutron_star`, `gwsignal_binary_black_hole`,
`lal_eccentric_binary_black_hole_no_spins` and `EFPE_binary_black_hole` now.
You can plugin any waveform generator pretty easily.

Three named parameter sets are available, and picking the right one matters:
`DEFAULT_PARAMETERS` (aligned spin), `PRECESSING_PARAMETERS` and
`ECCENTRIC_PRECESSING_PARAMETERS`. Using the aligned-spin set with a precessing
waveform is a *silent* mistake — bilby's conversion turns `chi_1`/`chi_2` into
`a_1`/`a_2` with both tilts zero, so the waveform generates happily and you
forecast an aligned-spin binary.

`SEOBNRv5HM`/`SEOBNRv5PHM` need `pyseobnr` (`pip install -e .[waveforms]`), and
`pyEFPEHM` is installed from
[git](https://github.com/gw-models/pyEFPEHM). `pyEFPEHM` ignores
`minimum_frequency`: set `f22_start` in the waveform arguments instead, and set
it *below* the analysis band, because its eccentric harmonics radiate below the
quadrupole start frequency. Similarly, for `TEOBResumS` install `teobresums`.

The segment length is chosen per source from
{func}`GWForge.conversion.get_safe_signal_durations`, the same estimate
`gwforge_optimal_snr` uses, so a forecast and an SNR see the same signal.

## Step sizes (Claudio speaking)

Finite-difference steps are calibrated per source so that one step changes the
waveform by a fixed fraction, $\lVert\Delta h\rVert/\lVert h\rVert \approx 10^{-2}$.

This matters more than it sounds. LAL places each higher mode's onset at a
mass-dependent frequency, so differencing across that edge contributes a
spurious term growing like $1/\text{step}$. With a step an order of magnitude
too small, $\sigma(\eta)$ for `IMRPhenomXHM` is wrong by **150%**; with the
calibration on, rescaling the starting guess by a factor of thirty moves every
$\sigma$ by less than half a percent.

### Measuring the waveform's own numerical noise

A finite difference cannot see beneath the noise of the function it differences.
GWForge measures that noise directly, before computing any derivative, using the
ECNoise procedure of [Moré & Wild
(2011)](https://epubs.siam.org/doi/10.1137/100786125): sample a smooth scalar
functional of the waveform at closely spaced points and build the finite
difference table. A smooth function's $k$-th differences shrink like $h^k$;
noise does not shrink at all.

The functional is the noise-weighted overlap $\langle h(x)\mid h_0\rangle$ —
smooth in the parameter (unlike $\lVert\Delta h\rVert$, which has a kink at
zero) and PSD-weighted exactly as the Fisher matrix is. It is kept complex,
because the real part is stationary at the fiducial point for a phase-like
parameter.

`WaveformDerivatives.noise_ratios` reports the noise divided by the change one
step produces — the approximate relative error of that derivative — and it is
stored in the results pickle so the caveat travels with the forecast. Measured
worst case per model:

| model | worst relative error |
|---|---|
| `IMRPhenomXHM` | 4e-6 — machine precision |
| `SEOBNRv5HM` | 3e-3 on masses and spins; **5e-14** on `theta_jn` and `phase` |
| `SEOBNRv5PHM` | 2e-4 on the precession angles |
| `pyEFPEHM` | 2e-4 on `chirp_mass` |

The SEOBNRv5HM row is the one that shows the estimator is measuring something
real: masses and spins pass through the adaptive EOB ODE solve and pick up
scatter, while `theta_jn` and `phase` are applied analytically outside it and
come out ten orders of magnitude cleaner.

The estimate is used only to *raise* a step that is too small to clear the
noise. For a closed-form waveform it never binds.

```{note}
Noise is not the dominant limit in practice. Scanning IMRPhenomXHM's `chi_1`
step from 3e-4 to 1e-2 moves `sigma(chi_1)` by 0.03%; the same scan for
SEOBNRv5HM moves it by 70%, with only a narrow flat region around 3e-4 to 1e-3.
That is ten times larger than the measured noise, so it is truncation — the
derivative genuinely varying with the step — not scatter. A Fisher forecast with
an EOB waveform is correspondingly less well determined than one with a
closed-form model, and no step-selection strategy changes that.
```

Override the calibration if you need to:

```ini
[Fisher]
calibrate-steps = False
step-sizes = {'chirp_mass': 1e-4, 'chi_1': 1e-3}
target-change = 1e-2
```

## Response corrections

```ini
[Fisher]
earth-rotation = True
finite-size = True
```

Both default to `True`, matching how `gwforge_inject` projects the same signal.
Turning both off recovers bilby's static long-wavelength response.

## Priors

```ini
[Priors]
prior-file = my.prior          ; optional at order 1, REQUIRED at order 2
enforce-physicality = False
```

* **Order 1** — optional. The reported covariance is always the pure-likelihood
  $\Gamma^{-1}$. A *Gaussian* prior additionally contributes $1/\sigma^2$ to the
  diagonal; a uniform prior carries no Fisher information in its interior.
* **Order 2** — required. The DALI posterior is a truncated Taylor expansion,
  which is only meaningful inside bounds; without them the cubic and quartic
  terms produce spurious modes that are artefacts of the truncation.

`enforce-physicality` additionally clips priors to physical ranges
($q \le 1$, $|\chi| < 1$, angles in range).

```{warning}
A prior narrower than the Fisher $\sigma$ silently *sets* the posterior width
rather than bounding it. GWForge warns when a prior reaches fewer than three
sigma from the fiducial value.
```

## Output

The results pickle is a list of one record per source, each holding the
fiducial parameters, the network and per-detector Fisher matrices, the
covariance, the sigmas, the optimal SNRs, the condition number and inversion
error, the calibrated steps — plus `doublet`, `quartic` and posterior `samples`
at order 2.

Order 2 posteriors are non-Gaussian, so they have no covariance matrix and no
closed-form marginals. A corner plot *is* a set of marginals, so GWForge samples
`log L_DALI` with emcee through `bilby.run_sampler`. All the waveform
evaluations happen while building the tensors; evaluating the DALI likelihood
afterwards is pure tensor contraction.

## Two waveforms, one binary

`validation/compare_waveforms.py` forecasts the same aligned-spin binary with
`IMRPhenomXHM` and `SEOBNRv5HM` and overlays the two error ellipses. Both models
describe the source equally well, so the difference between them is a lower
bound on the waveform systematic in any Fisher forecast.

```bash
python validation/compare_waveforms.py --outdir waveform_comparison
```

For a 28 M☉ chirp mass source at network SNR 329 on CE40+CE20+ET:

| quantity | agreement |
|---|---|
| network SNR | **0.07%** |
| `luminosity_distance`, `ra`, `dec`, `theta_jn`, `psi` | 0.6 – 1.6% |
| `chirp_mass` | 35% |
| `geocent_time`, `phase` | factor 1.8 – 2.4 |
| `chi_1`, `chi_2` | factor 4 – 5 |

The SNRs agree almost exactly and the distance/orientation block overlays
closely; the mass–spin–time–phase block, which is where the two models actually
differ, does not. This is a genuine waveform systematic, not a numerical
artefact — the measured noise in both models is three orders of magnitude too
small to explain it.