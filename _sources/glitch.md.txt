# Adding glitches
Real detector data carry non-Gaussian transients aka glitches that are mostly local to one detector. `gwforge_glitch` simulates such glitches as sine-Gaussian farts

$$h(t) = A \, e^{-(t - t_0)^2 / \tau^2} \cos\left(2\pi f_0 (t - t_0) + \phi\right), \qquad \tau = \frac{Q}{\sqrt{2}\,\pi f_0},$$

 following Kat Chatziioannou et al. 2021, [arXiv:2101.01200](https://arxiv.org/abs/2101.01200), and can add them incoherently in every detector. Currently, Claudio and I have shipped only two varieties --- blips and tomtes:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']

[Blip]
rate = 1
frequency-range = [50, 500]
q-range = [2, 6]
snr-range = [5, 200]
snr-index = 2.3
phases = [0, 3.141592653589793]
seed = 0

[Tomte]
rate = 3
frequency-range = [30, 65]
q-range = [4, 12]
snr-range = [5, 200]
snr-index = 2.2
phases = [3.141592653589793]
seed = 0
```

These glitches are Poisson distributed with mean `rate` per hour times the file duration (`rate` may be a per-detector dict, `{'CE40': 1.3, 'ET1': 0.65, ...}`); $f_0$ is log-uniform in `frequency-range`; the optimal SNR follows the power law $p(\rho) \propto \rho^{-\alpha}$ with $\alpha$ = `snr-index` on `snr-range` (`snr-index = 1` is log-uniform); $Q$ is uniform in `q-range`; $\phi$ is drawn from `phases` when given, else uniform; and $t_0$ is uniform but kept six envelope widths inside the file so no burst is cut by a file edge. The amplitude $A$ is then set so that the optimal SNR against the detector's PSD equals the drawn value.

## Where the numbers come from
They are ~~stolen~~ motivated by public Gravity Spy O3a Omicron tables ([Zenodo 5649212](https://zenodo.org/records/5649212), classifications with confidence above 0.9) and the glitch literature:

| class | rate [per hour] | peak frequency, 5–95 % [Hz] | Omicron SNR power-law index | duration |
|---|---|---|---|---|
| Blip | 1.3 at H1, 0.65 at L1 (about 2 in O1/O2) | 93–412 at H1, 111–630 at L1 | 2.2–2.4 | $\mathcal{O}(10)$ ms |
| Tomte | 0.22 at H1, 6.2 at L1 | 32–61 | 2.0–2.4 | about 0.25 s |


You can execute it as follows:
```bash
gwforge_glitch --config-file blip_tomte.ini --data-directory data --gps-start-time 1893024018 --gps-end-time 1893187858
```
The GPS range may also be written as `gps-start-time` and `gps-end-time` in `[IFOS]`; the command line wins when both are given. The draws are keyed on `(seed, class, detector, file start)`, so the same seed and range reproduce the same glitches, and the per-chunk jobs of the {doc}`workflow` (pass `--glitch-configuration-file`) produce the same glitches as one long job. What went in is written per detector to `glitches-{start}.h5` in the data directory, with the datasets `class, time, frequency, q, snr, phase, amplitude`.

```{note}
Glitches are incoherent: there is no sky position and no projection, and the detectors of a network, including the three ET interferometers, each draw their own. Run `gwforge_glitch` after `gwforge_inject`, as the workflow does; the order does not change the data, but it keeps the SNRs `gwforge_inject` records free of glitch power.
```
