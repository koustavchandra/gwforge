# Detector response

Bilby's `Interferometer.get_detector_response` — which GWForge used everywhere until now — makes two approximations that are harmless for LIGO/Virgo and wrong for XG detectors:

1. **The long-wavelength approximation.** The detector tensor is treated as frequency-independent, i.e. the wave is assumed not to change appreciably while light traverses an arm. The correction scales with $fL/c$ and becomes important near the *free spectral range* $c/2L$ — **3.75 kHz for CE's 40 km arms**, which is inside the analysis band.

2. **A static antenna pattern.** $F_{+,\times}$ and the geocentre-to-vertex delay are frozen at `geocent_time`. ET and CE see BNS inspirals lasting *hours*, over which the Earth rotates through a large angle, so both are really functions of frequency through the time-frequency relation $t(f)$.

{mod}`GWForge.ifo.antenna` drops both. It is **on by default**.

## Using it

Both corrections are controlled by two flags, available wherever GWForge projects a signal onto a detector.

For injections, in the `[Injections]` section of the config file:
```ini
[Injections]
injection-method = bilby
earth-rotation = True
finite-size = True
```

For `gwforge_optimal_snr`, on the command line:
```bash
gwforge_optimal_snr --injection-file injections.h5 --output-file snrs.h5 \
                    --no-earth-rotation --no-finite-size   # restores the old behaviour
```

From Python:
```python
from GWForge.ifo import antenna

signal = antenna.detector_response(
    interferometer, polarizations, parameters,
    earth_rotation=True, finite_size=True,
)
```

```{warning}
Because these default to `True`, injections and optimal SNRs computed with this version differ from those produced by earlier versions. Set both to `False` to reproduce the old numbers exactly — that combination is asserted to reproduce bilby to a relative $2.5\times10^{-16}$.
```

```{note}
The `pycbc` injection method is unaffected: it uses LAL's own projection via `Detector.project_wave`, which is a separate code path.
```

## The physics

The implementation follows [Baral et al. (arXiv:2304.09889)](https://arxiv.org/abs/2304.09889), with polarisation tensors as defined in [Nishizawa et al. (arXiv:0903.0528)](https://arxiv.org/abs/0903.0528) — the same convention bilby uses.

### Finite arm length

Each arm gets its own transfer function, from [Rakhmanov, Romano & Whelan (arXiv:0808.3805)](https://arxiv.org/abs/0808.3805):

$$
D(x, y) = \tfrac{1}{2}\left[
    e^{-i\pi x(1+y)}\,\mathrm{sinc}\big(x(1-y)\big)
  + e^{ i\pi x(1-y)}\,\mathrm{sinc}\big(x(1+y)\big)\right]
$$

with $x = fL/c$ and $y = -\hat\Omega\cdot\hat{a}$. The two half-outer-products $\tfrac{1}{2}\hat{x}\otimes\hat{x}$ and $\tfrac{1}{2}\hat{y}\otimes\hat{y}$ — whose difference is exactly bilby's `detector_tensor` — are contracted separately so each carries its own $D$.

Two things worth knowing about the size of this effect:

* $D(0, y) = 1$ exactly, so switching `finite_size` off is a genuine limit and not an approximation.
* To first order $D \approx 1 - i\pi x y$. The leading correction is a **phase**, not an amplitude change — it is the extra light travel to the arm midpoint. That is why it grows linearly in $f$ rather than quadratically, and why much of it is degenerate with the coalescence time. The genuinely new content appears at $\mathcal{O}(x^2)$.
* At $y = 0$ the response is $\mathrm{sinc}(x)\cos(\pi x)$, which vanishes at $f = c/2L$.

### Earth rotation

The frequency bin at $f$ was emitted $\tau(f)$ seconds before merger, so the Earth's orientation there is its orientation at `geocent_time - tau(f)`. GWForge uses the 3.5PN time-to-coalescence of [arXiv:0907.0700](https://arxiv.org/abs/0907.0700) eq. (3.8b). Mode $m$ reaches frequency $f$ when the quadrupole is at $2f/m$, i.e. $\tau_m(f) = \tau_2(2f/m)$.

The sidereal angle is advanced at the **sidereal** rate, $2\pi/86164.09\,\mathrm{s}$.

## Validation

`tests/test_antenna_response.py` layers its checks outward from claims that cannot be wrong.

| Check | Result |
|---|---|
| Both corrections off $\equiv$ bilby's `get_detector_response`, CE40 and all three ET arms | $2.5\times10^{-16}$ relative |
| $D(x,y)$ vs brute-force integration of the round-trip light path (shares no code) | $<10^{-9}$ |
| $D(0,y) = 1$; null at $f = c/2L$ | exact |
| First-order term is $-i\pi xy$, residual scales as $x^2$ | confirmed |

The first row is the anchor: it pins every sign, the handedness of $\psi$, and the delay convention against code GWForge already relies on, so a convention error surfaces there rather than being absorbed into the corrections layered on top.
