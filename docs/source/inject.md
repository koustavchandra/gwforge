# Injecting a signal
If you intend to inject signals into the data, follow these steps:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']

[Injections]
injection-file = bbh.h5
injection-type = bbh
waveform-approximant = IMRPhenomXO4a
waveform-minimum-frequency = 3
```
Similar to the {doc}`noise`, begin by defining the detector network. The data are the HDF5 files `gwforge_noise` writes, read from the `{IFO}:INJ` channels by default; if your files carry other channel names, add `channel-dict = {'CE20': 'CE20:STRAIN', ...}` to `[IFOS]`. The sampling frequency is taken from the data.

You must specify the path to `injection-file` and `waveform-approximant`; `injection-type` defaults to `bbh`, and the available types are `bbh, bns, bhns, imbhb, imbbh, pbh` (`nsbh` is accepted as an alias of `bhns`). The `waveform-minimum-frequency` is where every waveform starts, and the signals are generated at the reference frequency the population was drawn at, so the spins mean what `gwforge_population` meant by them.

You can execute it as follows:
```bash
gwforge_inject --config-file injections.ini --data-directory data --gps-start-time 1893024018 --gps-end-time 1893187858
```
The GPS range may also be written as `gps-start-time` and `gps-end-time` in `[Injections]`; the command line wins when both are given, which is what lets the workflow run one job per chunk from a single ini. Every signal of the population whose merger falls inside the range is injected, and the metadata of what went in — parameters and, for the bilby method, the optimal and matched-filter SNRs — are written to `injections-{type}-{start}.h5` in the data directory.

````{note}
The `waveform-approximant` must be implemented in `lalsimulation`, in the time or the frequency domain. The easiest way to check is:
```python
from pycbc.waveform import td_approximants, fd_approximants
print(td_approximants() + fd_approximants())
```
````

## Two ways of adding the signal

By default GWForge builds each signal in the frequency domain with `bilby` and adds it to the data coherently with its own detector response (see below). The signal is built over a window of consecutive 4096 s files sized from the longest signal in the population, so a signal is never longer than the window it is placed in; if the files before your start time that the window needs are missing, `gwforge_inject` stops and tells you how many it wanted. BNS and BHNS signals use the tidal source model, with the black hole of a BHNS given $\Lambda_1 = 0$.

Alternatively, specify `injection-method = pycbc` in `[Injections]` to generate the time-domain waveform and add it through LAL's `SimAddInjection`, as in [`pycbc.inject`](https://github.com/gwastro/pycbc/blob/master/pycbc/inject/inject.py). This is the realistic choice for long signals: the whole requested range is one window, every signal that overlaps it is added — including those that merge after the range ends — and whatever part of a waveform lies outside the data is simply cut, never wrapped around. The pycbc method reads `fft-scheme` (`numpy`, `mkl` or `cuda`) from `[IFOS]`, default `numpy`, and projects with LAL's `Detector.project_wave`.

## Detector response

The `bilby` injection method projects signals onto the detectors using GWForge's
frequency- and time-dependent antenna response, which drops the long-wavelength
and static-pattern approximations. Both corrections are **on by default**:

```ini
[Injections]
injection-method = bilby
earth-rotation = True
finite-size = True
```

Set either to `False` to recover bilby's historical behaviour exactly. The
`pycbc` injection method is unaffected — it uses LAL's own projection. See
{doc}`antenna` for the physics and its validation.
