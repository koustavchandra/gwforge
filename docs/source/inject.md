# Injecting a signal
If you intend to inject signals into the data, follow these steps:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']
channel-dict = {'CE20':'CE20:INJ', 'CE40':'CE40:INJ', 'ET':'ET:INJ'}
sampling-frequency = 8192

[Injections]
injection-file = bbh.h5
injection-type = bbh
waveform-approximant = IMRPhenomXO4a
waveform-minimum-frequency = 3
```
Similar to the [Noise](doc:noise), begin by defining the detector network. However, this time, provide the channel name of the Frame files and include the `sampling-frequency` (matching the detector data's sampling frequency) along with a `waveform-minimum-frequency` for the signals.

You must specify the path to `injection-file`, `injection-type` and `waveform-approximant` as extras.

You can execute it as follows:
```bash
gwforge_inject --config-file injections.ini --data-directory data --gps-start-time 1893024018 --gps-end-time 1893187858
```

```{note}
The `waveform-approximant` must be implemented in `lalsimulation`. The easiest way to check the waveform availability is to execute:
```bash
from pycbc.waveform import td_approximants
print(td_approximants())
```

```{warning}
By default, GWForge adds the signal coherently in the frequency domain with its own detector response (see below). Alternatively, specify `injection-method = pycbc` in `[Injections]` to add a time-domain waveform through LAL's `SimAddInjection`, as in [`pycbc.inject`](https://github.com/gwastro/pycbc/blob/master/pycbc/inject/inject.py): a signal longer than the data is then cut at the segment edges rather than wrapped, and all signals that overlap the segment are added. The pycbc method reads `fft-scheme` (`numpy`, `mkl` or `cuda`) from the `[IFOS]` section.
```

Available `injection-type` (for the moment) are `bbh, bns, bhns, imbhb, imbbh, pbh`; `nsbh` is accepted as an alias of `bhns`.
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
[Detector response](antenna.md) for the physics and its validation.
