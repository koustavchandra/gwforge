# Injecting a signal
If you intend to inject signals into the data, follow these steps:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']
channel-dict = {'CE20':'CE20:INJ', 'CE40':'CE40:INJ', 'ET':'ET:INJ'}
sampling-frequency = 8192
minimum-frequency = 6

[Injections]
injection-file = bbh.h5
injection-type = bbh
waveform-approximant = IMRPhenomXPHM
fft-scheme = numpy
```
Similar to the [Noise](doc:noise), begin by defining the detector network. However, this time, provide the channel name of the Frame files and include the `sampling-frequency` (matching the detector data's sampling frequency) along with a `minimum-frequency`.

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
By default, GWForge uses Bilby's [`inject_signal`](https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.detector.networks.InterferometerList.html#bilby.gw.detector.networks.InterferometerList.inject_signal) to add the signal in the data coherently. Alternatively, you can just specify `injection-method = pycbc` in `[Injections]` to use the PyCBC style of injecting signal. I have incorporated this from [`pycbc.inject`](https://github.com/gwastro/pycbc/blob/master/pycbc/inject/inject.py).

Again, like before you can choose `fft-scheme` to be either `numpy`, `mkl` or `cuda`
```

Available `injection-type` (for the moment) are `bbh, bns, nsbh, imbhb, imbbh, pbh`.