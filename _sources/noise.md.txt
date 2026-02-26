# Simulating the detector noise

GWForge (at the moment) can generate a realisation of coloured Gaussian noise from a power spectral density, and it also generates strain data without any noise at all. 

For example, here is a noise configuration file for XG detectors:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']
sampling-frequency = 8192
noise = Gaussian
fft-scheme = numpy
```
Which you can run as follows:
```bash
gwforge_noise --gps-start-time 1893024018 --gps-end-time 1893187858 --config-file xg.ini --output-directory data
```
This will generate roughly a day's worth of data for an XG detector network at a sampling frequency of 8192Hz. The strain data will be a realisation of coloured Gaussian noise, and the noise power spectrum is read from here:`~/.conda/envs/gwforge-test/lib/python3.9/site-packages/GWForge/ifo/noise_curves`. The generated strain data will be stored as `*.gwf` in the output directory with names `{IFO}-{GPS-TIME}.gwf` with channel names `{IFO}:INJ`. 

Here is a short-cut in case you forget:
```bash
FrChannels <frame-file>
```
will list all the channels in the frame files
and 
```bash
FrCheck -i <frame-file>
```
Will list the gps start and end time of the frame file.

Alternatively, you can define the configuration file as follows:
```ini
[IFOS]
detectors = ['CE20', 'CE40', 'ET']
sampling-frequency = 8192
noise = Gaussian
fft-scheme = numpy
gps-start-time = 1893024018
duration = 86400
```
And remove the gps option.

```{note}
By default, `gwforge_noise` will use [`pycbc.noise.reproducable`](https://pycbc.org/pycbc/latest/html/_modules/pycbc/noise/reproduceable.html) to generate the strain data. You can use Bilby's method by defining `noise-type=bilby`.
```

```{warning}
If you are interested in generating strain data without any noise, please define `noise-type=bilby`. 
```

```{warning}
If you have MKL properly configured or you want to use CUDA for FFT you can pass `fft-scheme=mkl` or `fft-scheme=cuda`. It should work, but no promises.
```

Finally, you can define your configuration file as follows:
```
[IFOS]
detectors = ['H1', 'L1', 'V1']
psd-dict = {'H1' : 'H1.txt', 'L1': 'L1.txt', 'V1':'V1.txt'}
sampling-frequency = 8192
noise = Gaussian
fft-scheme = numpy
gps-start-time = 1893024018
duration = 86400
```
by providing a Dictionary of PSD files to use.