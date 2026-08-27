# Usage and Examples

1. **Population Generation:**
   Generate a gravitational wave source population by specifying parameters such as the local merger rate, distribution functions, and additional keyword arguments.

2. **Generating Detector Noise:**
   Simulate coloured Gaussian or zero noise using a provided or default power spectrum to represent the detector noise.

3. **Injecting Signals:**
   Inject gravitational wave signal(s) into the generated detector data using the previously generated population and a chosen waveform model.

4. **Fisher Forecasts:**
   Forecast parameter-estimation uncertainties for a source or a whole population with the Fisher matrix, optionally with the second-order DALI correction for non-Gaussian posteriors.

5. **Detector Response:**
   Project signals onto the detectors with a frequency- and time-dependent antenna response, which drops the long-wavelength and static-pattern approximations that break for XG detectors.

I have curated some examples below for reference. Please give them a try!
```{toctree}
:caption: 'Contents:'
:maxdepth: 2

population
noise
inject
antenna
fisher
population_fisher
workflow
```

To generate the documentation, just run:
```bash
sphinx-autobuild docs/source/ docs/build/html/
```