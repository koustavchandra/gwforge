# Usage and Examples

1. **Population Generation:**
   Generate a gravitational wave source population by specifying parameters such as the local merger rate, distribution functions, and additional keyword arguments.

2. **Generating Detector Noise:**
   Simulate coloured Gaussian or zero noise using a provided or default power spectrum to represent the detector noise.

3. **Injecting Signals:**
   Inject gravitational wave signal(s) into the generated detector data using the previously generated population and a chosen waveform model.

4. **Adding Glitches:**
   Add incoherent sine-Gaussian glitches, drawn independently for every detector, to the generated data.

5. **Fisher Forecasts:**
   Forecast parameter-estimation uncertainties for a source or a whole population with the Fisher matrix, optionally with the second-order DALI correction for non-Gaussian posteriors.

6. **Detector Response:**
   Project signals onto the detectors with a frequency- and time-dependent antenna response, which drops the long-wavelength and static-pattern approximations that break for XG detectors.

7. **Population Forecasts:**
   Forecast how well a network measures the population itself -- the mass function, the merger-rate history, the spin distribution and, through spectral sirens, the cosmology -- with the population Fisher matrix.

8. **Workflows:**
   Chain noise generation, injections and glitches over long stretches of data as an HTCondor DAG.

I have curated some examples below for reference. Please give them a try!
```{toctree}
:caption: 'Contents:'
:maxdepth: 2

population
noise
inject
glitch
antenna
fisher
population_fisher
workflow
```

To generate the documentation, install the `docs` extra and run:
```bash
make -C docs html
```