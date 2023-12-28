# Usage and Examples

1. **Population Generation:**
   Generate a gravitational wave source population by specifying parameters such as the local merger rate, distribution functions, and additional keyword arguments.

2. **Generating Detector Noise:**
   Simulate coloured Gaussian or zero noise using a provided or default power spectrum to represent the detector noise.

3. **Injecting Signals:**
   Inject gravitational wave signal(s) into the generated detector data using the previously generated population and a chosen waveform model.

I have curated some examples below for reference. Please give them a try!
```{toctree}
:caption: 'Contents:'
:maxdepth: 2

population
noise
inject
workflow
```

To generate the documentation, just run:
```bash
sphinx-autobuild docs/source/ docs/build/html/
```