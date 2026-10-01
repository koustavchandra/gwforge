# Installation:

GWForge is currently available only from source. Follow the steps below for the recommended installation for development or use.

## 1. Create Conda Environment:

```bash
conda create --name gwforge-venv python=3.11
```
This command sets up a Conda environment named gwforge-venv with Python 3.11 and basic packages. To activate the Conda environment, use:
```bash
conda activate gwforge-venv
```

## 2. Install Optional dependencies

It is recommended to install lalsuite, lalsimulation, etc. This can be done using:
```bash
conda install -c conda-forge fftw lalsimulation lalsimulation-data lalsuite
```
Adjust the installation based on your specific needs.

## 3. Finally install the package itself

Proceed to install gwforge and its dependencies:
```bash
pip install git+https://github.com/koustavchandra/gwforge.git
```
This installs gwforge along with its dependencies: numpy, scipy, h5py, matplotlib, pandas, astropy, [bilby](https://lscsoft.docs.ligo.org/bilby/), lalsuite, emcee, corner, [pycbc](https://pycbc.org/pycbc/latest/html/index.html), [gwpy](https://gwpy.github.io/docs/) and rich. Ensure you are using the correct pip version; check by running:

```bash
which pip
```
You should see output similar to:
```bash
~/.conda/envs/gwforge-venv/bin/pip
```

### Development install

```bash
git clone https://github.com/koustavchandra/gwforge.git
cd gwforge
pip install -e .
```
Three extras are available: `pip install -e ".[docs]"` for building these pages, `".[test]"` for pytest and `".[waveforms]"` for `pyseobnr`. The HTCondor workflow generator needs `condor_submit_dag` on the cluster only; nothing else is platform specific.
