# Installation:

GWForge is currently available only from source. Follow the steps below for the recommended installation for development or use.

## 1. Create Conda Environment:

```bash
conda create --name gwforge-venv python=3.9.18
```
This command sets up a Conda environment named gwforge-venv with Python 3.9.18 and basic packages. To activate the Conda environment, use:
```bash
conda activate gwforge-venv
```

## 2. Install GWForge
Proceed to install gwforge and its dependencies:
```bash
pip install git+https://github.com/koustavchandra/gwforge.git
```
This installs gwforge along with necessary dependencies such as [bilby](https://lscsoft.docs.ligo.org/bilby/), [gwpy](https://gwpy.github.io/docs/), [pycbc](https://pycbc.org/pycbc/latest/html/index.html) [gwpopulation](https://colmtalbot.github.io/gwpopulation/), etc. Ensure you are using the correct pip version; check by running:

```bash
which pip
```
You should see output similar to:
```bash
~/.conda/envs/gwforge-venv/bin/pip
```

## 3. Optional: Install Additional Packages
It is recommended to install lalsuite, lalsimulation, etc. This can be done using:
```bash
conda install -c conda-forge fftw lalsimulation lalsimulation-data lalsuite lalframe lalapps gitpython jupyterlab wget
```
Adjust the installation based on your specific needs.

