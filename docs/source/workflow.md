# A workflow for gwforge
~~It can be useful~~It is definitely useful to use `gwforge_workflow` to generate long data periods with/without signals rather than manually generating them individually.

Here is the page devoted to step-wise:
1. Generate the source population by following the instructions in {doc}`population`
2. Define your noise-configuration file by following the instructions in {doc}`noise`
3. Define the injection-configuration file(s) by following the instructions in {doc}`inject`.

Once you have defined them, you simply run:
```bash
# Activate the environment gwforge is installed in
conda activate gwforge-venv || { echo "Failed to activate Conda environment." >&2; exit 1; }

# Set output directory
output_directory=${HOME}/projects/XG/test-workflow

# Workflow submission script
gwforge_workflow \
    --gps-start-time 1893024018 \
    --gps-end-time 1893187858 \
    --output-directory "${output_directory}/output" \
    --noise-configuration-file "${output_directory}/xg.ini" \
    --bbh-configuration-file "${output_directory}/injections.ini" \
    --accounting-group ligo.dev.o4.cbc.pe.lalinference \
    --workflow-name trial \
    --submit-now
```
and that's it! It will create a directory called `output` (or `output_1`, `output_2`, ... if it already exists) with the following structure.
```
output/
├── trial.condor
├── logs
│   ├── noise-0.log  noise-0.err  noise-0.out
│   └── bbh-0.log    bbh-0.err    bbh-0.out
└── submit
    ├── noise-0.sub
    └── bbh-0.sub
```
and submit all the jobs in the submit directory to HTCondor. The data directory `output/data/<IFO>/` is created by the noise jobs, which write the HDF5 files `<IFO>-<start>-<duration>.h5`; the injection jobs then modify them in place. HTCondor stores the log files in `logs` (or wherever `--log-directory` points). [I love to call this a donkey-sus.]

The range is split into chunks of 24576 s, one noise job and one injection job per source type per chunk, numbered `-0`, `-1`, ... Neighbouring chunks share files, so within each source type the odd-numbered jobs run first and the even-numbered ones depend on them; across source types the chain is noise, then BBH, then BHNS, then BNS, skipping whatever you did not ask for.

In case you want to add BNS and/or BHNS, pass the options:
```bash
--bns-configuration-file <bns-file.ini> --bhns-configuration-file <bhns-file.ini>
```
(`--nsbh-configuration-file` is still accepted for the latter.)
