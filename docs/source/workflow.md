# A workflow for gwforge
~~It can be useful~~It is definitely useful to use `gwforge_workflow` to generate long data periods with/without signals rather than manually generating them individually.

Here is the page devoted to step-wise:
1. Generate the source population by following the instructions in [population](doc:population)
2. Define your noise-configuration file by following the instructions in [noise](doc:noise)
3. Define the injection-configuration file(s) by following the instructions in [inject](doc:inject).

Once you have defined them, you simply run:
```bash
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
# Activate Conda environment
conda activate gwforge-venv || { echo "Failed to activate Conda environment." >&2; exit 1; }

# Set environment variables
export LAL_DATA_PATH=/ligo/home/ligo.org/koustav.chandra/soft/hdf5_data/

# Set output directory
output_directory=/ligo/home/ligo.org/koustav.chandra/projects/XG/test-workflow

# Workflow submission script
gwforge_workflow \
    --gps-start-time 1893024018 \
    --gps-end-time 1893187858 \
    --output-directory "${output_directory}/output" \
    --noise-configuration-file "${output_directory}/xg.ini" \
    --bbh-configuration-file "${output_directory}/injections.ini" \
    --workflow-name trial \
    --submit-now
```
and that's it! It will create a directory called `output` with the following structure.
```
output/
├── [ 4.0K]  data
│   ├── [ 4.0K]  CE20
│   ├── [ 4.0K]  CE40
│   ├── [ 4.0K]  ET1
│   ├── [ 4.0K]  ET2
│   └── [ 4.0K]  ET3
├── [ 4.0K]  logs
└── [ 4.0K]  submit
```
and submit all the jobs in the submit directory to HTCondor. HTCondor will store the resulting log files in logs and frame files in the respective IFO directory. [I love to call this a donkey-sus.]


In case you want to add BNS and/or NSBH, pass the options:
```
--bns-configuration-file <bns-file.ini> --nsbh-configuration-file <nsbh-file.ini>
```
