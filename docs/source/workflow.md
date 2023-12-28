# A workflow for gwforge
~~It can be useful~~It is definitely useful to use `gwforge_workflow` to generate long period of data with/without signals than manually generating them one by one.

Here is the page devoted to step-wise:
1. Generate the source population by the following the instructions in [population](doc:population)
2. Define your noise-configuration file by following the instructions in [noise](doc:noise)
3. Define the injection-configuration file(s) by following the instructions in [inject](doc:inject).

Once you have defined them, you simply run:
```bash
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
# Activate Conda environment
conda activate gwforge-test || { echo "Failed to activate Conda environment." >&2; exit 1; }

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
and that's it!. It will create a directory called `output` with the following structure
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
and submit all the jobs in submit directory to HTCondor. The resulting log files will be stored in logs, frame files in the respective ifo directory. [I love to call this a donkey-sus.]


In case you want to add bns and/or nsbh, pass the options:
```
--bns-configuration-file <bns-file.ini> --nsbh-configuration-file <nsbh-file.ini>
```
