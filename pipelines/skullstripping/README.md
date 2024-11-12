# Skullstripping Pipeline

The skullstripping pipeline provides several options (MASS and HD-BET) to remove skull and non-brain matter from brain images.

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/skullstripping/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline allows for two skullstripping options: mass (default) and hdbet. To use hdbet, set `-t "hdbet"`

Detailed examples are provided below (all in individual mode):

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/skullstripping/code/bash/skullstripping.sh -m /path/to/project -p sub-001 --ses ses-01 -f "*MPRAGE*.nii.gz" --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/skullstripping/code/bash/ skullstripping.sh -m /path/to/project -p sub-001 --ses ses-01 -f "*MPRAGE*.nii.gz" --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/skullstripping/code/bash/ skullstripping.sh -m /path/to/project -p sub-001 --ses ses-01 -f "*MPRAGE*.nii.gz" --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`

```bash
docker pull pennsive/neuror
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/skullstripping/code/bash/ skullstripping.sh -m /path/to/project -p sub-001 --ses ses-01 -f "*MPRAGE*.nii.gz" --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip
```

**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

## Output Data Structure
![Output](/pipelines/skullstripping/figure/output.png)

