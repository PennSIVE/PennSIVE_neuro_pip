# T1/T2 Pipeline

The T1T2 pipeline generates the ratio of T1-weighted to T2-weighted signal intensity (T1/T2). The T2 sequence can be specified as a T2-weighted or FLAIR (default) image.

## Diagram
![T1/T2 Workflow](/pipelines/t1t2/figure/t1t2_pipeline.png)

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/t1t2/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline contains two stages: 1) Estimation: calculates each participants' T1/T2 ratio and 2) Consolidation: consolidates all participants' results into a single .csv file.

Detailed examples are provided below (all in individual mode):

### Estimation

(If you don't have a brain mask derived from the skull-stripping pipeline, please set `-s TRUE`.)

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w.nii.gz" -f "*FLAIR.nii.gz" --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w.nii.gz" -f "*FLAIR.nii.gz" --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w.nii.gz" -f "*FLAIR.nii.gz" --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`

```bash
docker pull pennsive/neuror
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w.nii.gz" -f "*FLAIR.nii.gz" --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip 
```

**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

### Consolidation

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project --step consolidation -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project --step consolidation -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project --step consolidation -c singularity --toolpath /path/to/PennSIVE_neuro_pip
```


-   `docker`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/t1t2/code/bash/t1t2.sh -m /path/to/project --step consolidation -c docker --toolpath /path/to/PennSIVE_neuro_pip
```

