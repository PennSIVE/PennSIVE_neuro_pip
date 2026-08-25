# Automated Detection of Paramagnetic Rim Lesion (APRL) Pipeline

The APRL pipeline integrates an automated technique for paramanetic rim lesion (PRL) detection, developed by [Dr. Lou](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8503902/). It provides processed T1-weighted, T2-FLAIR, and T2star-Phase images, as well as white matter lesion masks and the probability of each lesion being a PRL. (The current pipeline uses a pre-trained MIMoSA model and APRL model)

## Diagram
![APRL Workflow](/pipelines/PRL/figure/aprl_pipeline.png)

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/PRL/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline contains three stages: 1) Preprocessing: processes MRI images to prepare for PRL probability calculation, 2) PRL Probability Calculation: calculates the probability of each lesion being a PRL, and 3) Consolidation: consolidates all participants' PRL results.

Detailed examples are provided below (all in individual mode):

### Preprocessing

(If you don't have a brain mask derived from the skull-stripping pipeline, please set `-s TRUE`.)

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --phase "*phase_T2star_UNWRAPPED.nii.gz" --mode individual -c "local" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --phase "*phase_T2star_UNWRAPPED.nii.gz" --mode individual -c "cluster" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://russellshinohara/pennsive_amd64_cputorch
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --phase "*phase_T2star_UNWRAPPED.nii.gz" --mode individual -c "singularity" --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`

```bash
docker pull russellshinohara/pennsive_amd64_cputorch
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --phase "*phase_T2star_UNWRAPPED.nii.gz" --mode individual -c "docker" --toolpath /path/to/PennSIVE_neuro_pip 
```

**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

### PRL Probability Calculation

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual -c "local" --toolpath /path/to/PennSIVE_neuro_pip --step PRL_run
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual -c "cluster" --toolpath /path/to/PennSIVE_neuro_pip --step PRL_run
```

-   `singularity` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual -c "singularity" --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path --step PRL_run
```


-   `docker`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual -c "docker" --toolpath /path/to/PennSIVE_neuro_pip --step PRL_run
```

### Consolidation

-   `local`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project --toolpath /path/to/PennSIVE_neuro_pip --step consolidation -c "local"
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/PRL/code/bash/PRL.sh -m /path/to/project --toolpath /path/to/PennSIVE_neuro_pip --step consolidation -c "cluster"
```

## Output Data Structure
![Output](/pipelines/PRL/figure/output.png)

