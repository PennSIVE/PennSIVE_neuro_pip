# Central Vein Sign (CVS) Pipeline

The CVS pipeline integrates an automated technique for the detection of the central vein sign in white matter lesions, developed by [Dr. Dworkin](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6177309/). It provides processed T1-weighted, T2-FLAIR, and T2*-EPI images, as well as subject-level CVS probabilities.

Note: An older version of this pipeline is also available as a [BIDS app](https://github.com/PennSIVE/cvs-bids-app) on [docker hub](https://hub.docker.com/r/pennsive/cvs). For instructions on using that version of the pipeline, see our [PennSIVE wiki](http://localhost:5869/pipelines.html#cvs).

## Diagram
![CVS Workflow](/pipelines/cvs/figure/cvs_pipeline.png)

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/cvs/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline contains three stages: 1) Preprocessing and CVS Probability Calculation: calculates the probability of each participant's probability of having cvs lesions, and 2) Consolidation: consolidates all participants' CVS results.

Detailed examples are provided below (all in individual mode):

### Preprocessing & CVS Calculation

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -e "*_T2star.nii.gz" -s TRUE --mode individual -c "local" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -e "*_T2star.nii.gz" -s TRUE --mode individual -c "cluster" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -e "*_T2star.nii.gz" -s TRUE --mode individual -c "singularity" --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```

-   `docker`

```bash
docker pull pennsive/neuror
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -e "*_T2star.nii.gz" -s TRUE --mode individual -c "docker" --toolpath /path/to/PennSIVE_neuro_pip 
```
**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

### Consolidation

After this stage, we also offer a second stage: consolidation, which consolidates all participants' CVS estimates into a single CSV file.

-   `local`

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project --step consolidation -c "local" --toolpath /path/to/PennSIVE_neuro_pip 
```

-   `cluster`

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/cvs/code/bash/cvs.sh -m /path/to/project --step consolidation -c "cluster" --toolpath /path/to/PennSIVE_neuro_pip 
```

## Output Data Structure
![Output](/pipelines/cvs/figure/output.png)