# BIDS Pipeline

The BIDS pipeline converts DICOM images to NIfTI and organizes the files in BIDS format. It uses the [heudiconv](https://heudiconv.readthedocs.io/en/latest/) DICOM converter and an RShiny app for heuristic customization.

## Data Structure
This pipeline requires all DICOMs to be in a folder called 'original_data' within the main data directory. An example is provided below:

![Data Structure](/pipelines/bids/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide two types of scenarios for running the pipeline: `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline contains three stages: 1) Heuristic: prepares heuristic template, 2) Customization: launches RShiny app for heuristic customization, and 3) BIDS: runs DICOM to NIfTI conversion and format into BIDS structure.

Detailed examples are provided below (all in individual mode):

### Heuristic

-   `singularity` 

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /path/to/project -p sub-001 -s ses-01 --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip

bash /home/ehorwath/projects/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /home/ehorwath/projects/quy_data -p sub-001 -s ses-01 --mode individual -c singularity --toolpath /home/ehorwath/projects/PennSIVE_neuro_pip
```
**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path. Otherwise, you will need to pull the below Singularity image before running the above command.

```bash
singularity pull -F $sin_path docker://nipy/heudiconv
```

-   `docker`

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip
```
**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `dockerpath`, which has been given a default path. Otherwise, you will need to pull the below Docker image before running the above command.

```bash
docker pull nipy/heudiconv
```


### Customization

In this step, an RShiny app will launch to customize the heuristic template created in the last step. **This step only runs in batch mode and does not need a container specification.** If you are unable to connect to the app from your terminal, try running this step in VSCode. 

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /path/to/project --step customization --toolpath /path/to/PennSIVE_neuro_pip
```

The Shiny app allows you to edit the heuristic file for all subjects in the original_data folder or for each subject individually. 

To begin, under Choose Python Script, load in the heuristic.py file in the template folder. 

To review each subjects' DICOM info and edit the heuristic on a **subject-level basis**, the DICOM Info Review will load each subject's info by clicking Next and Previous in the DICOM Selection. Edits can be made in the Update Heuristic Script section and finalized by clicking **Update Script**. 

**Group-level changes** to the heuristic can be made by edits to the Update Heuristic Script section, and when finished, clicking **Update All Scripts**. This will apply those changes to all subjects in the folder.


### BIDS

-   `singularity` 

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual --step bids -c singularity --toolpath /path/to/PennSIVE_neuro_pip
```
**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path. Otherwise, you will need to pull the below Singularity image before running the above command.

```bash
singularity pull -F $sin_path docker://nipy/heudiconv
```

-   `docker`

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/bids/code/bash/bids_curation.sh -m /path/to/project -p sub-001 --ses ses-01 --mode individual --step bids -c docker --toolpath /path/to/PennSIVE_neuro_pip
```
**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `dockerpath`, which has been given a default path. Otherwise, you will need to pull the below Docker image before running the above command.

```bash
docker pull nipy/heudiconv
```


## Output Data Structure
![Output](/pipelines/bids/figure/output.png)