# Automatic WM Lesion Segmentation (MIMoSA) Pipeline

The MIMoSA pipeline integrates an automated technique for white matter lesion segmentation, developed by [Dr. Valcarcel](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030441/). It provides processed T1-weighted and T2-FLAIR images, as well as a white matter lesion mask. (The current pipeline uses a pre-trained MIMoSA model, which was trained using 3T T1-weighted and T2-FLAIR images as the input)

Note: An older version of this pipeline is also available as a [BIDS app](https://github.com/PennSIVE/mimosa) on [docker hub](https://hub.docker.com/r/pennsive/mimosa). For instructions on using that version of the pipeline, see our [PennSIVE wiki](http://localhost:5869/pipelines.html#mimosa).

## Diagram
![MIMoSA Workflow](/pipelines/mimosa/figure/mimosa_pipeline.png)

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/mimosa/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). Detailed examples are provided below (all in individual mode):

(If you don't have a brain mask derived from the skull-stripping pipeline, please set `-s TRUE`.)

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/mimosa/code/bash/mimosa.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --mode individual -c "local" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/mimosa/code/bash/mimosa.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --mode individual -c "cluster" --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/mimosa/code/bash/mimosa.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --mode individual -c "singularity" --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`

```bash
docker pull pennsive/neuror
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/mimosa/code/bash/mimosa.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" --mode individual -c "docker" --toolpath /path/to/PennSIVE_neuro_pip 
```

**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.


## Output Data Structure
![Output](/pipelines/mimosa/figure/output.png)

