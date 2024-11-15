# Lesion Count Pipeline

The lesion_count pipeline provides two options to count the number of lesions present in MRI images using the "DworCount" method, developed by [Dr. Jordan Dworkin](https://www.ajnr.org/content/early/2018/02/22/ajnr.A5556) and the connected components method.

## Diagram
![Lesion Count Workflow](/pipelines/lesion_count/figure/lesion_count_pipeline.png)

## Data Structure
This pipeline requires all neuroimages to be organized in the BIDS format. An example is provided below:

![Data Structure](/pipelines/lesion_count/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline allows for three count options: DworCount (set --method dworcount), connected components (set --method cc), or both (default; set --method both). 

The pipeline contains three stages: 1) Preparation: prepares data for count by running N4 bias correction, skullstripping, registration, WhiteStripe normalization, and mimosa segmentation, 2) Count: counts number of lesions using specified method, 3) Consolidation: consolidates all participants’ results into a single .csv file.

Detailed examples are provided below (all in individual mode):

### Preparation

(If you don't have a brain mask derived from the skullstripping pipeline, please set `-s TRUE`.)

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --t1 "*T1w*.nii.gz" --flair "*FLAIR*.nii.gz" -s TRUE --method dworcount --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --t1 "*T1w*.nii.gz" --flair "*FLAIR*.nii.gz" -s TRUE --method dworcount  --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --t1 "*T1w*.nii.gz" --flair "*FLAIR*.nii.gz" -s TRUE --method dworcount --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```

-   `docker`
```bash
docker pull pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --t1 "*T1w*.nii.gz" --flair "*FLAIR*.nii.gz" -s TRUE --method dworcount --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip 
```


**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.


### Count

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --step count --method dworcount --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --step count --method dworcount --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --step count --method dworcount --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```

-   `docker`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project -p sub-001 --ses ses-01 --step count --method dworcount --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip 
```

### Consolidation

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project --step consolidation --method dworcount -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/lesion_count/code/bash/lesion_count.sh -m /path/to/project --step consolidation --method dworcount -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```


