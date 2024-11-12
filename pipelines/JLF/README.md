# JLF Pipeline

The JLF pipeline produces a high-resolution anatomical segmentation using the ANTs Joint Label Fusion algorithm. A T1-weighted image is the only input needed. 

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide four types of scenarios for running the pipeline: `local` (running the pipeline locally), `cluster` (running the pipeline on High Performance Computing Cluster), `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the docker container). 

The pipeline contains three stages: 1) Registration: registers atlas into participants' T1-weighted space, 2) antsjointfusion: segments T1 images using multi-atlas segmentation with joint label fusion, 3) Extraction: extracts ROI and lesion volumes.

Detailed examples are provided below (all in individual mode):

### Registration

(If you don't have a brain mask derived from the skull-stripping pipeline, please set `-s TRUE`.)

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --step registration --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --step registration --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
singularity pull -F $sin_path docker://pennsive/neuror
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --step registration --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`

```bash
docker pull pennsive/neuror
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --step registration --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip 
```

**Note**: If you are using the `takim` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

### ANTs Joint Fusion

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --mode individual -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --mode individual -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project -p sub-001 --ses ses-01 -t1 "*T1w*.nii.gz" --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip
```

### Extraction

-   `local` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project --step extraction -c local --toolpath /path/to/PennSIVE_neuro_pip
```

-   `cluster`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project --step extraction -c cluster --toolpath /path/to/PennSIVE_neuro_pip
```

-   `singularity` 
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project --step extraction -c singularity --toolpath /path/to/PennSIVE_neuro_pip --sinpath $sin_path
```


-   `docker`
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/JLF/code/bash/JLF.sh -m /path/to/project --step extraction -c docker --toolpath /path/to/PennSIVE_neuro_pip 
```
