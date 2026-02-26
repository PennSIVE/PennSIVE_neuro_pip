# Automated Lesion, PRL, and CVS Analysis (ALPaCA) Pipeline

The ALPaCA pipeline integrates an automated technique for white matter lesion, PRL, and CVS segmentation, developed by [Dr. Hu](https://direct.mit.edu/imag/article/doi/10.1162/IMAG.a.932/133343/Automated-segmentation-of-multiple-sclerosis). It provides processed T1-weighted, T2-FLAIR, T2*-magnitude, and T2*-phase images, as well lesion, CVS, and PRL masks and lesion-level probabilities.

## Diagram
![MIMoSA Workflow](/pipelines/alpaca/figure/alpaca_pipeline.png)

## Data Structure
This pipeline requires all images to be organized in BIDS format. An example is provided below:

![Data Structure](/pipelines/alpaca/figure/data_structure.png)

## Pipeline Options
We offer two modes for the pipeline: individual and batch. If users want to run the pipeline for different participants one at a time, the participant's ID and session ID should be specified. Users can also run the pipeline in batch mode. The participant's ID and session ID can be skipped. Additionally, we provide two types of scenarios for running the pipeline: `singularity` (running the pipeline on High Performance Computing Cluster using the Singularity container), and `docker` (running the pipeline locally using the Docker container). Please pull the container below corresponding to your system requirements:

- **AMD, CPU**: russellshinohara/pennsive_amd64_cputorch:v1.1
- **AMD, GPU**: russellshinohara/pennsive_amd64_gputorch:v1.1
- **ARM, CPU**: russellshinohara/pennsive_arm64_cputorch:v1.1

Detailed examples for running the pipeline are provided below (all in individual mode):

### Estimation

(If you don't have a brain mask derived from the skullstripping pipeline, please set `-s TRUE`. By default, images will be skullstripped with hd-bet.)


-   `singularity` 
```bash
singularity pull -F $sin_path docker://russellshinohara/pennsive_amd64_gputorch:v1.1
```
```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/alpaca/alpaca.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -ema "*_part-mag_T2star.nii.gz" -eph "*_part-phase_T2star_UNWRAPPED.nii.gz" -s TRUE --mode individual -c singularity --toolpath /path/to/PennSIVE_neuro_pip
```


-   `docker`

```bash
docker pull russellshinohara/pennsive_arm64_cputorch:v1.1
```

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/alpaca/alpaca.sh -m /path/to/project -p sub-001 --ses ses-01 -t "*_T1w.nii.gz" -f "*_FLAIR.nii.gz" -ema "*_part-mag_T2star.nii.gz" -eph "*_part-phase_T2star_UNWRAPPED.nii.gz" -s TRUE --mode individual -c docker --toolpath /path/to/PennSIVE_neuro_pip
```

**Note**: If you are using the `takim2` cluster within the PennSIVE group, you do not need to specify `sinpath`, which has been given a default path.

### Consolidation

-   `singularity` 

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/alpaca/alpaca.sh -m /path/to/project -c singularity --toolpath /path/to/PennSIVE_neuro_pip
```


-   `docker`

```bash
bash /path/to/PennSIVE_neuro_pip/pipelines/alpaca/alpaca.sh -m /path/to/project -c docker --toolpath /path/to/PennSIVE_neuro_pip
```

## Output Data Structure
![Output](/pipelines/alpaca/figure/output.png)

Lesion-level probabilities are available in "probabilities.csv". Lesions in "alpaca_mask.nii.gz" are labeled as follows:

- 1: PRL- CVS-
- 3: PRL+ CVS-
- 5: PRL- CVS+
- 7: PRL+ CVS+