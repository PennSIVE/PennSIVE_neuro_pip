#!/bin/bash

# Define a function to display the help message
show_help() {
  echo "Usage: cvs.sh [option]"
  echo "Options:"
  echo "  -h, --help    Show help message"
  echo "  -m, --mainpath    Look for files in the mainpath"
  echo "  -p, --participant    Specify the participant id"
  echo "  --ses    Specify the session id"
  echo "  -t, --t1    Specify the T1 sequence name"
  echo "  -f, --flair   Specify the FLAIR sequence name"
  echo "  -ema, --epimag   Specify the EPI magnitude image name"
  echo "  -eph, --epiphase   Specify the unwrapped EPI phase image name"
  echo "  -n, --n4   Specify whether to run bias correction step. Default is TRUE"
  echo "  -s, --skullstripping   Specify whether to run skull stripping step. Default is FALSE"
  echo "  -r, --registration   Specify whether to run registration step. Default is TRUE"
  echo "  -w, --whitestripe   Specify whether to run whitestripe step. Default is TRUE"
  echo "  --mimosa   Specify whether to run mimosa segmentation step. Default is TRUE"
  echo "  --threshold   Specify the threshold used to generate mimosa mask. Default is 0.05"
  echo "  --step   Specify the step of pipeline. estimation or consolidation. Default is estimation"
  echo "  --mode   Specify whether to run the pipeline individually or in a batch: individual or batch. Default is batch"
  echo "  -c, --container   Specify the container to use: singularity, docker, local, cluster. Default is cluster"
  echo "  --sinpath   Specify the path to the singularity image if a singularity container is used. A default path is provided: /project/singularity_images/pennsive_amd64_cputorch.sif"
  echo "  --dockerpath   Specify the path to the docker image if a docker container is used. Note - running on MacOS, be sure the memory limit of the VM (set in Docker Desktop) is sufficiently large. A default path is provided: russellshinohara/pennsive_arm64_cputorch:v1.1"
  echo "  --dockermem   Specify the memory and swap allocated to the docker image if a docker container is used. A default is provided: 48g"
  echo "  --hdbetpath   Specify the path to the HD-BET binary. For the Docker container, this is hard-coded in this version. If pre-downloading model, this should be in ~/hd-bet_params/release_2.0.0/. A default path is provided: ~/.local/bin/hd-bet"
  echo "  --toolpath   Specify the path to the saved pipeline folder, eg: /path/to/folder"
}

# Check if any argument is provided
if [ $# -eq 0 ]; then
  echo "Error: No arguments provided."
  show_help
  exit 1
fi

# Initialize variables
main_path=""
p=""
ses=""
t1=""
flair=""
ema=""
eph=""
n4=TRUE
skullstripping=FALSE
registration=TRUE
whitestripe=TRUE
mimosa=TRUE
threshold=0.05
step=estimation
mode=batch
c=cluster
sin_path="/project/singularity_images/pennsive_amd64_cputorch.sif"
tool_path=""
docker_path=russellshinohara/pennsive_arm64_cputorch:v1.1
hdbet_path="~/.local/bin/hd-bet"
docker_mem=48g
hdbet_model="~/hd-bet_params/"

# Parse command-line arguments
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)
      show_help
      exit 0
      ;;
    -m|--mainpath)
      shift
      main_path=$1
      ;;
    -p|--participant)
      shift
      p=$1
      ;;
    --ses)
      shift
      ses=$1
      ;;
    -t|--t1)
      shift
      t1=$1
      ;;
    -f|--flair)
      shift
      flair=$1
      ;;
    -ema|--epimag)
      shift
      epi_mag=$1
      ;;
    -eph|--epiphase)
      shift
      epi_phase=$1
      ;;
    -n|--n4)
      shift
      n4=$1
      ;;
    -s|--skullstripping)
      shift
      skullstripping=$1
      ;;
    -r|--registration)
      shift
      registration=$1
      ;;
    -w|--whitestripe)
      shift
      whitestripe=$1
      ;;
    --mimosa)
      shift
      mimosa=$1
      ;;
    --threshold)
      shift
      threshold=$1
      ;;
    --step)
      shift
      step=$1
      ;;
    --mode)
      shift
      mode=$1
      ;;
    -c|--container)
      shift
      c=$1
      ;;
    --sinpath)
      shift
      sin_path=$1
      ;;
    --dockerpath)
      shift
      docker_path=$1
      ;;
    --dockermem)
      shift
      docker_mem=$1
      ;;
    --hdbetpath)
      shift
      hdbet_path=$1
      ;;
    --hdbetmodel)
      shift
      hdbet_model=$1
      ;;
    --toolpath)
      shift
      tool_path=$1
      ;;
    *)
      echo "Error: Invalid option '$1'."
      show_help
      exit 1
      ;;
  esac
  shift
done

# Check if required options are provided
if [ -z "$main_path" ]; then
  echo "Error: Main path not specified."
  show_help
  exit 1
fi

mkdir -p $main_path/log/output
mkdir -p $main_path/log/error


if [ "$step" = "estimation" ]; then

  if [ -z "$t1" ]; then
    echo "Error: T1 MPRAGE not specified."
    show_help
    exit 1
  fi

  if [ -z "$flair" ]; then
    echo "Error: FLAIR not specified."
    show_help
    exit 1
  fi

  if [ -z "$epi_mag" ]; then
    echo "Error: EPI magnitude not specified."
    show_help
    exit 1
  fi

  if [ -z "$epi_phase" ]; then
    echo "Error: EPI phase not specified."
    show_help
    exit 1
  fi
  if [ "$mode" = "batch" ]; then


    echo "Available Patients dir: $main_path/data"

    for p in "$main_path/data"/*/; do
      
      [ -d "$p" ] || continue
      p="$(basename "${p%/}")"
      echo "Running Patient: $p"

      dir="$main_path/data/$p"
      echo "Available Sessions dir: $dir"

      for s in "$dir"/*/; do

        [ -d "$s" ] || continue
        s="$(basename "${s%/}")"
        echo "Running Session: $s"
          
          t1_r=`find $main_path/data/$p/$s/anat -name $t1 -type f | xargs -I {} basename {}`
          flair_r=`find $main_path/data/$p/$s/anat -name $flair -type f | xargs -I {} basename {}`
          epimag_r=`find $main_path/data/$p/$s/anat -name $epi_mag -type f | xargs -I {} basename {}`
          epipha_r=`find $main_path/data/$p/$s/anat -name $epi_phase -type f | xargs -I {} basename {}`
          if [ "$c" = "cluster" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
          elif [ "$c" = "local" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
          elif [ "$c" = "singularity" ]; then
            module load apptainer
            bsub -J "alpaca" -oo $main_path/log/output/alpaca_output_${p}_${s}.log -eo $main_path/log/error/alpaca_error_${p}_${s}.log singularity run --cleanenv \
               -B $main_path \
               -B $tool_path \
               -B /scratch $sin_path \
               Rscript $tool_path/pipelines/alpaca/alpaca_mod.R --mainpath $main_path \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --hdbetpath $hdbet_path --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
                elif [ "$c" = "docker" ]; then
            docker run --memory=$docker_mem --memory-swap=$docker_mem --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/alpaca/alpaca_mod.R --mainpath /home/main \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --hdbetpath /opt/fsl-6.0.7.19/bin/hd-bet --lesioncenter /home/tool/lesioncenter --mpath /home/tool/pipelines/mimosa/model/mimosa_model.RData --helpfunc /home/tool/help_functions > $main_path/log/output/alpaca_output_${p}_${s}.log 2> $main_path/log/error/alpaca_error_${p}_${s}.log
          fi
        done
    done
  elif [ "$mode" = "individual" ]; then
    if [ -z "$p" ]; then
      echo "Error: Participant id not provided for individual processing."
      show_help
      exit 1
    fi

    if [ -z "$ses" ]; then
      echo "Error: Session id not provided for individual processing."
      show_help
      exit 1
    fi
    t1_r=`find $main_path/data/$p/$ses/anat -name $t1 -type f | xargs -I {} basename {}`
    flair_r=`find $main_path/data/$p/$ses/anat -name $flair -type f | xargs -I {} basename {}`
    epimag_r=`find $main_path/data/$p/$ses/anat -name $epi_mag -type f | xargs -I {} basename {}`
    epipha_r=`find $main_path/data/$p/$ses/anat -name $epi_phase -type f | xargs -I {} basename {}`
    if [ "$c" = "cluster" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
            elif [ "$c" = "local" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1         
            elif [ "$c" = "singularity" ]; then
             module load apptainer
              bsub -J "alpaca" -oo $main_path/log/output/alpaca_output_${p}_${ses}.log -eo $main_path/log/error/alpaca_error_${p}_${ses}.log singularity run --cleanenv \
               -B $main_path \
               -B $tool_path \
               -B /scratch $sin_path \
               Rscript $tool_path/pipelines/alpaca/alpaca_mod.R --mainpath $main_path \
            --participant $p --session $ses --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --hdbetpath $hdbet_path --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
          elif [ "$c" = "docker" ]; then
            docker run --memory=$docker_mem --memory-swap=$docker_mem --rm -it -v $main_path:/home/main -v $tool_path:/home/tool -v $HOME/hd-bet_params:/home/$USER/hd-bet_params -e HOME=/home/$USER $docker_path Rscript /home/tool/pipelines/alpaca/alpaca_mod.R --mainpath /home/main \
            --participant $p --session $ses --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --hdbetpath /opt/fsl-6.0.7.19/bin/hd-bet --lesioncenter /home/tool/lesioncenter --mpath /home/tool/pipelines/mimosa/model/mimosa_model.RData --helpfunc /home/tool/help_functions > $main_path/log/output/alpaca_output_${p}_${ses}.log 2> $main_path/log/error/alpaca_error_${p}_${ses}.log
      fi
  fi
fi 

if [ "$step" = "consolidation" ]; then
    if [ "$c" = "cluster" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
    elif [ "$c" = "local" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1   
    elif [ "$c" = "singularity" ]; then
      module load apptainer
      bsub -J "alpaca" -oo $main_path/log/output/alpaca_output_consolidation.log -eo $main_path/log/error/alpaca_error_consolidation.log singularity run --cleanenv \
       -B $main_path \
       -B $tool_path \
       -B /scratch $sin_path \
       Rscript $tool_path/pipelines/alpaca/alpaca_mod.R  --mainpath $main_path --step $step
    elif [ "$c" = "docker" ]; then
      docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool -v $HOME/hd-bet_params:/home/$USER/hd-bet_params -e HOME=/home/$USER $docker_path Rscript /home/tool/pipelines/alpaca/alpaca_mod.R --mainpath /home/main --step $step > $main_path/log/output/alpaca_output_consolidation.log 2> $main_path/log/error/alpaca_error_consolidation.log
    fi
fi
  



