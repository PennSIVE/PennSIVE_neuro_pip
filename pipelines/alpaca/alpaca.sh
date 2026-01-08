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
  echo "  -s, --skullstripping   Specify whether to run skull stripping step. Default is FALSE"
  echo "  --mode   Specify whether to run the pipeline individually or in a batch: individual or batch. Default is batch"
  echo "  -c, --container   Specify the container to use: singularity, docker, local, cluster. Default is cluster"
  echo "  --sinpath   Specify the path to the singularity image if a singularity container is used. A default path is provided: /project/singularity_images/pennsive_amd64_cputorch.sif"
  echo "  --dockerpath   Specify the path to the docker image if a docker container is used. A default path is provided: russellshinohara/pennsive_amd64_cputorch"
  echo "  --hdbetpath   Specify the path to the HD-BET binary. A default path is provided: ~/.local/bin/hd-bet"
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
epi=""
skullstripping=FALSE
mode=batch
c=cluster
sin_path="/project/singularity_images/pennsive_amd64_cputorch.sif"
tool_path=""
docker_path=russellshinohara/pennsive_amd64_cputorch
hdbet_path="~/.local/bin/hd-bet"

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
    --hdbetpath)
      shift
      hdbet_path=$1
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
  if [ "$mode" == "batch" ]; then
    
    patient=`ls $main_path/data`

    # CVS Pipeline
    for p in $patient;
    do 
        ses=`ls $main_path/data/$p`
        for s in $ses;
        do
          t1_r=`find $main_path/data/$p/$s/anat -name $t1 -type f | xargs -I {} basename {}`
          flair_r=`find $main_path/data/$p/$s/anat -name $flair -type f | xargs -I {} basename {}`
          epimag_r=`find $main_path/data/$p/$s/anat -name $epi_mag -type f | xargs -I {} basename {}`
          epipha_r=`find $main_path/data/$p/$s/anat -name $epi_phase -type f | xargs -I {} basename {}`
          if [ "$c" == "cluster" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
          elif [ "$c" == "local" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
          elif [ "$c" == "singularity" ]; then
            module load apptainer
            bsub -J "alpaca" -oo $main_path/log/output/alpaca_output_${p}_${s}.log -eo $main_path/log/error/alpaca_error_${p}_${s}.log singularity run --cleanenv \
               -B $main_path \
               -B $tool_path \
               -B /scratch $sin_path \
               Rscript $tool_path/pipelines/alpaca//alpaca.R --mainpath $main_path \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --skullstripping $skullstripping --hdbetpath $hdbet_path \
            --lesioncenter $tool_path/lesioncenter --helpfunc $tool_path/help_functions
                elif [ "$c" == "docker" ]; then
            docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript $tool_path/pipelines/alpaca/alpaca.R --mainpath $main_path \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --skullstripping $skullstripping --hdbetpath $hdbet_path \
            --lesioncenter $tool_path/lesioncenter --helpfunc $tool_path/help_functions > /home/main/log/output/alpaca_output_${p}_${s}.log 2> /home/main/log/error/alpaca_error_${p}_${s}.log
          fi
        done
    done
  elif [ "$mode" == "individual" ]; then
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
    echo "session number is "
    echo $ses
    echo "hdbetpath is"
    echo $hdbetpath
    echo $hdbet_path
    t1_r=`find $main_path/data/$p/$ses/anat -name $t1 -type f | xargs -I {} basename {}`
    flair_r=`find $main_path/data/$p/$ses/anat -name $flair -type f | xargs -I {} basename {}`
    epimag_r=`find $main_path/data/$p/$ses/anat -name $epi_mag -type f | xargs -I {} basename {}`
    epipha_r=`find $main_path/data/$p/$ses/anat -name $epi_phase -type f | xargs -I {} basename {}`
    if [ "$c" == "cluster" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1          
            elif [ "$c" == "local" ]; then
                  echo "Error: ALPaCA is only available as a container."
                  show_help
                  exit 1         
            elif [ "$c" == "singularity" ]; then
             module load apptainer
              bsub -J "alpaca" -oo $main_path/log/output/alpaca_output_${p}_${ses}.log -eo $main_path/log/error/alpaca_error_${p}_${ses}.log singularity run --cleanenv \
               -B $main_path \
               -B $tool_path \
               -B /scratch $sin_path \
               Rscript $tool_path/pipelines/alpaca/alpaca.R --mainpath $main_path \
              --participant $p --session $ses --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --skullstripping $skullstripping --hdbetpath $hdbet_path \
              --lesioncenter $tool_path/lesioncenter --helpfunc $tool_path/help_functions
          elif [ "$c" == "docker" ]; then
            docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript $tool_path/pipelines/alpaca/alpaca.R --mainpath $main_path \
            --participant $p --session $ses --t1 $t1_r --flair $flair_r --epi_mag $epimag_r --epi_phase $epipha_r --skullstripping $skullstripping --hdbetpath $hdbet_path \
            --lesioncenter $tool_path/lesioncenter --helpfunc $tool_path/help_functions > /home/main/log/output/alpaca_output_${p}_${ses}.log 2> /home/main/log/error/alpaca_error_${p}_${ses}.log
          fi
  fi


