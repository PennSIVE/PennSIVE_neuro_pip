#!/bin/bash

# Define a function to display the help message
show_help() {
  echo "Usage: lesion_count.sh [option]"
  echo "Options:"
  echo "  -h, --help    Show help message"
  echo "  -m, --mainpath    Look for files in the mainpath"
  echo "  -p, --participant    Specify the participant id"
  echo "  --ses    Specify the session id"
  echo "  -t, --t1    Specify the T1 sequence name"
  echo "  -f, --flair   Specify the FLAIR sequence name"
  echo "  -n, --n4   Specify whether to run bias correction step. Default is TRUE"
  echo "  -s, --skullstripping   Specify whether to run skull stripping step. Default is FALSE"
  echo "  -r, --registration   Specify whether to run registration step. Default is TRUE"
  echo "  -w, --whitestripe   Specify whether to run whitestripe step. Default is TRUE"
  echo "  --mimosa   Specify whether to run mimosa segmentation step. Default is TRUE"
  echo "  --threshold   Specify the threshold used to generate mimosa mask. Default is 0.2"
  echo "  --step   Specify the step of pipeline. preparation, count, or consolidation. Default is preparation"
  echo "  --method   Specify the method of lesion counting. cc, dworcount, or both. Default is both"
  echo "  --mode   Specify whether to run the pipeline individually or in a batch: individual or batch. Default is batch"
  echo "  -c, --container   Specify the container to use: singularity, docker, local, cluster. Default is cluster"
  echo "  --sinpath   Specify the path to the singularity image if a singularity container is used. A default path is provided: /project/singularity_images/neuror_latest.sif"
  echo "  --dockerpath   Specify the path to the docker image if a docker container is used. A default path is provided: pennsive/neuror"
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
n4=TRUE
skullstripping=FALSE
registration=TRUE
whitestripe=TRUE
mimosa=TRUE
threshold=0.2
step=preparation
method=both
mode=batch
c=cluster
sin_path="/project/singularity_images/neuror_latest.sif"
tool_path=""
docker_path=pennsive/neuror

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
    --method)
      shift
      method=$1
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


if [ "$step" == "preparation" ]; then
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

  if [ "$mode" == "batch" ]; then

    participant=`ls $main_path/data`

    # Lesion Count Pipeline
    for p in $participant;
    do 
        ses=`ls $main_path/data/$p`
        for s in $ses;
        do
          t1_r=`find $main_path/data/$p/$s/anat -name $t1 \( -type f -o -type l \) | xargs -I {} basename {}`
          flair_r=`find $main_path/data/$p/$s/anat -name $flair \( -type f -o -type l \) | xargs -I {} basename {}`
          if [ "$c" == "cluster" ]; then
            bsub -oo $main_path/log/output/lesion_count_output_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_${p}_${s}.log \
            Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
          elif [ "$c" == "local" ]; then
            Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
            --participant $p --session $s --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions > $main_path/log/output/lesion_count_output_${p}_${s}.log 2> $main_path/log/error/lesion_count_error_${p}_${s}.log
          elif [ "$c" == "singularity" ]; then
            module load singularity
            bsub -J "lesion_count" -oo $main_path/log/output/lesion_count_output_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_${p}_${s}.log singularity run --cleanenv \
                -B $main_path \
                -B $tool_path \
                -B /scratch $sin_path \
                Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
                    --participant $p --session $s --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
                    --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
                    --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
          elif [ "$c" == "docker" ]; then
            docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/lesion_count/code/R/lesion_count.R --mainpath /home/main \
                    --participant $p --session $s --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
                    --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
                    --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath /home/tool/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions > /home/main/log/output/lesion_count_output_${p}_${s}.log 2> /home/main/log/error/lesion_count_error_${p}_${s}.log
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

    t1_r=`find $main_path/data/$p/$ses/anat -name $t1 \( -type f -o -type l \) | xargs -I {} basename {}`
    flair_r=`find $main_path/data/$p/$ses/anat -name $flair \( -type f -o -type l \) | xargs -I {} basename {}`
    
    if [ "$c" == "cluster" ]; then
      bsub -oo $main_path/log/output/lesion_count_output_${p}_${ses}.log -eo $main_path/log/error/lesion_count_error_${p}_${ses}.log \
      Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
      --participant $p --session $ses --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
      --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
      --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
    elif [ "$c" == "local" ]; then
            Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
            --participant $p --session $ses --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
            --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
            --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions > $main_path/log/output/lesion_count_output_${p}_${ses}.log 2> $main_path/log/error/lesion_count_error_${p}_${ses}.log
    elif [ "$c" == "singularity" ]; then
      module load singularity
      bsub -J "lesion_count" -oo $main_path/log/output/lesion_count_output_${p}_${ses}.log -eo $main_path/log/error/lesion_count_error_${p}_${ses}.log singularity run --cleanenv \
          -B $main_path \
          -B $tool_path \
          -B /scratch $sin_path \
          Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
              --participant $p --session $ses --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
              --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
              --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath $tool_path/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions
    elif [ "$c" == "docker" ]; then
        docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/lesion_count/code/R/lesion_count.R --mainpath /home/main \
                    --participant $p --session $ses --t1 $t1_r --flair $flair_r --n4 $n4 --skullstripping $skullstripping \
                    --registration $registration --whitestripe $whitestripe --mimosa $mimosa --threshold $threshold \
                    --step $step --method $method --lesioncenter $tool_path/lesioncenter --mpath /home/tool/pipelines/mimosa/model/mimosa_model.RData --helpfunc $tool_path/help_functions > /home/main/log/output/lesion_count_output_${p}_${ses}.log 2> /home/main/log/error/lesion_count_error_${p}_${ses}.log
    fi

  fi

elif [ "$step" == "count" ]; then
  if [ "$mode" == "batch" ]; then

      participant=`ls $main_path/data`

      # Lesion Count Pipeline
      for p in $participant;
      do 
          ses=`ls $main_path/data/$p`
          for s in $ses;
          do
            if [ "$c" == "cluster" ]; then
              bsub -oo $main_path/log/output/lesion_count_output_count_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_count_${p}_${s}.log \
              Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
              --participant $p --session $s --step $step --method $method 
            elif [ "$c" == "local" ]; then
              Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
              --participant $p --session $s --step $step --method $method > $main_path/log/output/lesion_count_output_count_${p}_${s}.log 2> $main_path/log/error/lesion_count_error_count_${p}_${s}.log
            elif [ "$c" == "singularity" ]; then
              module load singularity
              bsub -J "lesion_count" -oo $main_path/log/output/lesion_count_output_count_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_count_${p}_${s}.log singularity run --cleanenv \
                  -B $main_path \
                  -B $tool_path \
                  -B /scratch $sin_path \
                  Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
                      --participant $p --session $s --step $step --method $method
            elif [ "$c" == "docker" ]; then
              docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/lesion_count/code/R/lesion_count.R --mainpath /home/main \
                      --participant $p --session $s --step $step --method $method > /home/main/log/output/lesion_count_output_count_${p}_${s}.log 2> /home/main/log/error/lesion_count_error_count_${p}_${s}.log
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

    
    if [ "$c" == "cluster" ]; then
      bsub -oo $main_path/log/output/lesion_count_output_count_${p}_${ses}.log -eo $main_path/log/error/lesion_count_error_count_${p}_${ses}.log \
      Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
      --participant $p --session $ses --step $step --method $method 
    elif [ "$c" == "local" ]; then
            Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
            --participant $p --session $ses --step $step --method $method > $main_path/log/output/lesion_count_output_count_${p}_${ses}.log 2> $main_path/log/error/lesion_count_error_count_${p}_${ses}.log
    elif [ "$c" == "singularity" ]; then
      module load singularity
      bsub -J "lesion_count" -oo $main_path/log/output/lesion_count_output_count_${p}_${ses}.log -eo $main_path/log/error/lesion_count_error_count_${p}_${ses}.log singularity run --cleanenv \
          -B $main_path \
          -B $tool_path \
          -B /scratch $sin_path \
          Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
              --participant $p --session $ses --step $step --method $method 
    elif [ "$c" == "docker" ]; then
        docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/lesion_count/code/R/lesion_count.R --mainpath /home/main \
                    --participant $p --session $ses --step $step --method $method > /home/main/log/output/lesion_count_output_count_${p}_${ses}.log 2> /home/main/log/error/lesion_count_error_count_${p}_${ses}.log
    fi

  fi


elif [ "$step" == "consolidation" ]; then
  if [ "$c" == "cluster" ]; then
    bsub -oo $main_path/log/output/lesion_count_output_consolidation_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_consolidation_${p}_${s}.log \
    Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path --step $step --method $method 
  elif [ "$c" == "local" ]; then
    Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path \
    --step $step --method $method > $main_path/log/output/lesion_count_output_consolidation_${p}_${s}.log 2> $main_path/log/error/lesion_count_error_consolidation_${p}_${s}.log
  elif [ "$c" == "singularity" ]; then
    module load singularity
    bsub -J "lesion_count" -oo $main_path/log/output/lesion_count_output_consolidation_${p}_${s}.log -eo $main_path/log/error/lesion_count_error_consolidation_${p}_${s}.log singularity run --cleanenv \
        -B $main_path \
        -B $tool_path \
        -B /scratch $sin_path \
        Rscript $tool_path/pipelines/lesion_count/code/R/lesion_count.R --mainpath $main_path --step $step --method $method
  elif [ "$c" == "docker" ]; then
    docker run --rm -it -v $main_path:/home/main -v $tool_path:/home/tool $docker_path Rscript /home/tool/pipelines/lesion_count/code/R/lesion_count.R --mainpath /home/main \
            --step $step --method $method > /home/main/log/output/lesion_count_output_consolidation_${p}_${s}.log 2> /home/main/log/error/lesion_count_error_consolidation_${p}_${s}.log
  fi

fi


