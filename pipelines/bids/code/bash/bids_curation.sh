#!/bin/bash

# Define a function to display the help message
show_help() {
  echo "Usage: bids_curation.sh [option]"
  echo "Options:"
  echo "  -h, --help    Show help message"
  echo "  -m, --mainpath    Specify the mainpath of DICOM data"
  echo "  -p, --participant    Specify the name of the participant"
  echo "  -s, --session    Specify the name of the session"
  echo "  -c, --container    Specify the type of container to use, docker or singularity. Default is docker"
  echo "  --sinpath    Specify the path to the saved singularity image if singularity container is being used. A default path is provided: /project/singularity_images/heudiconv_latest.sif"
  echo "  --dockerpath   Specify the path to the docker image if a docker container is used. A default path is provided: nipy/heudiconv"
  echo "  --toolpath    Specify the path to the saved help function scripts"
  echo "  --step   Specify the step of heudiconv. heuristic, customization, or bids. Default is heuristic"
  echo "  --mode   Specify whether to process data individually or through a batch, individual or batch. Default is individual"
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
s=""
mode="individual"
container="docker"
sin_path="/project/singularity_images/heudiconv_latest.sif"
docker_path="nipy/heudiconv"
tool_path=""
step="heuristic"


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
    -s|--session)
      shift
      s=$1
      ;;
    -c|--container)
      shift
      container=$1
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
    --mode)
      shift
      mode=$1
      ;;
    --step)
      shift
      step=$1
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

if [ "$step" == "customization" ]; then
    if [ -z "$tool_path" ]; then
      echo "Error: The path to the help function scripts is not provided."
      show_help
      exit 1
    fi
    Rscript $tool_path/pipelines/bids/code/R/heuristic_customize.R -m $main_path
fi

if [ "$mode" == "individual" ]; then
    
    if [ -z "$p" ]; then
        echo "Error: Participant's name not provided."
        show_help
        exit 1
    fi

    if [ -z "$s" ]; then
        echo "Error: Session name not provided."
        show_help
        exit 1
    fi

    if [ "$container" == "docker" ]; then
        if [ "$step" == "heuristic" ]; then
            docker run --rm -it -v $main_path:/home $docker_path -s $p -ss $s -d /home/original_data/{subject}/{session}/*/*.dcm -o /home/bids/ -f convertall -c none 
            mkdir -p $main_path/heuristic_script/$p/$s
            mkdir $main_path/heuristic_script/template
            mkdir -p $main_path/dicominfo/$p/$s
            cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/template
            cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/$p/$s
            dic_info=`find $main_path/bids/.heudiconv/$p -name dicominfo_ses-$s.tsv -type f`
            cp $dic_info $main_path/dicominfo/$p/$s
        fi

        if [ "$step" == "bids" ]; then
            rm -r $main_path/bids/.heudiconv
            docker run --rm -it -v $main_path:/home $docker_path -s $p -ss $s -d /home/original_data/{subject}/{session}/*/*.dcm -o /home/bids/ -f /home/heuristic_script/$p/$s/heuristic.py -c dcm2niix -b --overwrite
        fi

    fi 

    if [ "$container" == "singularity" ]; then
        module load singularity
        if [ -z "$sin_path" ]; then
            echo "Error: Path to the singularity image not provided."
            show_help
            exit 1
        fi
        
        if [ "$step" == "heuristic" ]; then
            bsub -J "heuristic" <<EOF 
            singularity run --cleanenv \
            -B $main_path $sin_path -s $p -ss $s -d $main_path/original_data/{subject}/{session}/*/*.dcm -o $main_path/bids/ -f convertall -c none && {
            echo "heudiconv process is done!" 
            mkdir -p $main_path/heuristic_script/$p/$s 
            mkdir $main_path/heuristic_script/template 
            mkdir -p $main_path/dicominfo/$p/$s 
            cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/template  
            cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/$p/$s
            cp $main_path/bids/.heudiconv/$p/info/dicominfo_ses-$s.tsv $main_path/dicominfo/$p/$s
            }
EOF
        fi

        if [ "$step" == "bids" ]; then
            rm -r $main_path/bids/.heudiconv
            bsub -J "bids_curation" singularity run --cleanenv -B $main_path $sin_path -s $p -ss $s -d $main_path/original_data/{subject}/{session}/*/*.dcm -o $main_path/bids/ -f $main_path/heuristic_script/$p/$s/heuristic.py -c dcm2niix -b --overwrite
        fi

    fi

fi

if [ "$mode" == "batch" ]; then
    patient=`ls $main_path/original_data | sed 's/\/$//'`
    mkdir -p $main_path/heuristic_script/template
    cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/template
    for p in $patient;
    do
        session=`ls $main_path/original_data/$p | grep -v '\.zip$' | sed 's/\/$//'`
        for s in $session;
        do
          if [ "$container" == "docker" ]; then
            if [ "$step" == "heuristic" ]; then
              docker run --rm -it -v $main_path:/home $docker_path -s $p -ss $s -d /home/original_data/{subject}/{session}/*/*.dcm -o /home/bids/ -f convertall -c none 
              mkdir -p $main_path/heuristic_script/$p/$s
              mkdir -p $main_path/dicominfo/$p/$s
              cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/$p/$s/
              dic_info=`find $main_path/bids/.heudiconv/$p -name dicominfo_ses-$s.tsv -type f`
              cp $dic_info $main_path/dicominfo/$p/$s/
            fi

            if [ "$step" == "bids" ]; then
              rm -r $main_path/bids/.heudiconv
              docker run --rm -it -v $main_path:/home $docker_path -s $p -ss $s -d /home/original_data/{subject}/{session}/*/*.dcm -o /home/bids/ -f /home/heuristic_script/$p/$s/heuristic.py -c dcm2niix -b --overwrite
            fi
          fi

          if [ "$container" == "singularity" ]; then
              module load singularity
              if [ "$step" == "heuristic" ]; then
                  bsub -J "heuristic" <<EOF 
                  singularity run --cleanenv \
                  -B $main_path $sin_path -s $p -ss $s -d $main_path/original_data/{subject}/{session}/*/*.dcm -o $main_path/bids/ -f convertall -c none && {
                  mkdir -p $main_path/heuristic_script/$p/$s
                  mkdir $main_path/heuristic_script/template
                  mkdir -p $main_path/dicominfo/$p/$s
                  cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/template
                  cp $tool_path/pipelines/bids/code/python/heuristic.py $main_path/heuristic_script/$p/$s
                  cp $main_path/bids/.heudiconv/$p/info/dicominfo_ses-$s.tsv $main_path/dicominfo/$p/$s
                  }
EOF
              fi
              if [ "$step" == "bids" ]; then
                  rm -r $main_path/bids/.heudiconv
                  bsub -J "bids_curation" singularity run --cleanenv -B $main_path $sin_path -s $p -ss $s -d $main_path/original_data/{subject}/{session}/*/*.dcm -o $main_path/bids/ -f $main_path/heuristic_script/$p/$s/heuristic.py -c dcm2niix -b --overwrite
              fi
          fi 
        done
    done
fi

