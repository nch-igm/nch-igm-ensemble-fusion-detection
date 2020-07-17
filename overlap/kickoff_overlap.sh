#!/bin/bash

# Description: Print script usage message
# Arguments: None
# Returns: None
script_usage() {
  echo "
[USAGE]: This script kicks off an Rscript that overlaps, prioritizes, and \
filters fusion detection results from individual fuison detectors. Fusion \
detection algorithms must be run first. Please see README.md for more information.

-h * [None] Print this help message.
-p * [path1,path2,...] full path where the overlap \
script will be kicked off.
-s * [samples_file1,samples_file2,...] Files under top_dir1,top_dir2,... that \
contain the lists of samples to kickoff scripts for. The number of samples \
files provided must be equal to the number of top directories, and they must \
be listed in the same respective order. [Default: samples,samples,...].
"
}

# Description: Process command line inputs
# Arguments: See script_usage for details
# Globals: path_array, samples_files_array
# Returns: None
get_command_line_input() {
  if [ "$#" = "0" ]; then
    script_usage
    exit 1
  fi

  while getopts ":p:s:h" opt; do
    case $opt in
      h)
        script_usage
        exit 0
        ;;
      p)
        local path
        path="$OPTARG"
        ;;
      s)
        local samples_files
        samples_files="$OPTARG"
        ;;
      :)
        echo "[ERROR]: Flag -$OPTARG requires an argument. See -h for details."
        exit 1
        ;;
      ?)
        echo "[ERROR]: Invalid option -$OPTARG. See -h for valid options."
        exit 1
        ;;
    esac
  done

  if [[ -z "$path" ]]; then
    echo "[ERROR]: Please specify a path. See -h for details."
    exit 1
  fi
  path_array=($(echo $path | sed "s/,/ /g"))
  local num_path=${#path_array[@]}

  if [[ -z "$samples_files" ]]; then
    for (( i=0; i<${num_path}; i++ )); do
      samples_files_array[$i]=samples
    done
  else
    samples_files_array=($(echo $samples_files | sed "s/,/ /g"))
  fi
  local num_samples_files
  num_samples_files=${#samples_files_array[@]}

  if [[ $num_path -ne $num_samples_files ]]; then
    echo "[ERROR]: The number of top directories ($num_path) is not equal to"
    echo "the number of samples files ($num_samples_files). See -h for details."
    exit 1
  fi
}

# Description: Run the overlap script
# Arguments: None
# Globals: script_dir, path_array, samples_files_array
# Returns: None
run_overlap() {
  for (( i=0; i<${#path_array[@]}; i++ ))
  do
    path=${path_array[$i]}
    samples_file=${samples_files_array[$i]}
    projectdir=$path
    for sample in $(cat $projectdir/$samples_file)
    do
      echo $projectdir/$sample
      cd $projectdir/$sample
      $script_dir/../R/assemble_results.R --sample $sample
    done
  done
}

main() {
  script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  get_command_line_input "$@"
  run_overlap
}

main "$@"
