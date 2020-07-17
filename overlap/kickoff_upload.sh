#!/bin/bash

# Description: Print script usage message
# Arguments: None
# Returns: None
script_usage() {
  echo "
[USAGE]: Use this script to upload fusion detection results to the fusion \
detection database.

-h * [None] Print this help message.
-d * [results_dir1,results_dir2,...] List of directories (comma-separated) where the fusion results are located. Directories must be full paths. Ex: /path/to/fusion_results1,/path/to/fusion_results2.
-p * [patient1,patient2,...] List of patient IDs (comma-separated). The number of patient IDs must be equal to the number of results directories and they must be listed in the same respective order. Ex: patient1,patient2. (Also referred to as the 'subject ID'.)
-s * [sample1,sample2,...] List of sample IDs (comma-separated). The number of sample IDs must be equal to the number of results directories and they must be listed in the same respective order. Ex: sample1,sample2.
-r * [reads1,reads2,...] List of read pairs numbers (comma-separated). The number of read pairs numbers must be equal to the number of results directories and they must be listed in the same respective order. Optional. Ex: 40000000,52000000.

[Example 1]: kickoff_upload.sh -d /path/to/fusion_results1 -p p1 -s s1
[Example 2]: kickoff_upload.sh -d /path/to/fusion_results1,/path/to/fusion_results2 -p p1,p2 -s s1,s2
[Example 3]: kickoff_upload.sh -d /path/to/fusion_results1 -p p1 -s s1 -r 42543879
"
}

# Description: Process command line inputs
# Arguments: See script_usage for details
# Globals: sample_dirs_array, subject_names_array, sample_names_array, and
#          read_numbers_array
# Returns: None
get_command_line_input() {
  if [ "$#" = "0" ]; then
    script_usage
    exit 1
  fi
  
  while getopts ":d:p:s:r:h" opt; do
    case $opt in
      h)
        script_usage
        exit 0
        ;;
      d)
        local sample_dirs
        sample_dirs="$OPTARG"
        ;;
      p)
        local subject_names
        subject_names="$OPTARG"
        ;;
      s)
        local sample_names
        sample_names="$OPTARG"
        ;;
      r)
        local read_numbers
        read_numbers="$OPTARG"
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

  if [[ -z "$sample_dirs" ]]; then
    echo "[ERROR]: Please specify a results directory. See -h for details."
    exit 1
  fi
  sample_dirs_array=($(echo $sample_dirs | sed "s/,/ /g"))
  local num_sample_dirs=${#sample_dirs_array[@]}
  
  if [[ -z "$subject_names" ]]; then
    echo "[ERROR]: Please specify a patient ID. See -h for details."
    exit 1
  fi
  subject_names_array=($(echo $subject_names | sed "s/,/ /g"))
  local num_subject_names=${#subject_names_array[@]}

  if [[ $num_sample_dirs -ne $num_subject_names ]]; then
    echo "[ERROR]: The number of results directories ($num_sample_dirs) is not equal to the number of patient IDs ($num_subject_names). See -h for details."
    exit 1
  fi

  if [[ -z "$sample_names" ]]; then
    echo "[ERROR]: Please specify a sample ID. See -h for details."
    exit 1
  fi
  sample_names_array=($(echo $sample_names | sed "s/,/ /g"))
  local num_sample_names=${#sample_names_array[@]}

  if [[ $num_sample_dirs -ne $num_sample_names ]]; then
    echo "[ERROR]: The number of results directories ($num_sample_dirs) is not equal to the number of sample IDs ($num_sample_names). See -h for details."
    exit 1
  fi

  if [[ -z "$read_numbers" ]]; then
    for (( i=0; i<${num_sample_dirs}; i++ )); do
      read_numbers_array[$i]=""
    done
  else
    read_numbers_array=($(echo $read_numbers | sed "s/,/ /g"))
  fi
  local num_read_numbers=${#read_numbers_array[@]}

  if [[ $num_sample_dirs -ne $num_read_numbers ]]; then
    echo "[ERROR]: The number of results directories ($num_sample_dirs) is not equal to the number of read pairs ($num_read_numbers). See -h for details."
    exit 1
  fi
}

# Description: Run the upload script
# Arguments: None
# Globals: script_dir, sample_dirs_array, subject_names_array,
#          sample_names_array, and read_numbers_array
# Returns: None
run_upload() {
  for (( i=0; i<${#sample_dirs_array[@]}; i++ ))
  do
    sample_dir=${sample_dirs_array[$i]}
    subject_name=${subject_names_array[$i]}
    sample_name=${sample_names_array[$i]}
    read_number=${read_numbers_array[$i]}

    echo ""
    echo ""
    echo ""
    echo "Upload arguments -"
    echo "ResultsDir: $sample_dir"
    echo "Patient:    $subject_name"
    echo "Sample:     $sample_name"
    echo "Reads:      $read_number"
    echo ""
    read -p "Do you want to run the upload (Y/n)? " yn
    continue="y"
    case $yn in
      [Yy]* ) continue="y";;
      [Nn]* ) continue="";;
    esac
    echo ""

    if [[ ! -z "$continue" ]]; then
      cd $sample_dir
      if [[ -z "$read_number" ]]; then
        $script_dir/../R/upload_fusion_results.R --subject $subject_name --sample $sample_name
      else
        $script_dir/../R/upload_fusion_results.R --subject $subject_name --sample $sample_name --reads $read_number
      fi
    fi
  done
}

main() {
  script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  get_command_line_input "$@"
  run_upload
}

main "$@"
