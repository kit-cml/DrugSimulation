#!/bin/bash

create_drug_concentration_directories() {
  local RESULT_FOLDER="$1"
  local DRUG_NAME="$2"
  local CELL_MODEL="$3"
  shift 3
  local drug_concentrations=("$@")

  local CONTROL_CONC=0.00

  # Create control concentration directory
  mkdir -p "${RESULT_FOLDER}/${CONTROL_CONC}"

  # Loop through each drug concentration and create directory
  for conc in "${drug_concentrations[@]}"; do
    dir_name=$(printf "%.2f" "$conc")
    mkdir -p "${RESULT_FOLDER}/${dir_name}"
  done
}

