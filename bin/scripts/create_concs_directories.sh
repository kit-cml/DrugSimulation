#!/bin/bash

create_drug_concentration_directories() {
    local RESULT_FOLDER="$1"
    local USER_NAME="$2"
    local DRUG_NAME="$3"
    local CELL_MODEL="$4"
    shift 4
    local drug_concentrations=("$@")

    local CONTROL_CONC=0.00

    # Create base user directory and drug directory
    mkdir -p "${RESULT_FOLDER}/${USER_NAME}/${DRUG_NAME}_${CELL_MODEL}"

    # Create control concentration directory
    mkdir -p "${RESULT_FOLDER}/${USER_NAME}/${DRUG_NAME}_${CELL_MODEL}/${CONTROL_CONC}"

    # Loop through each drug concentration and create directory
    for conc in "${drug_concentrations[@]}"; do
        dir_name=$(printf "%.2f" "$conc")
        mkdir -p "${RESULT_FOLDER}/${USER_NAME}/${DRUG_NAME}_${CELL_MODEL}/${dir_name}"
    done
}

