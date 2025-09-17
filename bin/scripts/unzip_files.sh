#!/bin/bash

unzip_files() {
  local ZIP_FILE="$1"
  local UNZIP_DESTINATION="$2"
  local RESULTS_ROOT="$3"
  local SAMPLE_SIZE="$4"
  shift 4
  local drug_concentrations=("$@")
   
  unzip "${ZIP_FILE}" -d "${UNZIP_DESTINATION}"
  EXIT_CODE=$?
  if [ $EXIT_CODE -ne 0 ]; then
    echo "Cannot unzip ${ZIP_FILE}! Please check whether the file is existing and healthy!!"
    exit 1
  fi
  echo "Zip file ${ZIP_FILE} has been successfully extracted into ${UNZIP_DESTINATION}"

  # Multiple checks of the unzipped folders
  echo "Validating unzipped files..."
  if [[ ! -d "${RESULTS_ROOT}" ]]; then
    echo "Missing root folder: ${RESULTS_ROOT}"
    exit 1
  fi
  echo "Root folder ${RESULTS_ROOT} found."

  STATUS=0
  # Loop through each drug concentration and check the contents
  for conc in "${drug_concentrations[@]}"; do
    FOLDER=$(printf "%.2f" "${conc}")
    TARGET_DIR="${RESULTS_ROOT}/${FOLDER}"
    if [[ ! -d "${TARGET_DIR}" ]]; then
      echo "Missing concentration folder: ${TARGET_DIR}"
      STATUS=1
      continue
    fi

    ACTUAL_COUNT=$(find "${TARGET_DIR}" -type f | wc -l)
    if [[ "${ACTUAL_COUNT}" -ne "${SAMPLE_SIZE}" ]]; then
      echo "File count mismatch in ${TARGET_DIR}: expected ${SAMPLE_SIZE}, found ${ACTUAL_COUNT}"
      STATUS=1
    else
      echo "${TARGET_DIR}: ${ACTUAL_COUNT} files OK"
    fi   
  done

  if [[ "${STATUS}" -eq 0 ]]; then
    echo "All structure and file checks passed!"
  else
    echo "Some folder or file checks failed."
  fi
  
  return "${STATUS}"
}
