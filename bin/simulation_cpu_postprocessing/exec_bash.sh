#!/bin/bash

# Source the function script
source "../scripts/create_concs_directories.sh"
source "../scripts/unzip_files.sh"
source "../scripts/zip_files.sh"

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64"
echo "${LD_LIBRARY_PATH}"
export PATH="${PATH}"
echo "${PATH}"

NUMBER_OF_CPU="$(grep "^number_of_cpu" param.txt | cut -d'=' -f2 | sed 's/\/\/.*//' | xargs)"
MAX_CPU="$(nproc --all)"
if [ "${NUMBER_OF_CPU}" -gt  "${MAX_CPU}" ]; then
  echo "The processor requested exceed the number of maximum of CPU in this PC!"
  exit 1
fi

# to grab cell_model value from parameter file (thanks, ChatGPT).
# grep "^user_name": looks for the line starting with user_name
# cut -d'=' -f2: gets the right-hand side of =
# sed 's/\/\/.*//': removes any inline comment starting with //
# xargs: trims leading and trailing whitespace
CELL_MODEL="$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"

RESULT_FOLDER="./results"
PLOT_FOLDER="./plots"
USER_NAME="$(grep "^user_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
DRUG_NAME="$(grep "^drug_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
INITIAL_VALUES_ZIP_FILE="$(grep "^initial_values_zip_file" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
DRUG_CONCENTRATIONS="$(grep "^drug_concentrations" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
HILL_FILE="$(grep "^hill_file" param.txt | cut -d '=' -f2 | sed 's|//.*||' | xargs)"
echo "FILE HILL: ${HILL_FILE}"
SAMPLE_SIZE=$(( $(wc -l < "${HILL_FILE}") - 1 ))
echo "FILE LINE: ${SAMPLE_SIZE}"



FEATURES_SUBSTRING="*features*.csv"
FEATURES_ZIPNAME="./${DRUG_NAME}_${CELL_MODEL}_${USER_NAME}_features.zip"

# Split the string into an array
IFS=',' read -r -a "${drug_concentrations}" <<< "${DRUG_CONCENTRATIONS}"

# choose the binary based on the value of cell_model
if [[ "${CELL_MODEL}" == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE="../drugsim_CiPAORdv1.0_postprocessing"
elif [[ "${CELL_MODEL}" == *"ORd-static"* ]]; then
  BINARY_FILE="../drugsim_ORd-static_postprocessing"
elif [[ "${CELL_MODEL}" == *"ToR-ORd"* ]]; then
  BINARY_FILE="../drugsim_ToR-ORd_postprocessing"
elif [[ "${CELL_MODEL}" == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE="../drugsim_ToR-ORd-dynCl_postprocessing"
else
  echo "The cell model ${CELL_MODEL} is not specified to any simulations!!"
  exit 1
fi

EXISTING_TIME_SERIES_PLOT_FOLDER="${PLOT_FOLDER}/time_series/${DRUG_NAME}_${CELL_MODEL}"
EXISTING_FEATURES_PLOT_FOLDER="${PLOT_FOLDER}/features/${DRUG_NAME}_${CELL_MODEL}"
EXISTING_RESULT_FOLDER="${RESULT_FOLDER}/${DRUG_NAME}_${CELL_MODEL}"
EXISTING_REPORT_FILES="report_drug_${DRUG_NAME}_${CELL_MODEL}_${USER_NAME}*"
echo "Cleaning ${EXISTING_RESULT_FOLDER} folder, ${EXISTING_TIME_SERIES_PLOT_FOLDER}, ${EXISTING_FEATURES_PLOT_FOLDER} and ${EXISTING_REPORT_FILES} files..."
rm -rf "${EXISTING_TIME_SERIES_PLOT_FOLDER}" "${EXISTING_FEATURES_PLOT_FOLDER}" "${EXISTING_RESULT_FOLDER}" "${EXISTING_REPORT_FILES}" "${RESULT_FOLDER}/logfile" "*.zip"
echo "Cleaning successful!"
create_drug_concentration_directories "$RESULT_FOLDER" "$DRUG_NAME" "$CELL_MODEL" "${drug_concentrations[@]}"
unzip_files "${INITIAL_VALUES_ZIP_FILE}" "./" "${RESULT_FOLDER}" "${DRUG_NAME}_${CELL_MODEL}" "${SAMPLE_SIZE}" "${drug_concentrations[@]}"
echo "Run ${CELL_MODEL} cell model simulation with ${NUMBER_OF_CPU} cores."
( echo $$ > "${PIDFILE}"; exec mpiexec -np "${NUMBER_OF_CPU}" "${BINARY_FILE}" "-input_deck" "param.txt" >& "${RESULT_FOLDER}/logfile")
zip_files "${RESULT_FOLDER}" "${FEATURES_SUBSTRING}" "${FEATURES_ZIPNAME}"
echo "Simulation has finished! Check the logfile for more details."
