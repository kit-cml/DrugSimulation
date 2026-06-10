#!/bin/bash
source ../../../../../../venv/myenv/bin/activate

# Source the function script
source "../scripts/create_concs_directories.sh"
source "../scripts/zip_files.sh"
source "../scripts/progress_percentage.sh"

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}":/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo "${LD_LIBRARY_PATH}"
export PATH="${PATH}"
echo "${PATH}"

RESULT_FOLDER="./results"

# Clear any old PID file
PIDFILE="mpiexec.pid"
rm -f "${PIDFILE}"

NUMBER_OF_CPU="$(grep "^number_of_cpu" param.txt | cut -d'=' -f2 | sed 's/\/\/.*//' | xargs)"
MAX_CPU="$(nproc --all)"
if [[ "${NUMBER_OF_CPU}" -gt  "${MAX_CPU}" ]]; then
  echo "The processor requested exceed the number of maximum of CPU in this PC!" >> "${RESULT_FOLDER}/logfile" 2>&1
  exit 1
fi

# to grab cell_model value from parameter file (thanks, ChatGPT).
# grep "^user_name": looks for the line starting with user_name
# cut -d'=' -f2: gets the right-hand side of =
# sed 's/\/\/.*//': removes any inline comment starting with //
# xargs: trims leading and trailing whitespace
CELL_MODEL="$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
USER_NAME="$(grep "^user_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
DRUG_NAME="$(grep "^drug_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"
DRUG_CONCENTRATIONS="$(grep "^drug_concentrations" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)"

INITIAL_VALUES_SUBSTRING="*initial_values*.csv"
INITIAL_VALUES_ZIPNAME="${DRUG_NAME}_${CELL_MODEL}_initial_values.zip"
TIME_SERIES_SUBSTRING="*time_series*.csv"
TIME_SERIES_ZIPNAME="${DRUG_NAME}_${CELL_MODEL}_time_series.zip"
FEATURES_SUBSTRING="*features*.csv"
FEATURES_ZIPNAME="${DRUG_NAME}_${CELL_MODEL}_features.zip"
# Only uncomment these 2 lines when CMLDEBUG preprocessor is active
#LAST_PACES_SUBSTRING="*last_paces*.csv"
#LAST_PACES_ZIPNAME="${DRUG_NAME}_${CELL_MODEL}_last_paces.zip"

# Split the string into an array
IFS=',' read -r -a drug_concentrations <<< "${DRUG_CONCENTRATIONS}"

# choose the binary based on the value of cell_model
if [[ $CELL_MODEL == *"CiPAORdv1.0_Land"* ]]; then
  BINARY_FILE="drugsim_CiPAORdv1.0_Land"
elif [[ $CELL_MODEL == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE="drugsim_CiPAORdv1.0"
elif [[ $CELL_MODEL == *"ToR-ORd_Land"* ]]; then
  BINARY_FILE="drugsim_ToR-ORd_Land"
elif [[ $CELL_MODEL == *"ToR-ORd"* ]]; then
  BINARY_FILE="drugsim_ToR-ORd"
elif [[ $CELL_MODEL == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE="drugsim_ToR-ORd-dynCl"
elif [[ $CELL_MODEL == *"ORd-static_Land"* ]]; then
  BINARY_FILE="drugsim_ORd-static_Land"
elif [[ $CELL_MODEL == *"ORd-static"* ]]; then
  BINARY_FILE="drugsim_ORd-static"
elif [[ "${CELL_MODEL}" == *"Grandi"* ]]; then
  BINARY_FILE="drugsim_Grandi"
else
  echo "The cell model ${CELL_MODEL} is not specified to any simulations!!" >> "${RESULT_FOLDER}/logfile" 2>&1
  exit 1
fi

echo "Cleaning previous results and logs..." 2>&1
rm -rf "${RESULT_FOLDER}" "${PIDFILE}" "${DEPOL_FAILED_FILE}"
echo "Cleaning successful!"
create_drug_concentration_directories "${RESULT_FOLDER}" "${DRUG_NAME}" "${CELL_MODEL}" "${drug_concentrations[@]}"


# Logging starts from here
echo "DrugSimulation Log Start..." >& "${RESULT_FOLDER}/logfile" 2>&1
START_TIME=$(date +%s)
echo "Run ${CELL_MODEL} cell model simulation with ${NUMBER_OF_CPU} cores." >> "${RESULT_FOLDER}/logfile" 2>&1
echo "Start time: $(date)" >> "${RESULT_FOLDER}/logfile" 2>&1

# Setup the percentage watcher.
LOG_FILE="${RESULT_FOLDER}/logfile"
PERCENTAGE_FILE="${RESULT_FOLDER}/percentage.txt"
echo 0 > ${PERCENTAGE_FILE}
watch_progress_percentage "$LOG_FILE" "$PERCENTAGE_FILE" &
PROGRESS_WATCHER_PID=$!

# Execute the simulation
( echo $$ > "${PIDFILE}"; exec mpiexec -np "${NUMBER_OF_CPU}" ../${BINARY_FILE} -input_deck param.txt >> "${RESULT_FOLDER}/logfile" 2>&1 )

# Uncomment this for debugging, and comment the original one above.
#set -x
#exec mpiexec -np 1 gdb --args "../${BINARY_FILE}" -input_deck param.txt >> "${RESULT_FOLDER}/logfile" 2>&1
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Simulation program got some problems!!! Exiting..." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi

echo "Zipping folder..." >> "${RESULT_FOLDER}/logfile" 2>&1
emit_progress_percentage 75 "${LOG_FILE}"
zip_files "${RESULT_FOLDER}" "${INITIAL_VALUES_SUBSTRING}" "${INITIAL_VALUES_ZIPNAME}"
zip_files "${RESULT_FOLDER}" "${TIME_SERIES_SUBSTRING}" "${TIME_SERIES_ZIPNAME}"
zip_files "${RESULT_FOLDER}" "${FEATURES_SUBSTRING}" "${FEATURES_ZIPNAME}"
#zip_files "${RESULT_FOLDER}" "${LAST_PACES_SUBSTRING}" "${LAST_PACES_ZIPNAME}"
mv "${INITIAL_VALUES_ZIPNAME}" "${RESULT_FOLDER}/."
mv "${TIME_SERIES_ZIPNAME}" "${RESULT_FOLDER}/."
mv "${FEATURES_ZIPNAME}" "${RESULT_FOLDER}/."
#mv "${LAST_PACES_ZIPNAME}" "${RESULT_FOLDER}/."
echo "Zipping finished" >> "${RESULT_FOLDER}/logfile" 2>&1


echo "Generate Report using AI...." >> "${RESULT_FOLDER}/logfile" 2>&1
emit_progress_percentage 80 "${LOG_FILE}"
bash "./generate_report.sh" >> "${RESULT_FOLDER}/logfile" 2>&1
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Reporting program got some problems!!! Exiting..." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
fi
rm -rf "${PIDFILE}"
END_TIME=$(date +%s)
echo "End time: $(date)" >> "${RESULT_FOLDER}/logfile" 2>&1
ELAPSED_TIME=$(( ${END_TIME} - ${START_TIME} ))
ELAPSED_TIME_MINUTES=$(( ${ELAPSED_TIME} / 60 ))


emit_progress_percentage 100 "${LOG_FILE}"
echo "All process have finished and it took ${ELAPSED_TIME_MINUTES} minutes! Check the logfile for more details." >> "${RESULT_FOLDER}/logfile" 2>&1

kill "$PROGRESS_WATCHER_PID" 2>/dev/null
wait "$PROGRESS_WATCHER_PID" 2>/dev/null


exit 0;

deactivate
