#PBS -N drugsim_postprocessing_pbs_job
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20000:00:00
#PBS -e stderr.log
#PBS -o stdout.log
#Specific the shell types
#PBS -S /bin/bash
#Specific the queue type
#PBS -q dque

cd "${PBS_O_WORKDIR}"
NPROCS="$(wc -l < $PBS_NODEFILE)"
echo "This job has allocated ${NPROCS} nodes"

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64"
echo "${LD_LIBRARY_PATH}"
export PATH="${PATH}"
echo "${PATH}"

# Source the function script
# Need to be invoked after the
# cd $PBS_O_WORKDIR command
source "../scripts/create_concs_directories.sh"
source "../scripts/unzip_files.sh"
source "../scripts/zip_files.sh"

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

TIME_SERIES_SUBSTRING="*time_series*.csv"
TIME_SERIES_ZIPNAME="${DRUG_NAME}_${CELL_MODEL}_time_series.zip"
FEATURES_SUBSTRING="*features*.csv"
FEATURES_ZIPNAME="./${DRUG_NAME}_${CELL_MODEL}_${USER_NAME}_features.zip"

# Split the string into an array
IFS=',' read -r -a "${drug_concentrations}" <<< "${DRUG_CONCENTRATIONS}"

# Logging starts from here
RESULT_FOLDER="./results"
echo "Cleaning ${RESULT_FOLDER}"
rm -rf "${RESULT_FOLDER}"
echo "Cleaning successful!"
create_drug_concentration_directories "${RESULT_FOLDER}" "${DRUG_NAME}" "${CELL_MODEL}" "${drug_concentrations[@]}"
echo "DrugSimulationPostprocessing Log Start..." >& "${RESULT_FOLDER}/logfile"

# Clear any old PID file
PIDFILE="mpiexec.pid"
rm -f "${PIDFILE}"

# choose the binary based on the value of cell_model
if [[ "${CELL_MODEL}" == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE=drugsim_CiPAORdv1.0_postprocessing
elif [[ "${CELL_MODEL}" == *"ORd-static"* ]]; then
  BINARY_FILE=drugsim_ORd-static_postprocessing
elif [[ "${CELL_MODEL}" == *"ToR-ORd"* ]]; then
  BINARY_FILE=drugsim_ToR-ORd_postprocessing
elif [[ "${CELL_MODEL}" == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE=drugsim_ToR-ORd-dynCl_postprocessing
else
  echo "The cell model ${CELL_MODEL} is not specified to any simulations!!"
  exit 1
fi

echo "Unzipping files..." >> "${RESULT_FOLDER}/logfile" 2>&1
unzip_files "${INITIAL_VALUES_ZIP_FILE}" "./" "${RESULT_FOLDER}" "${SAMPLE_SIZE}" "${drug_concentrations[@]}" >> "${RESULT_FOLDER}/logfile" 2>&1
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Simulation program got some problems!!! Exiting..." >> "${RESULT_FOLDER}/logfile" 2>&1
  exit 1
fi
echo "Unzipping successful!!!" >> "${RESULT_FOLDER}/logfile" 2>&1
mpiexec "-machinefile" "${PBS_NODEFILE}" "-np" "${NPROCS}" "~/marcell/MetaHeart/DrugSimulationTest/bin/${BINARY_FILE}" -input_deck "param.txt" >> "${RESULT_FOLDER}/logfile" 2>&1
echo "Zipping folder..." >> "${RESULT_FOLDER}/logfile" 2>&1
zip_files "${RESULT_FOLDER}" "${TIME_SERIES_SUBSTRING}" "${TIME_SERIES_ZIPNAME}"
zip_files "${RESULT_FOLDER}" "${FEATURES_SUBSTRING}" "${FEATURES_ZIPNAME}"
mv "${TIME_SERIES_ZIPNAME}" "${RESULT_FOLDER}/."
mv "${FEATURES_ZIPNAME}" "${RESULT_FOLDER}/."
echo "Zipping finished" >> "${RESULT_FOLDER}/logfile" 2>&1
echo "Simulation has finished! Check the logfile for more details." >> "${RESULT_FOLDER}/logfile" 2>&1
