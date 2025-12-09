#!/bin/sh

RESULT_FOLDER="./results"
DRUG_NAME=$(grep "^drug_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)
USER_NAME=$(grep "^user_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)
DRUG_CONCENTRATIONS=$(grep "^drug_concentrations" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)
CELL_MODEL=$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

HILL_FILE=$(grep "^hill_file" param.txt | cut -d '=' -f2 | sed 's|//.*||' | xargs)
echo "FILE HILL: ${HILL_FILE}"
SAMPLE_SIZE=$(( $(wc -l < "${HILL_FILE}") - 1 ))
echo "FILE LINE: ${SAMPLE_SIZE}"

RESULT_DRUG_PATH="${RESULT_FOLDER}"
LATEX_FILE="report_drug_${DRUG_NAME}_${CELL_MODEL}.tex"

#Plot all the time-series result from the in-silico simulation
echo "Generate plots from the time-series result"
if [[ "${CELL_MODEL}" == *"Grandi"* ]]; then
  python3 ../scripts/plot_time_series_grandi.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${SAMPLE_SIZE}"
else
  python3 ../scripts/plot_time_series.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${SAMPLE_SIZE}"
fi
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Plot time series got some errors! Please check the logfile_report for more details." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for plot generating.
echo "DrugSimulation Report Progress: ${PROGRESS}%"

#Concat the separated feature data
echo "Unifying feature data"
if [[ "${CELL_MODEL}" == *"Grandi"* ]]; then
  python3 ../scripts/plot_features_grandi.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${USER_NAME}"
else
  python3 ../scripts/plot_features.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${USER_NAME}"
fi
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Plot features got some errors! Please check the logfile_report for more details." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for feature calculation.
echo "DrugSimulation Report Progress: ${PROGRESS}%"

#Calculating the average feature for each concentration
echo "Calculating average feature for each concentration."
python3 ../scripts/compute_dose_averages.py
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Compute dose average got some errors! Please check the logfile_report for more details." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for average feature.
echo "DrugSimulation Report Progress: ${PROGRESS}%"

#Generate report based on the pre-generated LaTEX file
echo "Generate PDF from LaTEX"
cd "${RESULT_FOLDER}"
pdflatex -interaction=nonstopmode -halt-on-error  "${LATEX_FILE}"
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "LaTEX to PDF got some errors! Please check the logfile_report for more details." >> "logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for PDF generating.
echo "DrugSimulation Report Progress: ${PROGRESS}%"
