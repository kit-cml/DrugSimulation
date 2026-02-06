#!/usr/bin/env bash

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
LATEX_FILE="report_template.tex"
MEDIAN_CI_FILE="${RESULT_FOLDER}/${DRUG_NAME}-median-ci.csv"

#Plot all the time-series result from the in-silico simulation
echo "Generate plots from the time-series result"
python3 ../scripts/plot_time_series.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${SAMPLE_SIZE}" "${CELL_MODEL}"
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
python3 ../scripts/plot_features.py "${RESULT_DRUG_PATH}" "${DRUG_NAME}" "${USER_NAME}" "${CELL_MODEL}"
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
python3 ../scripts/combine_and_average_features.py --base_dir "${RESULT_DRUG_PATH}"
echo "Calculating median and Confidence Interval 95% based on averaged sample."
python3 ../scripts/compute_mean_ci.py ./results/sample-averages-across-concentrations.csv "${MEDIAN_CI_FILE}" --n-boot "${SAMPLE_SIZE}"
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "Compute dose average got some errors! Please check the logfile_report for more details." >> "${RESULT_FOLDER}/logfile" 2>&1
  rm -rf "${PIDFILE}"
  exit 1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for average feature.
echo "DrugSimulation Report Progress: ${PROGRESS}%"

#Generate report using Perplexity API with another script.
echo "Generate PDF from LaTEX"
./exec_perplexity_report_generator.sh
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "LaTEX to PDF got some errors! Please check the logfile_report for more details." >> "logfile" 2>&1
fi
PROGRESS=$((PROGRESS + 25)) # 25% for PDF generating.
echo "DrugSimulation Report Progress: ${PROGRESS}%"
