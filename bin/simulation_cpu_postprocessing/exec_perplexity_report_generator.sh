#!/usr/bin/env bash
set -euo pipefail

# Load external functions
source ../scripts/create_escaped_b64_file.sh

# Variable Assignment Section
#--------------------------------------------
RESULT_FOLDER="results"
RETURN_CODE=1
NUMBER_OF_LOOPS=1

PARAM_FILE="param.txt"
PARAM_CAMEL_FILE="param-camel.txt"
PARAM_CAMEL_ESCAPED_FILE="param-camel-escaped.txt"
PARAM_CAMEL_ESCAPED_B64="param-camel-escaped.b64"

CELL_MODEL=$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)
DRUG_NAME=$(grep "^drug_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)
DRUG_CONCENTRATIONS=$(grep "^drug_concentrations" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

MEDIAN_CI_FILE="${RESULT_FOLDER}/${DRUG_NAME}-median-ci.csv"
MEDIAN_CI_ESCAPED_FILE="${RESULT_FOLDER}/${DRUG_NAME}-median-ci-escaped.csv"
MEDIAN_CI_ESCAPED_B64="${RESULT_FOLDER}/${DRUG_NAME}-median-ci-escaped.b64"

TIME_SERIES_CONTROL_FILE="${RESULT_FOLDER}/0.00/${DRUG_NAME}_0.00_time_series_smp0.csv"
TIME_SERIES_CONTROL_ESCAPED_FILE="${RESULT_FOLDER}/0.00/${DRUG_NAME}_0.00_time_series_smp0-escaped.csv"
TIME_SERIES_CONTROL_ESCAPED_B64="${RESULT_FOLDER}/0.00/${DRUG_NAME}_0.00_time_series_smp0-escaped.b64"

CONC_MAX=$(tr ',' '\n' <<< "$DRUG_CONCENTRATIONS" | sort -n | tail -1)
printf -v CONC_FORMAT '%.2lf' "${CONC_MAX}"
TIME_SERIES_DRUG_FILE="${RESULT_FOLDER}/${CONC_FORMAT}/${DRUG_NAME}_${CONC_FORMAT}_time_series_smp0.csv"
TIME_SERIES_DRUG_ESCAPED_FILE="${RESULT_FOLDER}/${CONC_FORMAT}/${DRUG_NAME}_${CONC_FORMAT}_time_series_smp0-escaped.csv"
TIME_SERIES_DRUG_ESCAPED_B64="${RESULT_FOLDER}/${CONC_FORMAT}/${DRUG_NAME}_${CONC_FORMAT}_time_series_smp0-escaped.b64"

LATEX_TEMPLATE_FILE="${RESULT_FOLDER}/report_template.tex"
LATEX_TEMPLATE_ESCAPED_FILE="${RESULT_FOLDER}/report_template-escaped.tex"
LATEX_TEMPLATE_ESCAPED_B64="${RESULT_FOLDER}/report_template-escaped.b64"

LATEX_FILE="report_drug_${DRUG_NAME}_${CELL_MODEL}.tex"

#--------------------------------------------

echo "Preparing escaped files with their respective Base64 files..."

# Create parameter files with camel case, 
# and convert it as a Base64 file.
../scripts/convert_params.sh "${PARAM_FILE}" "${PARAM_CAMEL_FILE}"

# Call this function to generate both escaped files and base64 of escaped files.
create_escaped_b64_file \
  "${PARAM_CAMEL_FILE}" \
  "${PARAM_CAMEL_ESCAPED_FILE}" \
  "${PARAM_CAMEL_ESCAPED_B64}"

create_escaped_b64_file \
  "${MEDIAN_CI_FILE}" \
  "${MEDIAN_CI_ESCAPED_FILE}" \
  "${MEDIAN_CI_ESCAPED_B64}"

create_escaped_b64_file \
  "${TIME_SERIES_CONTROL_FILE}" \
  "${TIME_SERIES_CONTROL_ESCAPED_FILE}" \
  "${TIME_SERIES_CONTROL_ESCAPED_B64}"

create_escaped_b64_file \
  "${TIME_SERIES_DRUG_FILE}" \
  "${TIME_SERIES_DRUG_ESCAPED_FILE}" \
  "${TIME_SERIES_DRUG_ESCAPED_B64}"

create_escaped_b64_file \
  "${LATEX_TEMPLATE_FILE}" \
  "${LATEX_TEMPLATE_ESCAPED_FILE}" \
  "${LATEX_TEMPLATE_ESCAPED_B64}"

# main loop for creating the AI report
while [[ "${RETURN_CODE}" -ne 0  ]] && [[ "${NUMBER_OF_LOOPS}" -lt 4 ]]
do
  echo "Call the report generator..."
  # Execute the report generator
  ../perplexity_report_generator -input_deck "${PARAM_FILE}"
  RETURN_CODE=$?
  if [ $RETURN_CODE -ne 0 ]; then
    echo "AI Report Generator Program got problem!!! Exiting..."
    exit 1
  fi

  # go to the RESULT_FOLDER
  cd "${RESULT_FOLDER}"

  # call the latex processor command
  # need to give set +e and set -e because there is error guard.
  set +e
  /opt/texlive/2025/bin/x86_64-linux/pdflatex -interaction=nonstopmode -halt-on-error  "${LATEX_FILE}"
  RETURN_CODE=$?
  set -e
  echo "LaTeX return code: ${RETURN_CODE}"
  # run another time to make sure the hypertext is properly attached
  if [[ "${RETURN_CODE}" -eq 0 ]]; then
    /opt/texlive/2025/bin/x86_64-linux/pdflatex -interaction=nonstopmode -halt-on-error  "${LATEX_FILE}"
  fi

  if [[ "${RETURN_CODE}" -ne 0 ]]; then
    echo "Return was ${RETURN_CODE} at attempt ${NUMBER_OF_LOOPS}, process is repeated"
    NUMBER_OF_LOOPS=$((NUMBER_OF_LOOPS + 1))
  fi

  cd ..  
done

if [[ "${RETURN_CODE}" -eq 0 ]]; then
  echo "Return was 0 after ${NUMBER_OF_LOOPS} times. Process done!!"
  exit 0
else
  echo "Report failed to generate even after $((${NUMBER_OF_LOOPS} - 1)) times. Please press the Generate Report button to re-generate the report!!"
  exit 1
fi
