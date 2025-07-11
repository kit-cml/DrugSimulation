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

RESULT_DRUG_PATH=${RESULT_FOLDER}/${DRUG_NAME}_${CELL_MODEL}/
LATEX_FILE=report_drug_${DRUG_NAME}_${CELL_MODEL}_${USER_NAME}.tex

echo "TESTING RESULT: $RESULT_DRUG_PATH"
#Plot all the time-series result from the in-silico simulation
echo Generate plots from the time-series result
python3 ../scripts/plot_time_series.py $RESULT_DRUG_PATH $SAMPLE_SIZE
PROGRESS=$((PROGRESS + 25)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $PROGRESS%"

#Concat the separated feature data
echo Unifying feature data
python3 ../scripts/plot_features.py $RESULT_DRUG_PATH $USER_NAME
PROGRESS=$((PROGRESS + 25)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $PROGRESS%"

#Generate report based on the pre-generated LaTEX file
echo "Generate PDF from LaTEX (on construction)"
pdflatex $LATEX_FILE
PROGRESS=$((PROGRESS + 50)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $PROGRESS%"
