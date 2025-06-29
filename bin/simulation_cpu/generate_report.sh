#!/bin/sh

if [ "$#" -lt 3 ]; then
    echo "Error: Provide result directory, number of samples, and LaTEX files sequentially"
    echo "Example: ./generate_report.sh results/quinidine/ 10 report_drug_quinidine_ORdstatic-Dutta_epi_marcell.tex"
    exit 1
fi

result_folder_path=$1

# to grab user_name value from parameter file (thanks, ChatGPT).
# grep "^user_name": looks for the line starting with user_name
# cut -d'=' -f2: gets the right-hand side of =
# sed 's/\/\/.*//': removes any inline comment starting with //
# xargs: trims leading and trailing whitespace
user_name=$(grep "^user_name" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

number_of_sample=$2

tex_file_template=$3

#Plot all the time-series result from the in-silico simulation
echo Generate plots from the time-series result
python3 ../scripts/plot_time_series.py $result_folder_path $2
progress=$((progress + 25)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $progress%"

#Concat the separated feature data
echo Unifying feature data
python3 ../scripts/plot_features.py $result_folder_path $user_name
progress=$((progress + 25)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $progress%"

#Generate report based on the pre-generated LaTEX file
echo "Generate PDF from LaTEX (on construction)"
pdflatex $tex_file_template
progress=$((progress + 50)) # 25% for plot generating.
echo "DrugSimulation Report Progress: $progress%"
