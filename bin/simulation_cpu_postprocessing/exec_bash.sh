#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Error: Provide the number of processor size to the script!"
    echo "Example: ./exec_bash.sh 10"
    exit 1
fi

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH

# Check whether the input is actually a number or not
if [ "$#" -eq 1 ] && echo "$1" | grep -Eq '^-?[0-9]+$';  then
  NUMBER_OF_CPU=$1
else
  echo "The processor size is not a number!!!"
  exit 1
fi

MAX_CPU=$(nproc --all)
if [ $NUMBER_OF_CPU -gt  $MAX_CPU ]; then
  echo "The processor requested exceed the number of maximum of CPU in this PC!"
  exit 1
fi

# to grab cell_model value from parameter file (thanks, ChatGPT).
# grep "^user_name": looks for the line starting with user_name
# cut -d'=' -f2: gets the right-hand side of =
# sed 's/\/\/.*//': removes any inline comment starting with //
# xargs: trims leading and trailing whitespace
cell_model=$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

# choose the binary based on the value of cell_model
if [[ $cell_model == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE=../drugsim_CiPAORdv1.0_postprocessing
elif [[ $cell_model == *"ORdstatic-Dutta"* ]]; then
  BINARY_FILE=../drugsim_ORdstatic-Dutta_postprocessing
elif [[ $cell_model == *"ToR-ORd"* ]]; then
  BINARY_FILE=../drugsim_ToR-ORd_postprocessing
elif [[ $cell_model == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE=../drugsim_ToR-ORd-dynCl_postprocessing
else
  echo "The cell model $cell_model is not specified to any simulations!!"
  exit 1
fi

rm -rf *.log results logfile
mkdir results
echo "Run $cell_model cell model simulation with $NUMBER_OF_CPU cores."
mpiexec -np $NUMBER_OF_CPU $BINARY_FILE -input_deck param.txt
echo "Simulation has finished! Check the logfile for more details."
