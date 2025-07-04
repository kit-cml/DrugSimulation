#!/bin/bash

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH

NUMBER_OF_CPU=$(grep "^number_of_cpu" param.txt | cut -d'=' -f2 | sed 's/\/\/.*//' | xargs)
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
CELL_MODEL=$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

# choose the binary based on the value of cell_model
if [[ $CELL_MODEL == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE=../drugsim_CiPAORdv1.0
elif [[ $CELL_MODEL == *"ORd-static"* ]]; then
  BINARY_FILE=../drugsim_ORd-static
elif [[ $CELL_MODEL == *"ToR-ORd"* ]]; then
  BINARY_FILE=../drugsim_ToR-ORd
elif [[ $CELL_MODEL == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE=../drugsim_ToR-ORd-dynCl
else
  echo "The cell model $CELL_MODEL is not specified to any simulations!!"
  exit 1
fi

rm -rf *.log results logfile
mkdir results
echo "Run $CELL_MODEL cell model simulation with $NUMBER_OF_CPU cores."
mpiexec -np $NUMBER_OF_CPU $BINARY_FILE -input_deck param.txt
echo "Simulation has finished! Check the logfile for more details."
