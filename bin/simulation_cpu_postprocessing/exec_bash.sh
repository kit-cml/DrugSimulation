#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Error: Provide the cell model name and number of processor size to the script!"
    echo "Example: ./exec_bash.sh ord 10"
    exit 1
fi

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH

# Check whether the input is actually a number or not
if [ "$#" -eq 2 ] && echo "$2" | grep -Eq '^-?[0-9]+$';  then
  NUMBER_OF_CPU=$2
else
  echo "The processor size is not a number!!!"
  exit 1
fi

MAX_CPU=$(nproc --all)
if [ $NUMBER_OF_CPU -gt  $MAX_CPU ]; then
  echo "The processor requested exceed the number of maximum of CPU in this PC!"
  exit 1
fi

if [ $1 == "ord" ]; then
  CAPTION="CiPAORdv1.0 DrugSim."
  BINARY_FILE=../drugsim_ord_postprocessing
elif [ $1 == "ordstatic" ]; then
  CAPTION="ORd2011-Dutta DrugSim."
  BINARY_FILE=../drugsim_ordstatic_postprocessing
elif [ $1 == "tomek" ]; then
  CAPTION="ToR-ORd DrugSim."
  BINARY_FILE=../drugsim_tomek_postprocessing
else
  echo "The cell model $1 is not specified to any simulations!!"
  exit 1
fi

rm -rf *.log results logfile
mkdir results
echo $CAPTION
echo "Run postprocessing simulation with $NUMBER_OF_CPU cores."
mpiexec -np $NUMBER_OF_CPU $BINARY_FILE -input_deck param.txt
echo "Simulation has finished! Check the logfile for more details."
