#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Error: Provide the processor size to the script!"
    exit 1
fi

# Use this to export the library path.
# Please change the directory according to your library's location.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH

# Check the maximum number of CPU and use them all.
# If you want to set it up manually, comment the line below
# until the 'fi' and uncomment the line after the 'fi'.
#NUMBER_OF_CPU=$(nproc 2>/dev/null)
#if [[ $? -ne 0 ]]; then
#    echo "Error: Command 'nproc' not found or failed to execute. Make sure to install nproc first."
#    exit 1
#fi
#NUMBER_OF_CPU=10
if [[ "$1" =~ ^-?[0-9]+$ ]]; then
  NUMBER_OF_CPU=$1
else
  echo "The processor size is not a number!!!"
  exit 1
fi


rm -rf *.log results logfile
mkdir results
echo "ToR-ORd DrugSim Euler."
echo "Run simulation with $NUMBER_OF_CPU cores."
mpirun -np $NUMBER_OF_CPU ../../drugsim_tomek -input_deck param.txt
echo "Simulation has finished! Check the logfile for more details."
