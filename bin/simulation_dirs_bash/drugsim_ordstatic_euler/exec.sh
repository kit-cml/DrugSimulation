#!/bin/bash


# Use this to export the library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH

# Check the maximum number of CPU and use them all.
# If you want to set it up manually, comment the line below
# until the 'fi' and uncomment the line after the 'fi'.
NUMBER_OF_CPU=$(nproc 2>/dev/null)
if [[ $? -ne 0 ]]; then
    echo "Error: Command 'nproc' not found or failed to execute. Make sure to install nproc first."
    exit 1
fi
#NUMBER_OF_CPU=1

#find . -name "*.plt" -type f -delete
rm -rf *.log result log
mkdir result
echo "Run simulation with $NUMBER_OF_CPU cores."
mpirun -np $NUMBER_OF_CPU ~/marcell/MetaHeart/DrugSimulation/bin/drugsim_ordstatic -input_deck param.txt > logfile
