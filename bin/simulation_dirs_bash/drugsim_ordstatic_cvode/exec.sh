#!/bin/bash


# Use this to export the library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundials/sundials-5.7.0/lib:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH


#find . -name "*.plt" -type f -delete
rm -rf *.log result log
mkdir result
mpirun -np 24 ~/marcell/MetaHeart/DrugSimulation/bin/drugsim_ordstatic -input_deck param.txt > logfile
