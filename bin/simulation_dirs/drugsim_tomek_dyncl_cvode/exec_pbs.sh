#PBS -N drugsim_tomek_dyncl_cvode
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20000:00:00
#PBS -e stderr.log
#PBS -o stdout.log
#Specific the shell types
#PBS -S /bin/bash
#Specific the queue type
#PBS -q dque

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

# Use this to export the library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/prog/sundial/lib64:/usr/local/lib64:/usr/lib64
echo $LD_LIBRARY_PATH
export PATH=$PATH
echo $PATH


rm -rf *.log result logfile
mkdir result
mpirun -machinefile $PBS_NODEFILE -np $NPROCS ~/marcell/MetaHeart/DrugSimulation/bin/drugsim_tomek_dyncl -input_deck param.txt > logfile
