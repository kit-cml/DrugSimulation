#PBS -N drugsim_pbs_job
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

# to grab cell_model value from parameter file (thanks, ChatGPT).
# grep "^user_name": looks for the line starting with user_name
# cut -d'=' -f2: gets the right-hand side of =
# sed 's/\/\/.*//': removes any inline comment starting with //
# xargs: trims leading and trailing whitespace
cell_model=$(grep "^cell_model" param.txt | cut -d'=' -f2 | cut -d'/' -f1 | cut -d'/' -f1 | sed 's/\/\/.*//' | xargs)

# choose the binary based on the value of cell_model
if [[ $cell_model == *"CiPAORdv1.0"* ]]; then
  BINARY_FILE=drugsim_CiPAORdv1.0
elif [[ $cell_model == *"ORd-static"* ]]; then
  BINARY_FILE=drugsim_ORd-static
elif [[ $cell_model == *"ToR-ORd"* ]]; then
  BINARY_FILE=drugsim_ToR-ORd
elif [[ $cell_model == *"ToR-ORd-dynCl"* ]]; then
  BINARY_FILE=drugsim_ToR-ORd-dynCl
else
  echo "The cell model $1 is not specified to any simulations!!"
  exit 1
fi

./clear_workspace.sh
mkdir results
mpiexec -machinefile $PBS_NODEFILE -np $NPROCS ~/marcell/MetaHeart/DrugSimulationTest/bin/$BINARY_FILE -input_deck param.txt > logfile
