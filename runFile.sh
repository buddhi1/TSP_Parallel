#!/bin/bash
#SBATCH -J test           # job name
#SBATCH -o output.out       # output file name (%j expands to jobID), this file captures standered output from the shell
#SBATCH -e output.err      # error file name (%j expands to jobID), this file captures standered errors genereted from the program
#SBATCH  --nodes  8      # total number of nodes requested
#SBATCH -p qCDER            # partition --qCDER (to find out available partitions please run 'sinfo' command)
#SBATCH -t 1:30:00       # run time (hh:mm:ss) - 0.5 hours
 
## Included any modules that may be required for your job. in this example load MPI module
module load Compilers/mvapich2_ACoRE
 
 
# execute any program
#to run MPI program
for i in 2 4 8; do
	for (( j = 2; j <= 5; j++ )); do
		srun --mpi=pmi2 -n $i ./mpi2.out	$j
	done
done
# srun --mpi=pmi2 -N 2 -n 4 ./mpi.out
# srun --mpi=pmi2 -N 2 -n 8 ./mpi.out
# srun --mpi=pmi2 -N 2 -n 16 ./mpi.out
# srun --mpi=pmi2 -N 4 -n 2 ./mpi.out
# srun --mpi=pmi2 -N 4 -n 4 ./mpi.out
# srun --mpi=pmi2 -N 4 -n 8 ./mpi.out
# srun --mpi=pmi2 -N 4 -n 16 ./mpi.out
# srun --mpi=pmi2 -N 8 -n 2 ./mpi.out
# srun --mpi=pmi2 -N 8 -n 4 ./mpi.out
# srun --mpi=pmi2 -N 8 -n 8 ./mpi.out
# srun --mpi=pmi2 -N 8 -n 16 ./mpi.out

# ./MPI_scatter is run using 8 nodes and 4 tasks per node
# srun --mpi=pmi2 -N 8 -n 4 ./MPI_scatter
