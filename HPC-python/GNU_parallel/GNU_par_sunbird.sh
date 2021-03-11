#!/bin/bash --login
###
#SBATCH --ntasks 10 # Number of processors we will use
#SBATCH --output fourier.out.%J # Output file location
#SBATCH --time 00:10:00 # Time limit for this job
#SBATCH --account=scw1389 # specify our current project, change this for your own work
#SBATCH --reservation=scw1389_46 # specify the reservation we have for the training workshop, remove this for your own work, replace XX with the code provided by your instructor
###

module load parallel # Ensure that parallel is available to us
module load anaconda/2020.07 # Load Python and activate the environment we will use for this work
source activate scw_test

export OMP_NUM_THREADS=1 # Only use one thread per copy of Python, since we are using GNU Parallel for parallelism

srun="srun --nodes 1 --ntasks 1" # Define srun arguments:
# --nodes 1 --ntasks 1         allocates a single core to each task

parallel="parallel --max-procs $SLURM_NTASKS --joblog parallel_joblog" # Define parallel arguments:
# --max-procs $SLURM_NTASKS  is the number of concurrent tasks parallel runs, so number of CPUs allocated
# --joblog name     parallel's log file of tasks it has run

# Run the tasks:
# $parallel "$srun python fourier_new.py {1} \
#     --fourier_restricted_output=fourier_restricted_\$(basename {1}).pdf \
#     --noise_isolation_output=noise_isolation_\$(basename {1}).pdf \
#     --phase_contrast_output=phase_contrast_\$(basename {1}).pdf" :::: files_to_process.txt

python write-Omeg_range.py "5.01,7.01,0.10"
$parallel "$srun python write-Omeg_range.py" :::: Omeg_range.txt
python plot.py

# python3 write-Omeg_range.py "5.01,7.02,0.10"
# parallel python3 script_GNU_parallel.py :::: Omeg_range.txt
# python3 plot.py

