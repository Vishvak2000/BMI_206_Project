#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -o /wynton/home/corces/vishvak/SCENT_dataset/job_logs/
#$ -e /wynton/home/corces/vishvak/SCENT_dataset/job_logs/
#$ -l mem_free=1G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=48:00:00   # job requires up to 24 hours of runtime
#$ -t 1-1000           # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted
#$ -pe smp 6

## If you array jobs (option -t), this script will run T times, once per task.
## For each run, $SGE_TASK_ID is set to the corresponding task index (here 1-10).
## To configure different parameters for each task index, one can use a Bash 
## array to map from the task index to a parameter string.
## Select the parameter for the current task index
## Arrays are indexed from 0, so we subtract one from the task index
# param="${params[$((SGE_TASK_ID - 1))]}"
conda activate scverse
num_cores=6
file_SCENT_obj=/wynton/home/corces/vishvak/SCENT_dataset/fibroblast_SCENT_no_covariates.rds
celltype=fibroblast
regr=poisson
bin=FALSE
output_dir=/wynton/home/corces/vishvak/SCENT_dataset/SCENT_no_covariates_output/

Rscript SCENT_parallelization.R $SGE_TASK_ID ${num_cores} ${file_SCENT_obj} ${celltype} ${regr} ${bin} ${output_dir}
