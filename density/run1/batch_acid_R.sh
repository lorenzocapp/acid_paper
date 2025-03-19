#!/bin/bash

#SBATCH --job-name=k			    ### Job Name
#SBATCH -o /homes/users/lcappello/scratch/output/k_%A_%a.out	  ### File in which to store job output/error
#SBATCH -e  /homes/users/lcappello/scratch/error/k_%A_%a.err        ### error file messages. 
#SBATCH --time=00-48:00:00         			 ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --requeue                  		 ### On failure, requeue for another try
#SBATCH --verbose               	    ### Increase informational messages
#SBATCH --array=1-100   	        ### Array index | %50: number of simultaneously tasks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lorenzo.cappello@upf.edu 
#SBATCH --partition=haswell

echo
echo "****************************************************************************"
echo "*                                                                          *"
echo "********************** sbatch script for array job *************************"
echo "*                                                                          *"
echo "****************************************************************************"
echo



current_dir=${PWD##*/}
echo "Current dir: $current_dir"
echo
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
pwd
echo

# First we ensure a clean running environment:
module purge
module load modulepath/noarch
# Load R
#module load R/3.6.0-foss-2018b
ml R/4.2.0-foss-2021b

## Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# Output file
outFile="output_parameter_${i}.txt"

# Pass line #i to a R script 
Rscript --vanilla master_kernel.R ${i} ${outFile}

echo
echo '******************** FINISHED ***********************'
echo
