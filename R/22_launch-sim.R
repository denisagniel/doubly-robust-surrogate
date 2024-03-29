#!/bin/bash
#SBATCH -n 1                    # Number of cores requested
#SBATCH -t 120:00:00                    # Runtime in minutes
# Or use HH:MM:SS or D-HH:MM:SS, instead of just number of minutes
#SBATCH -p medium                # Partition (queue) to submit to
#SBATCH --open-mode=truncate      # append adds to outfile, truncate deletes first
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o Rout/22_launch-sim.out               # Standard out goes to this file
#SBATCH --mail-type=END         # Mail when the job ends

module load gcc/6.2.0
R CMD BATCH --no-save 21_higher-dim-sim.R Rout/21_higher-dim-sim.Rout
