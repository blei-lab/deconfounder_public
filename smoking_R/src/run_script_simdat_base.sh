#!/bin/sh
#
#
#SBATCH --account=stats      # The account name for the job.
#SBATCH -c 4                     # The number of cpu cores to use.
#SBATCH --time=2:00:00             # The time the job will take to run.
#SBATCH --mem-per-cpu=4gb        # The memory the job will use per cpu core.

module load R


# simulate data

Rscript --verbose simdat.R ${NTRIAL} ${CONFSCALE} ${OUTCOME} ${SIMSET} ${DATASEED}
