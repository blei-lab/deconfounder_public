#!/bin/bash

#SBATCH --account=sml
#SBATCH -c 2
#SBATCH --time=11:00:00
#SBATCH --mem-per-cpu=16gb

source /proj/sml_netapp/opt/anaconda2-4.2.0/etc/profile.d/conda.sh
conda activate py2 

# simulate data

echo "python simdat.py -nc ${NCAUSES} -nu ${NUNITS} -sim ${SIMSET} -alpha ${ALPHA} -seed ${DATASEED}"

python simdat.py -nc ${NCAUSES} -nu ${NUNITS} -sim ${SIMSET} -alpha ${ALPHA} -seed ${DATASEED}
