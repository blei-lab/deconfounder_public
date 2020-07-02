#!/bin/bash

#SBATCH --account=sml
#SBATCH -c 4
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=32gb

source /proj/sml_netapp/opt/anaconda2-4.2.0/etc/profile.d/conda.sh
conda activate py2 

echo "python fitfactor.py -nc ${NCAUSES} -nu ${NUNITS} -sim ${SIMSET} -alpha ${ALPHA} -dataseed ${DATASEED} -factorseed ${FACTORSEED} -nitr ${NITR} -fm ${FACTORMODEL}"

python fitfactor.py -nc ${NCAUSES} -nu ${NUNITS} -sim ${SIMSET} -alpha ${ALPHA} -dataseed ${DATASEED} -factorseed ${FACTORSEED} -nitr ${NITR} -fm ${FACTORMODEL}

 