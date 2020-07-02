#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)

FACTORDSTD_SWEEP="2"
FACTORZSTD_SWEEP="0.1"
DATASEED_SWEEP="20133224"


# parameters to fit different factor models

SIMSET_SWEEP="indep"
FACTORK_SWEEP="1"
FACTOR_SWEEP="linear"
FACTORSEED_SWEEP="20144219"

# SIMSET_SWEEP="indep"
# FACTORK_SWEEP="1"
# FACTOR_SWEEP="quadratic"
# FACTORSEED_SWEEP="20144133"

# SIMSET_SWEEP="indep"
# FACTORK_SWEEP="2"
# FACTOR_SWEEP="quadratic"
# FACTORSEED_SWEEP="20144219"

# SIMSET_SWEEP="indep"
# FACTORK_SWEEP="3"
# FACTOR_SWEEP="quadratic"
# FACTORSEED_SWEEP="20135747"

# SIMSET_SWEEP="dep"
# FACTORK_SWEEP="1"
# FACTOR_SWEEP="linear"
# FACTORSEED_SWEEP="20144133"


# SIMSET_SWEEP="dep"
# FACTORK_SWEEP="1"
# FACTOR_SWEEP="quadratic"
# FACTORSEED_SWEEP="20144133"


OUT_SUFFIX=".out"
ROUT_SUFFIX=".Rout"

for SIMSETi in ${SIMSET_SWEEP}; do
	export SIMSET=${SIMSETi}
	for FACTORKi in ${FACTORK_SWEEP}; do
		export FACTORK=${FACTORKi}
		for FACTORi in ${FACTOR_SWEEP}; do
			export FACTOR=${FACTORi}
			for FACTORDSTDi in ${FACTORDSTD_SWEEP}; do
				export FACTORDSTD=${FACTORDSTDi}
				for FACTORZSTDi in ${FACTORZSTD_SWEEP}; do
					export FACTORZSTD=${FACTORZSTDi}
					for DATASEEDi in ${DATASEED_SWEEP}; do
						export DATASEED=${DATASEEDi}
						for FACTORSEEDi in ${FACTORSEED_SWEEP}; do
							export FACTORSEED=${FACTORSEEDi}
							# export FACTORSEED=$(date +'%S%6N')
							echo ${FACTORSEED}
					    	export NAME=fitfactor_${SIMSETi}_${FACTORKi}_${FACTORi}_${FACTORDSTDi}_${FACTORZSTDi}_${DATASEEDi}_${FACTORSEEDi}
					    	export OUTNAME=${NAME}_${TIMESTAMP}${OUT_SUFFIX}
					    	export ROUT=${NAME}_${TIMESTAMP}${ROUT_SUFFIX}
				            echo ${NAME}
				            sbatch --job-name=${NAME} \
				            --output=${OUTNAME} \
				            run_script_fitfactor_base.sh	
				        done
				    done
				done				
		    done
		done
	done
done