#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)

DATASEED_SWEEP="1549139415"
NTRIAL_SWEEP="100"
CONFSCALE_SWEEP="1"
OUTCOME_SWEEP="linear"
SIMSET_SWEEP="indep"

# DATASEED_SWEEP="20133224"
# NTRIAL_SWEEP="100"
# CONFSCALE_SWEEP="1"
# OUTCOME_SWEEP="linear"
# SIMSET_SWEEP="dep"

OUT_SUFFIX=".out"
ROUT_SUFFIX=".Rout"

for DATASEEDi in ${DATASEED_SWEEP}; do
	export DATASEED=${DATASEEDi}
	for SIMSETi in ${SIMSET_SWEEP}; do
		export SIMSET=${SIMSETi}
		for CONFSCALEi in ${CONFSCALE_SWEEP}; do
			export CONFSCALE=${CONFSCALEi}
			for NTRIALi in ${NTRIAL_SWEEP}; do
				export NTRIAL=${NTRIALi}
				for OUTCOMEi in ${OUTCOME_SWEEP}; do
					export OUTCOME=${OUTCOMEi}
			    	export NAME=simdat_${DATASEEDi}_${SIMSETi}_${CONFSCALEi}_${NTRIALi}_${OUTCOMEi}
			    	export OUTNAME=${NAME}_${TIMESTAMP}${OUT_SUFFIX}
			    	export ROUT=${NAME}_${TIMESTAMP}${ROUT_SUFFIX}
		            echo ${NAME}
		            sbatch --job-name=${NAME} \
		            --output=${OUTNAME} \
		            run_script_simdat_base.sh					
			    done
			done
		done
	done
done