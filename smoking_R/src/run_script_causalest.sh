#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)

DATASEED_SWEEP="20133224"

SIMSET_SWEEP="indep"
FACTOR_SWEEP="linear"
FACTORK_SWEEP="1"
ESTCTRL_SWEEP="Actrl Zctrl nctrl oracle" 
FACTORSEED_SWEEP="20144219" 

# SIMSET_SWEEP="indep"
# FACTOR_SWEEP="quadratic"
# FACTORK_SWEEP="1"
# ESTCTRL_SWEEP="Actrl Zctrl Acovctrl Zcovctrl"
# FACTORSEED_SWEEP="20144133" 

# SIMSET_SWEEP="indep"
# FACTOR_SWEEP="quadratic"
# FACTORK_SWEEP="2"
# ESTCTRL_SWEEP="Actrl Zctrl"
# FACTORSEED_SWEEP="20144219" 

# SIMSET_SWEEP="indep"
# FACTOR_SWEEP="quadratic"
# FACTORK_SWEEP="3"
# ESTCTRL_SWEEP="Actrl Zctrl"
# FACTORSEED_SWEEP="20135747" 

# SIMSET_SWEEP="dep"
# FACTOR_SWEEP="linear quadratic"
# FACTORK_SWEEP="1"
# ESTCTRL_SWEEP="Actrl Zctrl nctrl oracle" 
# FACTORSEED_SWEEP="20144133" 

TRIAL_SWEEP=$(seq 1 100)

RANDSEED_SWEEP=$(date +%d%H%M%S)
ALGO_SWEEP="meanfield"
PRIORSP_SWEEP="horseshoe" # alternatively use normal
OUTCOME_SWEEP="linear" 

OUT_SUFFIX=".out"
ROUT_SUFFIX=".Rout"

for DATASEEDi in ${DATASEED_SWEEP}; do
	export DATASEED=${DATASEEDi}
	for SIMSETi in ${SIMSET_SWEEP}; do
		export SIMSET=${SIMSETi}
		for FACTORi in ${FACTOR_SWEEP}; do
			export FACTOR=${FACTORi}
			for FACTORKi in ${FACTORK_SWEEP}; do
				export FACTORK=${FACTORKi}
				for ESTCTRLi in ${ESTCTRL_SWEEP}; do
					export ESTCTRL=${ESTCTRLi}
					for TRIALi in ${TRIAL_SWEEP}; do
						export TRIAL=${TRIALi}
						for FACTORSEEDi in ${FACTORSEED_SWEEP}; do
							export FACTORSEED=${FACTORSEEDi}
							for RANDSEEDi in ${RANDSEED_SWEEP}; do
								export RANDSEED=${RANDSEEDi}
								for ALGOi in ${ALGO_SWEEP}; do
									export ALGO=${ALGOi}
									for PRIORSPi in ${PRIORSP_SWEEP}; do
										export PRIORSP=${PRIORSPi}
										for OUTCOMEi in ${OUTCOME_SWEEP}; do
											export OUTCOME=${OUTCOMEi}
									    	export NAME=causalest_${DATASEEDi}_${SIMSETi}_${FACTORi}_${FACTORKi}_${ESTCTRLi}_${TRIALi}_${FACTORSEEDi}_${RANDSEEDi}_${ALGOi}_${PRIORSPi}_${OUTCOMEi}
									    	export OUTNAME=${NAME}_${TIMESTAMP}${OUT_SUFFIX}
									    	export ROUT=${NAME}_${TIMESTAMP}${ROUT_SUFFIX}
								            echo ${NAME}
								            sbatch --job-name=${NAME} \
								            --output=${OUTNAME} \
								            run_script_causalest_base.sh
								        done
								    done
								done
							done
						done
					done
			    done
			done
		done
	done
done