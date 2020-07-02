#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)


NCAUSES_SWEEP="5000 100000"
NUNITS_SWEEP="5000"
SIMSET_SWEEP="BN TGP HGDP"
ALPHA_SWEEP="10"

# NCAUSES_SWEEP="5000 100000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="PSD"
# ALPHA_SWEEP="1 10 50 100"

# NCAUSES_SWEEP="5000 100000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="SP"
# ALPHA_SWEEP="10 25 50 100"


FACTORSEED_SWEEP="29200422"
DATASEED_SWEEP="29200017"
NITR_SWEEP="10000"
FACTORMODEL_SWEEP="PCA PPCA PF DEF GMM LFA"

# FACTORSEED_SWEEP=$(date +%d%H%M%S)

RUN_SCRIPT="run_script_fitfactor_base_emayhem.sh"
OUT_SUFFIX=".out"

for DATASEEDi in ${DATASEED_SWEEP}; do
    export DATASEED=${DATASEEDi}
    for SIMSETi in ${SIMSET_SWEEP}; do
        export SIMSET=${SIMSETi}
        for NCAUSESi in ${NCAUSES_SWEEP}; do
            export NCAUSES=${NCAUSESi}
            for NUNITSi in ${NUNITS_SWEEP}; do
                export NUNITS=${NUNITSi}
                for ALPHAi in ${ALPHA_SWEEP}; do
                    export ALPHA=${ALPHAi}
                    for FACTORSEEDi in ${FACTORSEED_SWEEP}; do
                        export FACTORSEED=${FACTORSEEDi}
                        for NITRi in ${NITR_SWEEP}; do
                            export NITR=${NITRi}
                            for FACTORMODELi in ${FACTORMODEL_SWEEP}; do
                                export FACTORMODEL=${FACTORMODELi}
                                export NAME=fitfactor_${DATASEEDi}_${FACTORSEEDi}_${SIMSETi}_${ALPHAi}_${NUNITSi}_${NCAUSESi}_model${FACTORMODELi}_nitr${NITRi}
                                export OUTNAME=${NAME}_${TIMESTAMP}${OUT_SUFFIX}
                                echo ${NAME}
                                sbatch --job-name=${NAME} \
                                --output=${OUTNAME} \
                                ${RUN_SCRIPT}   
                            done
                        done
                    done
                done
            done
        done
    done
done