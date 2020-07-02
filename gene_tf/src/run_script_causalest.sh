#!/bin/bash
#

TIMESTAMP=$(date +%Y%m%d%H%M%S)

# high SNR

NCAUSES_SWEEP="5000"
NUNITS_SWEEP="5000"
SIMSET_SWEEP="BN TGP HGDP"
ALPHA_SWEEP="10"

# NCAUSES_SWEEP="5000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="PSD"
# ALPHA_SWEEP="1 10 50 100"

# NCAUSES_SWEEP="5000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="SP"
# ALPHA_SWEEP="10 25 50 100"

SNPSIG_SWEEP="40"
CONFINT_SWEEP="40"
CAUSALPROP_SWEEP="10"
ASLIN_SWEEP="4"
ASLOG_SWEEP="0"
AALIN_SWEEP="0"
AALOG_SWEEP="0.0001"


# low SNR

# NCAUSES_SWEEP="100000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="BN TGP HGDP"
# ALPHA_SWEEP="10"

# NCAUSES_SWEEP="100000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="PSD"
# ALPHA_SWEEP="1 10 50 100"

# NCAUSES_SWEEP="100000"
# NUNITS_SWEEP="5000"
# SIMSET_SWEEP="SP"
# ALPHA_SWEEP="10 25 50 100"

# SNPSIG_SWEEP="10"
# CONFINT_SWEEP="20"
# CAUSALPROP_SWEEP="10"
# ASLIN_SWEEP="100"
# ASLOG_SWEEP="10"
# AALIN_SWEEP="0"
# AALOG_SWEEP="0.0001"


# FACTORSEED_SWEEP=$(date +%d%H%M%S)

FACTORSEED_SWEEP="29200422"
DATASEED_SWEEP="29200017"
NITR_SWEEP="20000"
CV_SWEEP="0"


RUN_SCRIPT="run_script_causalest_base_emayhem.sh"
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
                            for CVi in ${CV_SWEEP}; do
                                export CV=${CVi}
                                for SNPSIGi in ${SNPSIG_SWEEP}; do
                                    export SNPSIG=${SNPSIGi}
                                    for CONFINTi in ${CONFINT_SWEEP}; do
                                        export CONFINT=${CONFINTi}
                                        for CAUSALPROPi in ${CAUSALPROP_SWEEP}; do
                                            export CAUSALPROP=${CAUSALPROPi}
                                            for ASLINi in ${ASLIN_SWEEP}; do
                                                export ASLIN=${ASLINi}
                                                for ASLOGi in ${ASLOG_SWEEP}; do
                                                    export ASLOG=${ASLOGi}
                                                    for AALINi in ${AALIN_SWEEP}; do
                                                        export AALIN=${AALINi}
                                                        for AALOGi in ${AALOG_SWEEP}; do
                                                            export AALOG=${AALOGi}
                                                            for i in $(seq 1 100); do
                                                                export NAME=causalest_${DATASEEDi}_${FACTORSEEDi}_${SIMSETi}_${ALPHAi}_${NUNITSi}_${NCAUSESi}_nitr${NITRi}_cv${CVi}_SNPSIG${SNPSIGi}_confint${CONFINTi}_cp${CAUSALPROPi}_aslin${ASLINi}_aslog${ASLOGi}_aalin${AALINi}_aalog${AALOGi}
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
                            done
                        done
                    done
                done
            done
        done
    done
done

