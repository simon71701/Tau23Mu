#!/bin/bash
# Usage:
#    submitAllJobs.sh <N_files> <AnalysisType>

helpstring="Usage:
submitAllJobs.sh [N_files] [AnalysisType] [OutputName]

N_files: Number of ntuples to be used per job
AnalysisType: Choose between 'tau3mu' and 'control'
"

N_FILES=$1
ANALYSISTYPE=$2
YEAR=$3
OUT_NAME=$4


# Check inputs
if [ -z ${3+x} ]
then
echo ${helpstring}
return
fi

# Control channel
if [ ${ANALYSISTYPE} == "control" ]
then
    echo -e "\nDsPhiPi analysis"
    DATATYPE='data_control'
elif [ ${ANALYSISTYPE} == "tau3mu" ]
then
    echo -e "\nTau3mu analysis"
    DATATYPE='data'
elif [ ${ANALYSISTYPE} == "phimunu" ]
then
    echo -e "\nPhiMuNu analysis"
    DATATYPE='data_phimunu'
fi


for i in 2022C_0 2022C_1 2022C_2 2022C_3 2022C_4 2022C_5 2022C_6 2022C_7 2022D_0 2022D_1 2022D_2 2022D_3 2022D_4 2022D_5 2022D_6 2022D_7 2022E_0 2022E_1 2022E_2 2022E_3 2022E_4 2022E_5 2022E_6 2022E_7 2022F_0 2022F_1 2022F_2 2022F_3 2022F_4 2022F_5 2022F_6 2022F_7 2022G_0 2022G_1 2022G_2 2022G_3 2022G_4 2022G_5 2022G_6 2022G_7
do
    echo -e "\nData $i"
    echo -e "\npython3 createRunFile_new.py --run $i --n ${N_FILES} ${DATATYPE} ${YEAR} ${ANALYSISTYPE} --outName ${OUT_NAME}"
    python3 createRunFile_new.py --run $i --n ${N_FILES} ${DATATYPE} ${YEAR} ${ANALYSISTYPE} --outName ${OUT_NAME}
    sleep 1
    echo -e "\nsubmit_analysis_${i}_${ANALYSISTYPE}_${OUT_NAME}.sh"
    echo "$(pwd)"
    source submit_analysis_${i}_${ANALYSISTYPE}_${OUT_NAME}.sh
    cd ..
    sleep 1
done
