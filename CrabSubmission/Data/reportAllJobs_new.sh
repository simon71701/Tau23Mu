#!/bin/sh
# Usage:
#    reportAllJobs.sh

cd SkimTau3mu
for i in 2022_RunC 2022_RunD-v1 2022_RunD-v2 2022_RunE 2022_RunF 2022_RunG
do
    echo -e "\nData $i"
    cd reco${i}
    source report_v6.sh
    cd ..
    sleep 1
done

cd ..
cd SkimPhiPi
for i in 2022_RunC 2022_RunD-v1 2022_RunD-v2 2022_RunE 2022_RunF 2022_RunG
do
    echo -e "\nData $i"
    cd reco${i}
    source report_v6.sh
    cd ..
    sleep 1
done
