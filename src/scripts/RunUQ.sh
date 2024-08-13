#!/bin/bash

declare -a Rs=("5e5" "1e6" "5e6")
declare -a OPPs=("0.8" "1.2" "1.0")

for R in "${Rs[@]}"
do
    for OPP in "${OPPs[@]}"
    do
	time python3.10 TODelete.py ../../Results/DataForPaper/FewerTermsInStage1/PopulationParameters.csv $R $OPP
    done
done
