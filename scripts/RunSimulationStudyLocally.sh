#!/bin/bash

echo "$(date) Run begins."

SIMUL_PATH=$1
CONFIG_FILE=$2

Rscript --vanilla GenerateInputForSimulatedStudy.R \
	${SIMUL_PATH} 
wait

/Users/seonghwanjun/PhylExAnalysis/PhylEx/run -c ${CONFIG_FILE}

echo "$(date) Run finished."
