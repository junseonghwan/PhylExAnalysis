#!/bin/bash

echo "$(date) Run begins."

SIMUL_PATH=$1
SC_THRESHOLD=$2

Rscript --vanilla Rscripts/GenerateInputForBScite.R \
	${SIMUL_PATH} ${SC_THRESHOLD}

echo "$(date) Run finished."
