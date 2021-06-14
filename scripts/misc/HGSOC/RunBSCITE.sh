#!/bin/bash

echo "$(date) Run begins."

DATA_PATH=$1
SC_MUT_MATRIX_PATH=$2

Rscript --vanilla Rscripts/GenerateInputForBScite.R \
	${DATA_PATH} ${SC_MUT_MATRIX_PATH}

echo "$(date) Run finished."
