#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J scRNA

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SIMUL_PATH=$1
CONFIG_FILE=$2

Rscript --vanilla GenerateInputForSimulatedStudy.R \
	${SIMUL_PATH} 
wait

/proj/sc_ml/users/x_seoju/PhylExAnalysis/PhylEx/run -c ${CONFIG_FILE}

echo "$(date) Run finished."
