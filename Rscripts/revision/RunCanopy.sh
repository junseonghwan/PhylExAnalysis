#!/bin/bash

#SBATCH -A snic2021-5-259
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem 8G
#SBATCH -J RunCanopy

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SEED=$1
SIMUL_PATH=$2

Rscript --vanilla RunCanopy.R \
	${SEED} ${SIMUL_PATH}

echo "$(date) Run finished."
