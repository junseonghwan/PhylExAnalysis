#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J ddClone

echo "$(date) Run begins."

#module load R/3.6.3-nsc1-gcc-7.3.0
#wait

SEED=$1
DATA_PATH=$2
MUT_MATRIX_PATH=$3
MCMC_ITER=$4

Rscript --vanilla Rscripts/RunDDClone.R \
	${SEED} ${DATA_PATH} ${MUT_MATRIX_PATH} ${MCMC_ITER}

echo "$(date) Run finished."
