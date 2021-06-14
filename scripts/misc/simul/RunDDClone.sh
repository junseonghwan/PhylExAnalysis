#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J ddClone

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SEED=$1
SIMUL_PATH=$2
MCMC_ITER=$3
SC_VAR_THRESHOLD=$4

Rscript --vanilla Rscripts/RunDDClone.R \
	${SEED} ${SIMUL_PATH} ${MCMC_ITER} ${SC_VAR_THRESHOLD}

echo "$(date) Run finished."
