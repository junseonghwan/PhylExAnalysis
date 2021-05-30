#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J Canopy

echo "$(date) Run begins."

#module load R/3.6.3-nsc1-gcc-7.3.0
#wait

SEED=$1
DATA_PATH=$2
K_BEGIN=$3
K_END=$4

Rscript --vanilla Rscripts/RunCanopy.R \
	${SEED} ${DATA_PATH} ${K_BEGIN} ${K_END}

echo "$(date) Run finished."
