#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem 1G
#SBATCH -J ProcessSingleRegion

echo "$(date) Run begins."

module load R/3.6.3-nsc1-gcc-7.3.0
wait

SIMUL_PATH=$1

Rscript --vanilla PrepareSingleRegionData.R \
	${SIMUL_PATH}

echo "$(date) Run finished."
