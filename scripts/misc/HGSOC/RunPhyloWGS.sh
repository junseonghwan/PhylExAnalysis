#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mem 4G
#SBATCH -J PhyloWGS

module load R/3.6.3-nsc1-gcc-7.3.0
module load Python/2.7.14-nsc1-gcc-2018a-eb

DATA_PATH=$1
PWGS_SSM_PATH=$2
PWGS_CNV_PATH=$3
OUTPUT_PATH=$4
BURN_IN=$5
MCMC_ITER=$6
MH_ITER=$7
SEED=$8

# Rscript --vanilla Rscripts/GenerateInputForPhyloWGS.R \
# 	${DATA_PATH} 
# wait

#python /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/phylowgs/evolve.py ${PWGS_SSM_PATH} \
#	${PWGS_CNV_PATH} -O ${OUTPUT_PATH} -B ${BURN_IN} -s ${MCMC_ITER} -i ${MH_ITER} -r ${SEED}

python /proj/sc_ml/users/x_seoju/ScRNACloneEvaluation/phylowgs/multievolve.py --num-chains 4 --ssms ${PWGS_SSM_PATH} \
	--cnvs ${PWGS_CNV_PATH} -O ${OUTPUT_PATH} -B ${BURN_IN} -s ${MCMC_ITER} -i ${MH_ITER}
