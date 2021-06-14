#!/bin/bash

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

python /Users/seonghwanjun/ScRNACloneEvaluation/phylowgs/multievolve.py --num-chains 4 --ssms ${PWGS_SSM_PATH} \
	--cnvs ${PWGS_CNV_PATH} -O ${OUTPUT_PATH} -B ${BURN_IN} -s ${MCMC_ITER} -i ${MH_ITER}
