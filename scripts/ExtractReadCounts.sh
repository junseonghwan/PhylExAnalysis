#!/bin/bash -l

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J ExtractSingleCellReads
#SBATCH --mem 4G

echo "$(date) Running on: $(hostname)"

module load R/3.6.3-nsc1-gcc-7.3.0
wait

LOCI_FILE=$1
BAM_FILE=$2
OUTPUT_FILE=$3
R_OUTPUT_PATH=$4

echo $LOCI_FILE
echo $BAM_FILE
echo $OUTPUT_FILE
echo $R_OUTPUT_PATH

Rscript --vanilla \
		Rscripts/ProcessSingleCellBAM.R \
		$LOCI_FILE $BAM_FILE $OUTPUT_FILE > $R_OUTPUT_PATH
wait

echo "$(date) Finished extracting reads from scRNA-seq."
