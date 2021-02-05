#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J PhylEx

echo "$(date) Run begins."

CONFIG_FILE=$1

cd ../PhylEx
./run -c ${CONFIG_FILE}

echo "$(date) Run finished."
