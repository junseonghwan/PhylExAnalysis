#!/usr/bin/env python

import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

np.random.seed(1)

# code to process command line arguments
parser = argparse.ArgumentParser(description='Script to process single cell BAMs.')
parser.add_argument('loci_file', help="Loci file name within the sample_path.", type=str)
parser.add_argument('sc_bam_path', help="Path to the scRNA BAM files.", type=str)
parser.add_argument('output_path', help="Path to output processed files.", type=str)
parser.add_argument('-s', '--suffix', help="Single cell BAM file suffix.", type=str, default=".bam")
args = parser.parse_args()

# This script is step 2 in the following data processing steps:
# 1. Process the VCF files to generate loci file: ExtractLoci.sh->Rscripts/ExtractLoci.R. Do it before using this script.
# 2. Process the single cell files at the identified loci from step 1. ExtractReadCounts.sh->Rscripts/ProcessSingleCellBAMScript.R
# 3. Combine single cell data and identify final set of loci.
# 4. Generate inputs for each of the methods to be used for comparison.

# 2. Process single cell BAMs.
if not os.path.exists(args.output_path):
	os.makedirs(args.output_path)
samples = glob.glob(args.sc_bam_path + "/*" + args.suffix)
for sample in samples:
	sample_name = os.path.basename(sample)
	cell_name = sample_name.split(args.suffix)[0]
	output_name = os.path.join(args.output_path, cell_name + ".reads.txt")
	r_output_path = os.path.join(args.output_path, cell_name + ".Rout")
	command = "sbatch ExtractReadCounts.sh " + args.loci_file + " " +  sample + " " + output_name + " " + r_output_path
	#print(command)
	os.system(command)
