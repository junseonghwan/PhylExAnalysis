# script to launch PhyloWGS experiments
import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

parser = argparse.ArgumentParser(description='Script to batch submit simulation studies.')
parser.add_argument('seed', help="Seed for random.", type=int)
parser.add_argument('data_path', help="Path to the data.", type=str)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 5000)
parser.add_argument('-p', '--mh', help="Number of MH iterations.", type=int, default = 2000)
parser.add_argument('-u', '--burn_in', help="Burn_in.", type=int, default = 1000)
parser.add_argument('--local', help="Run it locally.", action="store_true")
parser.add_argument('--generate_input', help="Generate input.", action="store_true")
args = parser.parse_args()

if args.local:
    PATH_TO_EXECUTABLE  =  "RunPhyloWGSLocally.sh" # path to executable
else:
    PATH_TO_EXECUTABLE  =  "RunPhyloWGS.sh" # path to executable
SEED = args.seed
DATA_PATH = args.data_path
N_MCMC_ITER = args.mcmc
N_MH_ITER = args.mh
BURN_IN = args.burn_in

np.random.seed(SEED)

seed = np.random.randint(100000000)
pwgs_bulk_path = os.path.join(DATA_PATH, "pwgs_snv.txt")
pwgs_cnv_path = os.path.join(DATA_PATH, "pwgs_cnv.txt")
output_path = os.path.join(DATA_PATH, "phylowgs/")

if not os.path.exists(output_path):
    os.makedirs(output_path)

if args.generate_input:
	run_command = "Rscript --vanilla Rscripts/GenerateInputForPhyloWGS.R " + DATA_PATH
	os.system(run_command)

if args.local:
    run_command = "./" + PATH_TO_EXECUTABLE + " "
else:
    run_command = "sbatch "
    run_command += PATH_TO_EXECUTABLE + " "

run_command += DATA_PATH + " "
run_command += pwgs_bulk_path + " "
run_command += pwgs_cnv_path + " "
run_command += output_path + " "
run_command += str(BURN_IN) + " "
run_command += str(N_MCMC_ITER) + " "
run_command += str(N_MH_ITER) + " "
run_command += str(seed) + " "

#print(run_command)
os.system(run_command)
