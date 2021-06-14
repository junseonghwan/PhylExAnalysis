import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "RunDDClone.sh" 
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('sim_path', help="Path to simulation case.", type=str)
parser.add_argument('-t', '--threshold', help="Threshold to use for calling variants.", type=int, default = 2)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 20)
parser.add_argument('--cbegin', help="Case begin. For ddClone the first case is 1 since case 0 does not have any scRNA-seq data.", type=int, default = 1)
parser.add_argument('--cend', help="Case end (not inclusive).", type=int, default = 4)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 20000)
args = parser.parse_args()

SIM_PATH = args.sim_path
REP_BEGIN = args.rbegin
REP_END = args.rend
CASE_BEGIN = args.cbegin
CASE_END = args.cend
SC_VAR_THRESHOLD = args.threshold

for rep_no in range(REP_BEGIN, REP_END):
	REP_PATH = os.path.join(SIM_PATH, "rep" + str(rep_no))
	for case_no in range(CASE_BEGIN, CASE_END):
		CASE_PATH = os.path.join(REP_PATH, "case" + str(case_no))
		seed = np.random.randint(100000000)

		run_command = "sbatch "
		run_command += PATH_TO_EXECUTABLE + " "
		run_command += str(seed) + " "
		run_command += CASE_PATH + " "
		run_command += str(args.mcmc) + " "
		run_command += str(SC_VAR_THRESHOLD)
		#os.system(run_command)
		print(run_command)
