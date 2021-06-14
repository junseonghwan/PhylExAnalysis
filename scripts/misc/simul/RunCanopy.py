import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np
 
PATH_TO_EXECUTABLE  =  "RunCanopy.sh" 
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('sim_path', help="Path to simulation data.", type=str)
parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 20)
parser.add_argument('--cbegin', help="Case begin.", type=int, default = 0)
parser.add_argument('--cend', help="Case end (not inclusive).", type=int, default = 4)
parser.add_argument('--kbegin', help="Minimum number of clones.", type=int, default = 3)
parser.add_argument('--kend', help="Minimum number of clones.", type=int, default = 12)
args = parser.parse_args()

SIM_PATH = args.sim_path
REP_BEGIN = args.rbegin
REP_END = args.rend
CASE_BEGIN = args.cbegin
CASE_END = args.cend

for rep_no in range(REP_BEGIN, REP_END):
	REP_PATH = os.path.join(SIM_PATH, "rep" + str(rep_no))
	for case_no in range(CASE_BEGIN, CASE_END):
		CASE_PATH = os.path.join(REP_PATH, "case" + str(case_no))
		seed = np.random.randint(100000000)

		run_command = "sbatch "
		run_command += PATH_TO_EXECUTABLE + " "
		run_command += str(seed) + " "
		run_command += CASE_PATH + " "
		run_command += str(args.kbegin) + " "
		run_command += str(args.kend)
		#os.system(run_command)
		print(run_command)
