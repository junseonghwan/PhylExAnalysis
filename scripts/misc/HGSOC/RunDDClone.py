import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "./RunDDClone.sh" 
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('data_path', help="Path to data.", type=str)
parser.add_argument('single_cell_matrix_path', help="Path to single cell mutation matrix.", type=str)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 20000)
args = parser.parse_args()

DATA_PATH = args.data_path

seed = np.random.randint(100000000)

run_command = PATH_TO_EXECUTABLE + " "
run_command += str(seed) + " "
run_command += DATA_PATH + " "
run_command += args.single_cell_matrix_path + " "
run_command += str(args.mcmc)
os.system(run_command)
#print(run_command)
