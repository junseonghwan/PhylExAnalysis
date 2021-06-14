import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_EXECUTABLE  =  "./RunCanopy.sh" 
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('data_path', help="Path to the data.", type=str)
# parser.add_argument('--sbegin', help="Sim begin.", type=int, default = 0)
# parser.add_argument('--send', help="Sim end (not inclusive).", type=int, default = 10)
# parser.add_argument('--rbegin', help="Rep begin.", type=int, default = 0)
# parser.add_argument('--rend', help="Rep end (not inclusive).", type=int, default = 10)
parser.add_argument('--kbegin', help="Minimum number of clones.", type=int, default = 3)
parser.add_argument('--kend', help="Minimum number of clones.", type=int, default = 12)
args = parser.parse_args()

rep_seed = np.random.randint(100000000)

run_command = PATH_TO_EXECUTABLE + " "
run_command += str(rep_seed) + " "
run_command += args.data_path + " "
run_command += str(args.kbegin) + " "
run_command += str(args.kend)
os.system(run_command)
#print(run_command)
