import argparse
import datetime
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

PATH_TO_PROCESSING_SCRIPT  =  "./RunBSCITE.sh" # path to B-SCITE executable
PATH_TO_EXECUTABLE  =  "/Users/seonghwanjun/ScRNACloneEvaluation/B-SCITE/src/bscite.exe" # path to B-SCITE executable
parser = argparse.ArgumentParser(description='Script to batch run simulation studies.')
parser.add_argument('data_path', help="Path to data.", type=str)
parser.add_argument('single_cell_matrix_path', help="Path to single cell mutation matrix.", type=str)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 20000)
parser.add_argument('-n', '--num_repeats', help="Number of repeats.", type=int, default = 4)
args = parser.parse_args()

DATA_PATH = args.data_path

fp = 0.01 # esimated false positive rate of SCS experiment
fn = 0.2     # estimated false negative rate of SCS experiment
r  = args.num_repeats       # number of repeats
l  = args.mcmc    # number of loops

run_command = PATH_TO_PROCESSING_SCRIPT + " "
run_command += DATA_PATH + " "
run_command += args.single_cell_matrix_path
os.system(run_command)

OUTPUT_FILES_PREFIX = os.path.join(DATA_PATH, "B-SCITE/bscite") # this is the prefix used for the three output files (.matrices, .gv and .newick) given as the output of B-SCITE
SCFile   = os.path.join(DATA_PATH, "B-SCITE/bscite.SC") # path input SC File
bulkFile = os.path.join(DATA_PATH, "B-SCITE/bscite.bulk")     # path input bulk file

cell_count = np.loadtxt(os.path.join(DATA_PATH, "B-SCITE", "cell_count.txt"), dtype=int, ndmin=1)[0]
snv_count = np.loadtxt(os.path.join(DATA_PATH, "B-SCITE", "snv_count.txt"), dtype=int, ndmin=1)[0]

run_command = ""
run_command += PATH_TO_EXECUTABLE + " "
run_command += "-i "   + SCFile   + " "
run_command += "-bulkFileLocation " + bulkFile + " "
run_command += "-n "  + str(snv_count)  + " "
run_command += "-m "  + str(cell_count)  + " "
run_command += "-fd " + str(fp) + " "
run_command += "-ad " + str(fn) + " "
run_command += "-r "  + str(r)  + " "
run_command += "-l "  + str(l)  + " "
run_command += "-o "  + OUTPUT_FILES_PREFIX + " "

start = datetime.datetime.now()
os.system(run_command)
end = datetime.datetime.now()
result = end - start
exec_time_in_seconds = result.total_seconds()
with open(os.path.join(DATA_PATH, "B-SCITE/timing.txt"), "w+") as f:
	f.write(str(exec_time_in_seconds) + " seconds\n")
