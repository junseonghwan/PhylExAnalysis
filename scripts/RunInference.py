import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

parser = argparse.ArgumentParser(description='Script to run PhylEx jobs.')
parser.add_argument('seed', help="Seed for random.", type=int)
parser.add_argument('bulk_data_path', help="Path to bulk data.", type=str)
parser.add_argument('sc_data_path', help="Path to single cell data.", type=str)
parser.add_argument('sc_hp_data_path', help="Path to single cell hyper parameters.", type=str)
parser.add_argument('output_path', help="Path to output the results of PhylEx.", type=str)
parser.add_argument('-n', '--num_reps', help="Number of chains to run.", type=int, default = 4)
parser.add_argument('-m', '--mcmc', help="Number of MCMC iterations.", type=int, default = 2000)
parser.add_argument('-p', '--mh', help="Number of MH iterations within each MCMC iteration.", type=int, default = 2000)
parser.add_argument('-t', '--thinning', help="Thinning interval.", type=int, default = 20)
parser.add_argument('-u', '--burn_in', help="Burn in.", type=int, default = 500)
parser.add_argument('-o', '--overwrite', help="Overwrite main.config if it exists.", action="store_true")
#parser.add_argument('--geometric', help="Use geometric mean for single cell likelihood.", action="store_true")
parser.add_argument('--local', help="Run it locally.", action="store_true")
parser.add_argument('--bulk_only', help="Run it without single cells.", action="store_true")
args = parser.parse_args()

PATH_TO_EXECUTABLE  =  "RunInference.sh"
SEED = args.seed
OUTPUT_PATH = os.path.abspath(args.output_path)
N_MCMC_ITER = args.mcmc
N_MH_ITER = args.mh
THINNING = args.thinning
BURN_IN = args.burn_in

# Default values.
ALPHA0_MAX = 1
LAMBDA_MAX = 0.2
GAMMA_MAX = 0.5
ALPHA0_MIN = 0
LAMBDA_MIN = 0
GAMMA_MIN = 0
SEQ_ERR = 0.001
VAR_CP_PROB = 0.25
SC_BURSTY_ALPHA0 = 0.01
SC_BURSTY_BETA0 = 0.01
#GEOMETRIC_MEAN = (1 if args.geometric else 0)
GEOMETRIC_MEAN = 0

np.random.seed(SEED)

for rep_no in range(args.num_reps):
    # write main.config file to REP_PATH
    if args.bulk_only:
        REP_PATH = OUTPUT_PATH + "/tssb/chain" + str(rep_no) + "/"
    else:
        REP_PATH = OUTPUT_PATH + "/phylex/chain" + str(rep_no) + "/"
    if not os.path.exists(REP_PATH):
        os.makedirs(REP_PATH)
    config_file_path = os.path.join(REP_PATH, "main.config")
    if not os.path.exists(config_file_path) or args.overwrite:
        print("Writing a new main.config file for " + REP_PATH)
        rep_seed = np.random.randint(100000000)

        f = open(config_file_path, "w")
        f.write("seed: " + str(rep_seed) + "\n")
        f.write("bulk_data_path: " + os.path.abspath(args.bulk_data_path) + "\n")
        if not args.bulk_only:
            f.write("sc_rna_data_path: " + os.path.abspath(args.sc_data_path) + "\n")
            f.write("hyperparams_path: " + os.path.abspath(args.sc_hp_data_path) + "\n")
            output_path = os.path.join(REP_PATH, "")
        f.write("output_path: " + REP_PATH + "\n")
        f.write("n_mcmc_iter: " + str(N_MCMC_ITER) + "\n")
        f.write("n_mh_iter: " + str(N_MH_ITER) + "\n")
        f.write("thinning: " + str(THINNING) + "\n")
        f.write("burn_in: " + str(BURN_IN) + "\n")
        f.write("alpha0_max: " + str(ALPHA0_MAX) + "\n")
        f.write("lambda_max: " + str(LAMBDA_MAX) + "\n")
        f.write("gamma_max: " + str(GAMMA_MAX) + "\n")
        f.write("alpha0_min: " + str(ALPHA0_MIN) + "\n")
        f.write("lambda_min: " + str(LAMBDA_MIN) + "\n")
        f.write("gamma_min: " + str(GAMMA_MIN) + "\n")
        f.write("seq_err: " + str(SEQ_ERR) + "\n")
        f.write("var_cp_prob: " + str(VAR_CP_PROB) + "\n")
        #f.write("sc_dropout_alpha0: " + str(SC_DROPOUT_ALPHA0) + "\n")
        #f.write("sc_dropout_beta0: " + str(SC_DROPOUT_BETA0) + "\n")
        f.write("sc_bursty_alpha0: " + str(SC_BURSTY_ALPHA0) + "\n")
        f.write("sc_bursty_beta0: " + str(SC_BURSTY_BETA0) + "\n")
        f.write("geometric_mean: " + str(GEOMETRIC_MEAN) + "\n")
        f.close()

    if args.local:
        run_command = "./" + PATH_TO_EXECUTABLE + " "
    else:
        run_command = "sbatch "
        run_command += PATH_TO_EXECUTABLE + " "

    run_command += config_file_path
    #print(run_command)
    os.system(run_command)
