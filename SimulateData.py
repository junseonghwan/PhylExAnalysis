import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

parser = argparse.ArgumentParser(description='Script to generate simulation data.')
parser.add_argument('seed', help="Seed for random.", type=int)
parser.add_argument('data_path', help="Path to store simulation data.", type=str)
parser.add_argument('-c', '--ncells', help="Number of cells to use, separated by comma.", type=str, default="0,100,200,400")
parser.add_argument('-r', '--nreps', help="Number of replicates.", type=int, default = 20)
parser.add_argument('-b', '--nbranches', help="Number of branches.", type=int, default = 2)
parser.add_argument('-d', '--max_depth', help="Max depth of the tree to simulate.", type=int, default = 3)
parser.add_argument('-s', '--nsites', help="Number of loci.", type=int, default = 100)
parser.add_argument('-R', '--nregions', help="Number of regions.", type=int, default = 1)
parser.add_argument('-f', '--min_cf', help="Min cell fraction.", type=float, default = 0.05)
parser.add_argument('--birth_rate', help="Birth rate.", type=float, default = 0)
parser.add_argument('--death_rate', help="Death rate.", type=float, default = 0)
parser.add_argument('--max_cn', help="Maximum copy number.", type=int, default = 10)
parser.add_argument('--variant_cp_prob', help="Probability of a copy being variant for CN simulation.", type=float, default = 0.5)
parser.add_argument('--var_clonal_cn_prob', help="Probability distn for variant copy number when using clonal copy number simulator. This value will be ignored if birth_rate > 0 is specified.", type=str, default ="0,1,0,0")
parser.add_argument('--ref_clonal_cn_prob', help="Probability distn for reference copy number when using clonal copy number simulator. This value will be ignored if birth_rate > 0 is specified.", type=str, default ="0,1,0,0")
parser.add_argument('--bursty_prob', help="The rate at which a site in a single cell data is bursty.", type=float, default = 0.8)
parser.add_argument('--dropout_rate', help="Dropout rate.", type=float, default = 0.8)
parser.add_argument('--randomize_cf', help="Randomize cell franction.", type=int, default = 0)
parser.add_argument('--randomize_branching', help="Randomize number of branches per node.", type=int, default = 0)
parser.add_argument('--randomize_dropout', help="Randomize dropout.", type=int, default = 1)
parser.add_argument('--bulk_mean_depth', help="Mean depth for bulk.", type=int, default = 1000)
parser.add_argument('--sc_mean_depth', help="Mean depth for scRNA-seq.", type=int, default = 20)
parser.add_argument('--seq_err', help="Sequencing error.", type=float, default = 0.01)
parser.add_argument('--sc_bursty_alpha', help="Alpha param for dropout distn.", type=float, default = 0.01)
parser.add_argument('--sc_bursty_beta', help="Beta param for dropout distn.", type=float, default = 0.01)
parser.add_argument('--beta_binomial_hp_max', help="Max value forBeta-Binomial hyperparam.", type=float, default = 10)
parser.add_argument('--bulk_sc_sparsity', help="Proportion of bulk sites to have single cell coverage.", type=float, default = 1.0)
parser.add_argument('--local', help="Run it locally.", action="store_true")
args = parser.parse_args()

if args.local:
	PATH_TO_EXECUTABLE  =  "PhylEx/simul" 
else:
	PATH_TO_EXECUTABLE  =  "/proj/sc_ml/users/x_seoju/PhylExAnalysis/PhylEx/simul" 

ncells = args.ncells.split(",")
# Assign common seed across all cases to ensure that 
# the tree and bulk data is the same across cases for same replicates.
for rep_no in range(args.nreps):
	rep_path = os.path.join(args.data_path, "rep" + str(rep_no))
	rep_seed = np.random.randint(100000000)
	for i, ncell in enumerate(ncells):
		sim_path = os.path.join(rep_path, "case" + str(i))
		if not os.path.exists(sim_path):
			os.makedirs(sim_path)
		config_file_path = os.path.join(sim_path, "simul.config")

		f = open(config_file_path, "w")
		f.write("seed: " + str(rep_seed) + "\n")
		f.write("num_branches: " + str(args.nbranches) + "\n")
		f.write("max_depth: " + str(args.max_depth) + "\n")
		f.write("n_sites: " + str(args.nsites) + "\n")
		f.write("n_cells: " + ncell + "\n")
		f.write("n_regions: " + str(args.nregions) + "\n")
		f.write("min_cf: " + str(args.min_cf) + "\n")
		if args.birth_rate > 0 and args.death_rate > 0:
			f.write("birth_rate: " + str(args.birth_rate) + "\n")
			f.write("death_rate: " + str(args.death_rate) + "\n")
		else:
			f.write("var_allele_copy_prob: " + str(args.var_clonal_cn_prob) + "\n")
			f.write("ref_allele_copy_prob: " + str(args.ref_clonal_cn_prob) + "\n")
		f.write("max_cn: " + str(args.max_cn) + "\n")
		f.write("var_cp_prob: " + str(args.variant_cp_prob) + "\n")
		f.write("bursty_prob: " + str(args.bursty_prob) + "\n")
		f.write("dropout_rate: " + str(args.dropout_rate) + "\n")
		f.write("randomize_cf: " + str(args.randomize_cf) + "\n")
		f.write("randomize_branching: " + str(args.randomize_branching) + "\n")
		f.write("randomize_dropout: " + str(args.randomize_dropout) + "\n")
		f.write("bulk_mean_depth: " + str(args.bulk_mean_depth) + "\n")
		f.write("sc_mean_depth: " + str(args.sc_mean_depth) + "\n")
		f.write("seq_err: " + str(args.seq_err) + "\n")
		f.write("sc_bursty_alpha0: " + str(args.sc_bursty_alpha) + "\n")
		f.write("sc_bursty_beta0: " + str(args.sc_bursty_beta) + "\n")
		f.write("beta_binomial_hp_max: " + str(args.beta_binomial_hp_max) + "\n")
		f.write("snv_sc_sparsity: " + str(args.bulk_sc_sparsity) + "\n")
		f.write("output_path: " + sim_path + "\n")
		f.close()

		run_command = PATH_TO_EXECUTABLE + " -c " + config_file_path
		print(run_command)
		os.system(run_command)

