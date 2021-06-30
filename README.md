# PhylExAnalysis

## Install PhylExR

We provide common functions used for data processing and analysis as an R package `PhylExR`. We recommend to install this package first.
From the command line,
```
R CMD build PhylExR
```
This command should create `PhylExR_1.0.tar.gz` file.
The following command will install the package.
```
R CMD install PhylExR_1.0.tar.gz
```

## Install PhylEx
Install dependencies for [PhylEx](https://github.com/junseonghwan/PhylEx).
```
git submodule init
git submodule update
cd PhylEx
mkdir build
cd build
cmake ..
make install
```

Installation should only take a few minutes.

## PhylEx input data

The bulk data has 5 columns separated by tabs: `ID`, `b`, `d`, `major_cn`, `minor_cn`.
The columns `b` and `d` are integer counts of the number of reads mapping to the variants and the total depth at a loci identified by `ID`.
If the data consists of multiple regions, the columns `b`, `d`, `major_cn`, `minor_cn` should contain the counts separated by comma. The user may also specify 9 columns including CHR, POS, REF, and ALT alleles in the following order: `ID`, `CHR`, `POS`, `REF`, `ALT`, `b`, `d`, `major_cn`, `minor_cn`. The chromosome, position, reference allele, and alternate allele columns are simply ignored by PhylEx but having them in one data file with the bulk data may be desirable for further analysis.

The single cell data contains 4 columns: `ID`, `Cell`, `a`, `d`.
The `ID` column corresponds to the `ID` column in the bulk data file. The `Cell` column identifies the cell name/ID. The column `a` is the number of reads mapping to the reference allele at `ID` for cell `Cell` and `d` is the total number of reads mapping to the loci identified by `ID`.

Finally, there is a file containing the hyperparameters to be used in the inference. We recommend to auto-generate this file using the scripts provided (see the subsequent sections for details).

The data files accepted by PhylEx are provided in [Zenodo](https://zenodo.org/record/4950446#.YNaC9BNKj6A). Download HER2_POS_SS3.tar.gz to `data/`, extract it and see `bulk.txt`, `sc.txt`, and `sc_hp.txt`.

## Real biological data processing

We have provided scripts to help generate the input files and the configuration file needed to run PhylEx in this section.

### Bulk data

Relevant script:
- `Rscripts/PrepareBulkData.R`

We first processes the VCF file to generate the bulk data. An Rscript `Rscripts/PrepareBulkData.R` is provided as a reference. This script accepts three arguments: `VCF_PATH`, `CNV_PATH`, and `OUTPUT_PATH`. The `VCF_PATH` is the path to `.vcf` file generated from a variant caller such as MuTect2. The `CNV_PATH` is to the path containing the outputs of TitanCNA. We will look for the optimal copy number profile selected by TitanCNA in the path: `${CNV_PATH}/cna/results/titan/hmm/optimalClusterSolution.txt`.
This script outputs `bulk.txt` and `loci.txt`. 

### scRNA-data

Relevant scripts:
- `scripts/ProcessSingleCells.py`
- `scripts/ExtractReadCounts.sh`
- `scripts/ProcessSingleCellBAM.R`
- `Rscripts/CombineSingleCellData.R`

We assume that scRNA-seq data has been put through the standard pre-processing, i.e., alignment, procesing of UMI barcodes, etc, and are available as one BAM file for each cell. We have an Rscript that takes `loci.txt` outputted from the previous step and location to the scRNA-seq BAM files to extract the read data (`scripts/ProcessSingleCellBAM.R`). We recommend batch run this R script on all BAM files using `scripts/ProcessSingleCells.py`. Then, use `Rscripts/CombineSingleCellData.R` to generate `sc.txt` and `sc_hp.txt` file.

### Estimating the hyperparameters

Relevant scripts:
- `Rscripts/HGSOCEstimateHyperParams.R`

Generating the bulk and scRNA data files can be accomplished without the above scripts as long as the final data meets the format described above. In such a case, we recommend to generate the hyper parameter file for the single cell data by following the example provided in `Rscripts/HGSOCEstimateHyperParams.R`. Essentially, the function `PhylExR::EstimateHyperParameters` performs the task of estimating the hyperparameters from the scRNA-seq data.

## Running PhylEx

Relevant scripts:
- `scripts/RunInference.py`
- `scripts/RunInference.sh`

Once the input data are generated, we recommend to use a script to `scripts/RunInference.py`, which generates the configuration file and runs PhylEx.
We recommend to use 4 chains and the number of MCMC iterations should depend on the number of SNVs. 
As a point of reference, we used 20,000 iterations for HER2+ data, which had 432 SNVs and 369 cells and 10,000 iterations for HGSOC data, which had 67 SNVs and 360 cells.
The run time for HGSOC should take about 10-30 minutes depending on the hardware using a single core. 
The run time is mostly affected by the number of SNVs under consideration; hence, the run time for HER2+ can take about 12 hours up to 24 hours.
In general, start with a smaller number of iterations and monitor the log likelihood.

**We plan to add a feature to start a run from a previous run to facilitate this process.**


## Reproducibility
The figures and tables presented in [PhylEx paper](https://www.biorxiv.org/content/10.1101/2021.02.16.431009v1) can be reproduced from the results deposited to [Zenodo](https://zenodo.org/record/4950446). We have added resource files needed in the analysis to `data/` directory.

1. Download the results files for HGSOC and HER2+: `HGSOC_SS3.tar.gz`, `HER2_POS_SS3.tar.gz`. Uncompress the files under `data/`.
2. Run `Rscripts/HER2_DGE.R` to reproduce figures for HER2+ breast cancer.
    - For HER2+ analysis, we downloaded [PanCancer pathway genes](https://www.nanostring.com/products/ncounter-assays-panels/panel-selection-tool/) to be downloaded from NanoString (Seattle, WA) -- made available as `data/LBL-10025_nCounter_PanCancer_Human_Pathways_Panel_Gene_List.xlsx`.
    - It also requires exons, we have a pre-processed the file as `data/exons.bed`.
    - Finally, we use MSigDB's Hallmark geneset available as `data/human_H_v5p2.rdata`.
3. The phylo-phenotypic analysis for ovarian cancer cell-line data can be reproduced by running `Rscripts/HGSOC_SS3_DGE.R` and `Rscripts/HGSOC_SS3_coclustering_plot.R`. 
    - For pathway analysis, we use MSigDB's ontology geneset available as `data/human_c5_v5p2.rdata`.
4. Different simulation cases and theirs results are also deposited to Zenodo. For example, `binary_cn.tar.gz` contains the data files and the result files for binary clonal tree and bulk data simulated with birth-death process acting on the copy numbers. The result files also include outputs from running other softwares. 

For full reproducibility, copy number and variant call results from the bulk data are also deposited to [Zenodo](https://zenodo.org/record/4950446). For example, `HER2_POS.tar.gz` and `HGSOC.tar.gz` files include copy number and variant calls from the bulk data. 

## Simulation data generation
The code for simulating the data is provided in the script: `scripts/SimulateData.py`.
The simulation options can be separated into three parts: 1) simulation conditions, 2) bulk data, 3) single cell data.

The simulation conditions are number of cells, number of replicates, tree size and minimum clone fraction:

- `-c`: Number of cells separated by comma. Default: `0,100,200,400`.
- `-r`: Number of replicates. Default: `20`.
- `-b`: Number of branches. Default: `2` for binary tree.
- `-d`: Depth of tree. Default: `3`.
- `-f`: Minimum clone fraction (cellular prevalence). Default: `0.05`.
- `--randomize_cf`: Randomize simulation of cellular fraction. Default: `0` for deterministical generation of cellular fraction where the cell fraction is halved at each node of the parent's remaining cellular fraction.
- `--randomize_branching`: Randomize the number of branches. Default: `0`. Setting this to `1` will randomly select the number of children for each node up to the value specified using `-d`.

The bulk data options control the number of SNVs, number of regions, and copy number simulation:

- `-s`: Number of sites (SNVs). Default: `100`.
- `-R`: Number of regions for simulating multi-regional bulk data. Default: `1` for single region.
- `--bulk_mean_depth`: The depth for the bulk data. Default: `1,000`.
- `-birth_rate`: Birth rate. Default: `0`.
- `-death_rate`: Death rate. Default: `0`.
- `-max_cn`: Maximum (total) copy number. Default: `10`.
- `-variant_cp_prob`: Probability of a copy being a variant. Default: `0.5`.
- `--var_clonal_cn_prob`: An alternative approach to generating copy number profiles. For variant copy number simulation. Default `0,1,0,0` for 1 copy of variant.
- `--ref_clonal_cn_prob`: An alternative approach to generating copy number profiles. For reference copy number simulation. Default `0,1,0,0` for 1 copy of reference copy.

The single cell data options are bursty probability controlling the probability of mono-allelic expression, dropout

- `--bursty_prob`: Probability of mono-allelic expression. Default: `0.8`.
- `--dropout_rate`: Number of sites to drop in the single cell expression. Default `0.8`, 80% of the number of sites will be dropped or 20% will be expressed.
- `--randomize_dropout`: Randomize dropout. The exact dropout rate will be randomly select from the value specified by `--dropout_rate` to `1`. Default: `1`.
- `--sc_mean_depth`: Mean depth for single cell data. Default: `20`.
- `--sc_bursty_alpha`: Hyperparameter $\alpha$ for mono-allelic distribution (Beta distribution). Default: `0.01`.
- `--sc_bursty_beta`: Hyperparameter $\beta$ for mono-allelic distribution (Beta distribution). Default: `0.01`.
- `--beta_binomial_hp_max`: The maximum value to use in quantifying the bi-allelic distribution across sites. Default: `10`.

There are options that relate to both the single cell and bulk data.
- `--snv_sc_sparsity`: The proportion of sites (SNVs) to have single cell coverage, specify a value int he range `[0, 1]`. Default: `1.0`, all sites are candidates to be covered. For a value below `1`, we will randomly select sites to never be expressed in the single cell simulation.
- `--seq_err`: The error probability for bulk data and for the error distribution of the single cell data. Default: `0.01`.

There is an option `--local` if it is specified, it will use the path to the executable under `PhylEx/simul`. If it is not specified, it uses an alternative path; you will need to modify `PATH_TO_EXECUTABLE` in the script to point to the correct path.

### Example

To simulate multifurcating tree with maximum of 4 branches and CN simulation with birth and death rates of `1` and `0.2` to be run locally, the below command can be used:
```
python scripts/SimulateData.py 123 simul/ -b 4 --randomize_branching 1 --birth_rate 1 --death_rate 0.2 --local 
```
The default values are used for the number of replicates (`20`), number of cells (`0,100,200,400`), and etc.

For binary tree,
```
python scripts/SimulateData.py 123 simul/ -b 2 --birth_rate 1 --death_rate 0.2 --local 
```

The simulated data generated using the script and analyzed in [PhylEx paper](https://www.biorxiv.org/content/10.1101/2021.02.16.431009v1) are deposited on [Zenodo](https://doi.org/10.5281/zenodo.4533670).

### Running PhylEx on simulated data
The script `scripts/RunSimulationStudy.py` contains the code for batch running the simulation studies.
It assumes that there are 20 replicates and 4 cases involving the number of cells. To run simulation study for specific replicates, use the options `--rbegin` and `--rend`. And to run it for specific cases, use `--cbegin` and `--cend`.

The number of MCMC iterations for TSSB and MH iterations for sampling the cellular prevalences can be specified using the options `-m` and `-p`. The thinning and burn in can be controlled using `-t` and `-u`. The option `--nchains` specifies how many parallel chains to deploy for each simulation rep+case.

Finally, we have option `--local` for doing a test run before a batch run.

This python script executes `scripts/RunSimulationStudy.sh`, which runs the script to generate the necessary input for PhylEx, which is found in `scripts/GenerateInputForSimulatedStudy.R` before running PhylEx. Be sure to check these scripts and set the correct path to the executables depending on where you run them from and your installation path.

Example:
```
python RunSimulationStudy.py 123 simul/ --nchains 4 -m 2000 -p 2000
```


## Comparison to other softwares
In the [PhylEx paper](https://www.biorxiv.org/content/10.1101/2021.02.16.431009v1), we have compared our method against [PhyloWGS](https://github.com/morrislab/phylowgs), [Canopy](https://github.com/yuchaojiang/Canopy), [ddClone](https://github.com/sohrabsa/ddclone), and [B-SCITE](https://github.com/smalikic/B-SCITE). The scripts are found under `scripts/misc/`. 

These scripts were run locally/NSC Tetralith and will need work to get them running on your local machine. We are working on providing the scripts within a container/environment for maximal reproducibility.

For evaluation, we wrote Jupyter notebooks, the files can be found under `notebook/HGSOC` for evaluation concerning HGSOC cell-line and `notebook/simul` for evaluation concerning simulated data.

## Contact
For questions, please email: sjun2@fredhutch.org. Please report bugs using Github issue.
