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
A script analyzing the ovarian cancer cell line data is provided in `Rscripts/HGSOCAnalysis.R`, which requires `PhylExR` package to be installed.

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

Examples for the input data is in `data/`. 
Check `data/HGSOC_bulk.txt` for bulk data input and `data/HGSOC_sc.txt` for single cell data.

## Example data processing step
Our data processing step processes the VCF file to generate the bulk data. An Rscript can be found in `Rscripts/PrepareBulkData.R`. This script outputs `bulk.txt` and `loci.txt`. The bulk data is initially filled in with major and minor copy numbers equal to 1. Separate script is needed to extract copy number profiles from software used for copy number analysis.

We assume that scRNA-seq data has been pre-processed, such as alignment and procesing UMI barcodes, and available as one BAM file for each cell. We have an Rscript that takes `loci.txt` and location to BAM file and extracts read data (see `scripts/ProcessSingleCellBAM.R`). We recommend batch run this script on all BAM files, i.e., cells (see `ProcessSingleCells.py`). Then, use `Rscripts/CombineSingleCellData.R` to generate `sc.txt` and `sc_hp.txt` file.

If different processes were followed to generate the bulk and single cell data files, we recommend to estimate the hyper parameters for the single cell data using `Rscripts/HGSOCEstimateHyperParams.R`.

## Run PhylEx
The scripts for running PhylEx is provided in `scripts/RunInference.py`.
We recommend to use 4 chains and the number of MCMC iterations should correspond to the number of SNVs. 
We used should use 20,000 iterations for HER2+ data, which had 432 SNVs and 369 cells.
We used 10,000 iterations for HGSOC data, which had 67 SNVs and 360 cells.

## Contact
For questions, please email: sjun2 at fredhutch.org.
