# PhylExAnalysis

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

Examples for the input data is in `data/`. Check `data/HGSOC_bulk.txt` for bulk data input and `data/HGSOC_sc.txt` for single cell data.

## Example data processing step
Our data processing step was to first process VCF file to generate the bulk data. An Rscript can be found in `Rscripts/PrepareBulkData.R`. This script outputs `bulk.txt` and `loci.txt`. The bulk data is initially filled in with major and minor copy numbers equal to 1. Separate script is needed to extract copy number profiles from software used for copy number analysis.

We assume that scRNA-seq data has been pre-processed, such as alignment, and available as a BAM file. We have an Rscript that takes `loci.txt` and location to BAM file and extracts read data. We recommend running this script on all cells and use script `Rscripts/CombineSingleCellData.R` to generate `sc.txt` and `sc_hp.txt` file.

If separate process was followed to generate the bulk and single cell data files, we recommend to estimate the hyper parameters for the single cell data using `Rscripts/HGSOCEstimateHyperParams.R`.

## Run PhylEx
The scripts for running PhylEx is provided in `scripts/RunInference.py`.
We recommend to use 4 chains and the number of MCMC iterations should correspond to the number of SNVs. 
We used should use 20,000 iterations for HER2+ data, which had 432 SNVs and 369 cells.
We used 10,000 iterations for HGSOC data, which had 67 SNVs and 360 cells.

