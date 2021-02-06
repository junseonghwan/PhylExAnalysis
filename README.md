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

File for estimating hyper parameters from the single cell data is provided in `Rscripts/HGSOCEstimateHyperParams.R`.
Run this script to generate `data/HGSOC_sc_hp.txt`.

The scripts for running PhylEx is provided in `scripts/RunInference.py`.
We recommend to use 4 chains and the number of MCMC iterations should correspond to the number of SNVs. 
Typically, should use 10,000 or more iterations when the number of SNVs in the range of hundreds.
