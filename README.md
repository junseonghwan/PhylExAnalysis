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

The scripts for running PhylEx is provided in `scripts/RunInference.py`.

Examples for the input data is in `data/`. Check `data/HGSOC_bulk.txt` for bulk data input and `data/HGSOC_sc.txt` for single cell data.

File for estimating hyper parameters from the single cell data is provided in `Rscripts/HGSOCEstimateHyperParams.R`.
