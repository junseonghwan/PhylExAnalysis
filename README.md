# PhylExAnalysis

```
cd PhylEx
git submodule init
git submodule update
mkdir build
cd build
cmake ..
make install
cd ../../
./PhylEx/run -c main.config
```

This will generate the posterior samples in the directory `_output/HGSOC`.

We will use the script in `Rscripts` to analyze and reproduce the results.