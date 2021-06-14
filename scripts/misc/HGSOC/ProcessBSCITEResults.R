library(Rcpp)

source("R/EvaluationFunctions.R")
source("R/FindClones.R")
sourceCpp("src/rcpp_hello_world.cpp")

#data_path <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/total_cn/"
data_path <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/"

ssms <- read.table(paste(data_path, "bulk.txt", sep="/"), header=T, as.is=TRUE)
vafs <- ssms$b/ssms$d
mutation_count <- length(vafs)

bscite_output <- paste(data_path, "B-SCITE", "bscite.matrices", sep="/")
if (file.exists(bscite_output)) {
    bscite_clones <- GetClones(vafs, bscite_output)
    df <- data.frame(ID=ssms$ID, Cluster=bscite_clones, VAF=vafs)
    write.table(data.frame(ID=ssms$ID, CloneName=bscite_clones), paste(data_path, "/B-SCITE/results.txt", sep=""), quote=F, row.names=F)
}
