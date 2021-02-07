rm(list=ls())
#install.packages("PhylExR", repos=NULL, type="source")
library(PhylExR)

DATA_PATH <- "~/PhylExAnalysis/data/"
BULK_PATH <- paste(DATA_PATH, "HGSOC_bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "HGSOC_sc.txt", sep="/")

# Load the data for evaluation.
dat <- read.table(BULK_PATH, header=T, sep="\t")
sc <- read.table(SC_PATH, header=T, as.is = TRUE)

hp <- EstimateHyperParameters(dat$ID, sc)
names(hp)
write.table(hp, paste(DATA_PATH, "HGSOC_sc_hp.txt", sep="/"), quote=F, col.names = T, row.names = F)
