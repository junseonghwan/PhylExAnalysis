rm(list=ls())
#install.packages("PhylExR", repos=NULL, type="source")
library(PhylExR)

DATA_PATH <- "../data/HGSOC_SS3/"
BULK_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "sc.txt", sep="/")

# Load the data for evaluation.
dat <- read.table(BULK_PATH, header=T, sep="\t")
sc <- read.table(SC_PATH, header=T, as.is = TRUE)

hp <- EstimateHyperParameters(dat$ID, sc)
names(hp)
write.table(hp, paste(DATA_PATH, "sc_hp.txt", sep="/"), sep="\t", quote=F, col.names = T, row.names = F)
