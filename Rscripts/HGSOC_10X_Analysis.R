library(PhylExR)

library(edgeR)
library(dplyr)
library(ggplot2)
library(Glimma)
library(limma)
library(mclust)
library(Rtsne)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(zinbwave)

#DATA_PATH <- "~/PhylExAnalysis/data/"
DATA_PATH <- "~/data/cell-line/10X2/"
#PHYLEX_OUTPUT_PATH <- "~/PhylExAnalysis/_output/HGSOC_10X/phylex"
PHYLEX_OUTPUT_PATH <- "~/data/cell-line/10X2/phylex/"
ANALYSIS_OUTPUT_PATH <- "~/data/cell-line/10X2/analysis"
BULK_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "sc.txt", sep="/")
SC_HP_PATH <- paste(DATA_PATH, "sc_hp.txt", sep="/")
GT_PATH <- paste(DATA_PATH, "gt.txt", sep="/")

# BULK_PATH <- paste(DATA_PATH, "HGSOC_10X_bulk.txt", sep="/")
# SC_PATH <- paste(DATA_PATH, "HGSOC_10X_sc.txt", sep="/")
# SC_HP_PATH <- paste(DATA_PATH, "HGSOC_10X_sc_hp.txt", sep="/")
# GT_PATH <- paste(DATA_PATH, "HGSOC_10X_gt.txt", sep="/")

if (!dir.exists(ANALYSIS_OUTPUT_PATH)) {
    dir.create(ANALYSIS_OUTPUT_PATH, recursive = T)
}

MIN_SNV_COUNT <- 1
chrs <- c(1:22, "X", "Y")
bursty_hp <- list(alpha=0.01, beta=0.01)
chains <- 0:3

# Load the data for evaluation.
dat <- read.table(BULK_PATH, header=T, sep="\t")
sc <- read.table(SC_PATH, header=T, as.is = TRUE)
sc_hp <- read.table(SC_HP_PATH, header=T)
valid_clones <- c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "A_B", "C_D", "A", "B", "C", "D", "E_F_G_H_I", "E_F", "E", "F")
gt <- read.table(GT_PATH, header=T, as.is = T)
validation_idx <- which(gt$CloneName %in% valid_clones)

# Find the best replicate to use for evaluation
best_chain <- FindBestRep(chain_paths = PHYLEX_OUTPUT_PATH, chains = chains)
best_chain_path <- paste(PHYLEX_OUTPUT_PATH, "/chain", best_chain, sep="")

datum2node <- read.table(paste(best_chain_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)
datum2node_ <- CollapseClones(datum2node, MIN_SNV_COUNT)
table(datum2node_$Node)

temp <- left_join(datum2node_, gt)
temp <- temp[,c("ID", "Node", "CloneName")]
temp <- temp[validation_idx,]
head(temp)
tree <- temp[order(nchar(temp$Node), temp$Node),]
head(tree)
for (node in unique(tree$Node))
{
    print(paste("Node:", node))
    print(table(subset(tree, Node == node)$CloneName))
}
table(tree$Node)
table(tree$CloneName)
subset(tree, CloneName == "A_B")
subset(tree, CloneName == "A_B_C_D")
subset(tree, CloneName == "A_B_C_D_E_F_G_H_I")
subset(tree, CloneName == "E")
subset(tree, CloneName == "E_F")
subset(tree, CloneName == "E_F_G_H_I")

