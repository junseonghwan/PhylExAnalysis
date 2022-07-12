rm(list=ls())
library(PhylExR)

library(biomaRt)
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

DATA_PATH <- "~/PhylExAnalysis/data/HGSOC_SS3/"
PHYLEX_OUTPUT_PATH <- "~/PhylExAnalysis/data/HGSOC_SS3/phylex/"
ANALYSIS_OUTPUT_PATH <- "~/PhylExAnalysis/_figures/HGSOC/"
FEATURE_COUNTS_PATH <- paste(DATA_PATH, "featureCounts.txt", sep="/")
BULK_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "sc.txt", sep="/")
SC_HP_PATH <- paste(DATA_PATH, "sc_hp.txt", sep="/")
GT_PATH <- paste(DATA_PATH, "gt.txt", sep="/")

if (!dir.exists(ANALYSIS_OUTPUT_PATH)) {
  dir.create(ANALYSIS_OUTPUT_PATH, recursive = T)
}

topK <- 1000
MIN_CELLS <- 10
MIN_READS <- 5
MIN_SNV_COUNT <- 1
chrs <- c(1:22, "X", "Y")
bursty_hp <- list(alpha=0.01, beta=0.01)
FDR_THRESHOLD <- 0.1
chains <- 0:19

# Load the data for evaluation.
dat <- read.table(BULK_PATH, header=T, sep="\t")
sc <- read.table(SC_PATH, header=T, as.is = TRUE)
sc$b <- sc$d - sc$a
sc_hp <- read.table(SC_HP_PATH, header=T)
fc <- read.table(FEATURE_COUNTS_PATH, header=T, row.names = 1, check.names = F)
valid_clones <- c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "A_B", "C_D", "A", "B", "C", "D", "E_F_G_H_I", "E_F", "E", "F")
gt <- read.table(GT_PATH, header=T, as.is = T)
validation_idx <- which(gt$CloneName %in% valid_clones)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Find the best replicate to use for evaluation.
best_chain <- FindBestRep(chain_paths = PHYLEX_OUTPUT_PATH, chains = chains)
best_chain_path <- paste(PHYLEX_OUTPUT_PATH, "/chain", best_chain, sep="")
datum2node <- read.table(paste(best_chain_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)

datum2node_ <- left_join(datum2node, gt)
datum2node_ <- datum2node_[,c("ID", "Node", "CloneName")]
datum2node_ <- datum2node_[validation_idx,]
tree <- datum2node_[order(nchar(datum2node_$Node), datum2node_$Node),]
table(tree$CloneName)
table(tree$Node)
datum2node_$Node <- as.character(datum2node_$Node)
datum2node.sorted <- datum2node_[order(datum2node_$Node),]
# We can map the inferred nodes to ground truth clones:
# 0_0: Ancestral
# 0_0_0: ABCD
# 0_0_1: EF/EFGHI
# 0_0_0_0: CD
# 0_0_0_0_0: C
# Note: node names may vary depending on the number of chains and length of MCMC iterations.
unique(tree[,-1])

cell.df <- AssignCellsBursty(sc, datum2node_, bursty_hp, sc_hp, include_normal_clone = FALSE)
cell.df$Node <- as.character(cell.df$Node)
table(cell.df$Node)

cell.df$Clone <- cell.df$Node
cell.df[cell.df$Clone == "0_0",]$Clone <- "Ancestral"
cell.df[startsWith(cell.df$Clone, "0_0_1"),]$Clone <- "EF"
cell.df[cell.df$Clone == "0_0_0",]$Clone <- "ABCD"
cell.df[cell.df$Clone == "0_0_0_0",]$Clone <- "CD"
cell.df[cell.df$Clone == "0_0_0_0_0",]$Clone <- "C"
cell.df$Clone <- factor(cell.df$Clone, levels = c("Ancestral", "ABCD", "CD", "C", "EF"))
cell_order <- cell.df[order(cell.df$Clone), "Cell"]

sc_join <- left_join(sc, cell.df, by = "Cell")
sc_join <- subset(sc_join, ID %in% datum2node.sorted$ID)
sc_join$Cell <- factor(sc_join$Cell, levels = cell_order)
sc_join$ID <- factor(sc_join$ID, levels = datum2node.sorted$ID)

base_size <- 11
p <- ggplot(subset(sc_join, b >0), aes(ID, Cell, fill=Clone)) + geom_tile(colour = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cell")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.text = element_text(size=7, face="bold")) + guides(fill=guide_legend(title="Clone", size = 7))
p
ggsave(p, filename = paste(ANALYSIS_OUTPUT_PATH, "co-clustering.pdf", sep="/"), height = 8, width = 3.5, units = "in")

pl1 <- PlotTotalCounts(sc)
ggsave(pl1, filename = paste(ANALYSIS_OUTPUT_PATH, "raw_expression.pdf", sep="/"), height = 8, width = 3.5, units = "in")
