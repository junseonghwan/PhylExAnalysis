args = commandArgs(trailingOnly=TRUE)

print(args)

REP_PATH <- args[1]
THRESHOLD <- as.numeric(args[2])
#THRESHOLD <- 1
#REP_PATH <- "/Users/seonghwanjun/data/simul/quadternary_multiregion/rep18/case1/"
SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
SC_PATH <- paste(REP_PATH, "simul_sc.txt", sep="/")
OUTPUT_PATH <- paste(REP_PATH, "bscite", sep="/")

library(dplyr)
library(matrixStats)
library(reshape2)
library(ScRNAClone)
library(TailRank)

bulk <- read.table(SNV_PATH, header=T, as.is = T)
sc <- read.table(SC_PATH, header=T, as.is = T)
snv_count <- dim(bulk)[1]

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}

# Bulk data preparation.
if (!is.numeric(bulk$b) & !is.numeric(bulk$d)) {
    bulk$b <- as.character(bulk$b)
    bulk$d <- as.character(bulk$d)
    V <- matrix(as.numeric(unlist(strsplit(bulk$b, ","))), nrow = snv_count, byrow = T)
    D <- matrix(as.numeric(unlist(strsplit(bulk$d, ","))), nrow = snv_count, byrow = T)
    R <- D - V
    bulk$a <- apply(R, 1, function(row) {
        paste(row, collapse=",")
    })
} else {
    bulk$a <- bulk$d - bulk$b
}
bulk_bscite <- bulk[,c("ID", "CHR", "POS", "b", "a")]
names(bulk_bscite) <- c("ID", "Chromosome", "Position", "MutantCount", "ReferenceCount")
bulk_bscite$INFO <- "NA"
bulk_bscite_file <- paste(OUTPUT_PATH, "/bscite.bulk", sep="")
write.table(bulk_bscite, file=bulk_bscite_file, sep="\t", quote=F, row.names = F, col.names = T)

# Single cell data preparation.
sc$b <- sc$d - sc$a
sc$ID <- factor(sc$ID, levels = bulk$ID)
cells <- unique(sc$Cell)
sc$Cell <- factor(sc$Cell, levels = cells)

sc_var <- reshape2::dcast(sc, ID ~ Cell, value.var = "b")
sc_depth <- reshape2::dcast(sc, ID ~ Cell, value.var = "d")
sc_var[is.na(sc_var)] <- 0
sc_depth[is.na(sc_depth)] <- 0

sc_var <- as.matrix(sc_var[,-1])
sc_depth <- as.matrix(sc_depth[,-1])
n_cells <- length(cells)

# Call variant using threshold.
# Rows: mutations.
# Cols: cells.
sc_mut_matrix <- matrix(0, nrow = length(bulk$ID), ncol = length(cells))
sc_mut_matrix[sc_var >= THRESHOLD] <- 1
sc_mut_matrix[sc_depth == 0] <- 3
colnames(sc_mut_matrix) <- as.character(cells)
rownames(sc_mut_matrix) <- as.character(bulk$ID)
dim(sc_mut_matrix)

sc_bscite_file <- paste(OUTPUT_PATH, "bscite.SC", sep="/")
write.table(sc_mut_matrix, sc_bscite_file, col.names=F, quote=F, row.names = F)
cell_count_file <- paste(OUTPUT_PATH, "cell_count.txt", sep="/")
write.table(n_cells, cell_count_file, quote=F, row.names = F, col.names = F)
