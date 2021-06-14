#rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

print(args)

DATA_PATH <- args[1]
SC_MUT_MATRIX_PATH <- args[2]

#DATA_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/"
#SC_MUT_MATRIX_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/mut_matrix.txt"

SNV_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_READS_PATH <- paste(DATA_PATH, "sc.txt", sep="/")
OUTPUT_PATH <- paste(DATA_PATH, "B-SCITE", sep="/")
if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH)
}

library(dplyr)
library(matrixStats)
library(reshape2)
library(ScRNAClone)
library(TailRank)

bulk <- read.table(SNV_PATH, header=T)
sc_mut_matrix <- read.table(SC_MUT_MATRIX_PATH, header=T)
sc_reads <- read.table(SC_READS_PATH, header=T, as.is = T)

# Bulk data preparation.
bulk$a <- bulk$d - bulk$b
bulk_bscite <- bulk[,c("ID", "CHR", "POS", "b", "a")]
names(bulk_bscite) <- c("ID", "Chromosome", "Position", "MutantCount", "ReferenceCount")
bulk_bscite$INFO <- "NA"
bulk_bscite_file <- paste(OUTPUT_PATH, "/bscite.bulk", sep="")
write.table(bulk_bscite, file=bulk_bscite_file, sep="\t", quote=F, row.names = F, col.names = T)

# Single cell data preparation.
cells <- rownames(sc_mut_matrix)
snvs <- colnames(sc_mut_matrix)
n_cells <- length(cells)
for (i in 1:n_cells) {
    temp <- subset(sc_reads, SampleName == cells[i] & ID %in% snvs)
    sc_mut_matrix[i,temp$d == 0] <- 3
}
sc_bscite_file <- paste(OUTPUT_PATH, "bscite.SC", sep="/")
write.table(t(sc_mut_matrix), sc_bscite_file, col.names=F, quote=F, row.names = F)

n_snvs<-length(snvs)
print(paste("Number of SNVs:", n_snvs))
print(paste("Number of cells:", n_cells))

cell_count_file <- paste(OUTPUT_PATH, "cell_count.txt", sep="/")
write.table(n_cells, cell_count_file, quote=F, row.names = F, col.names = F)

snv_count_file <- paste(OUTPUT_PATH, "snv_count.txt", sep="/")
write.table(n_snvs, snv_count_file, quote=F, row.names = F, col.names = F)
