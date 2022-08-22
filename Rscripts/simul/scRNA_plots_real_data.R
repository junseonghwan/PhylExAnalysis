library(biomaRt)
library(GenomicRanges)
library(PhylExR)

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", GRCh = 37)
grch37_exons <- read.table("~/Dropbox/seong/data/GRCh37_exons.bed", header = T)
head(listAttributes(ensembl)[,1], 20)
bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'external_gene_name', "chromosome_name", "start_position", "end_position"),
            filters = 'ensembl_gene_id',
            values = unique(grch37_exons$ensembl_gene_id),
            mart = ensembl)
bm_ranges <- ConstructGranges(bm$chromosome_name, start = bm$start_position, width = bm$end_position - bm$start_position)

# Let's plot the scRNA-seq data for Laks data.
PlotFunction <- function(p, base_size = 12) {
    p <- p + theme(axis.title = element_text(size = base_size * 2))
    p <- p + theme(axis.text.x = element_text(size = base_size * 2), axis.text.y = element_text(size = base_size * 2))
    p <- p + theme(plot.title = element_text(size = base_size * 2))
    return(p)
}


library(dplyr)
library(ggplot2)
library(matrixStats)
#library(ScRNAClone)
library(PhylExR)
library(TailRank)

# Process simulation data to get sc hyperparameters.
#data_path <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/"
data_path <- "data/HGSOC_SS3/"
ssm <- read.table(paste(data_path, "bulk.txt", sep="/"), header=T, sep="\t")
sc <- read.table(paste(data_path, "sc.txt", sep="/"), header=T)

ssm_ranges <- ConstructGranges(chr = ssm$CHR, start = ssm$POS, width = 0)

# Now, we estimate the hyperparameters for the single cell reads.
n_snvs <- dim(ssm)[1]
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
p <- ggplot(temp, aes(b/d)) + geom_histogram(fill="black") + theme_bw()
p <- p + ylab("Count") + xlab("Var reads / Depth") + ggtitle("HGSOC scRNA-seq histogram")
p <- PlotFunction(p)
ggsave(filename = "_figures/supp/Laks_scRNA.pdf", p, height = 8, width = 8, units = "in")

temp2 <- temp %>% subset(b > 0 & b/d < 1)
p <- ggplot(temp2, aes(b/d)) + geom_histogram(fill="black") + theme_bw()
p <- p + ylab("Count") + xlab("Var reads / Depth") + ggtitle("HGSOC scRNA-seq histogram")
p <- PlotFunction(p)
ggsave(filename = "_figures/supp/Laks_scRNA_biallelic.pdf", p, height = 8, width = 8, units = "in")

# Supplementary Figure 7d.
col_names <- c("ID", "Cell", "a", "b", "d")
write.csv(temp[,col_names], file = "data/NatComm/SupplementaryFigure7d.csv", row.names = F, quote = F)
# Supplementary Figure 7f.
write.csv(temp2[,col_names], file = "data/NatComm/SupplementaryFigure7f.csv", row.names = F, quote = F)

# Make a figure for HGSOC using gene names and co-occurrence for each cell.
overlaps <- findOverlaps(ssm_ranges, bm_ranges)
idx <- !duplicated(overlaps@from)
ssm_ <- ssm[overlaps@from[idx],]
ssm_$Gene <- bm[overlaps@to[idx],"external_gene_name"]

# x-axis: genes
# y-axis: cells
sc_ <- dplyr::left_join(sc, ssm_, by ="ID")
head(sc_)
sc_ <- sc_[,c("Cell", "b.x", "d.x", "Gene")]
sc_.melted <- reshape2::melt(sc_)
sc_.melted <- subset(sc_.melted, "b.x" > 0)
dim(sc_.melted)
ggplot(sc_.melted, aes(Gene, Cell, b.x)) + geom_tile()
sc_.melted

# Make a similar figure for HER2+.
#data_path <- "/Users/seonghwanjun/data/TNBC/BCSA2/MultiRegion/"
data_path <- "data/HER2_POS_SS3/"
ssm <- read.table(paste(data_path, "bulk.txt", sep="/"), header=T, sep="\t")
sc <- read.table(paste(data_path, "sc.txt", sep="/"), header=T)

n_snvs <- dim(ssm)[1]
sc$b <- sc$d - sc$a
temp <- subset(sc, d > 0)
p <- ggplot(temp, aes(b/d)) + geom_histogram(fill="black") + theme_bw()
p <- p + ylab("Count") + xlab("Var reads / Depth") + ggtitle("HER2+ scRNA-seq histogram")
p <- PlotFunction(p)
ggsave(filename = "_figures/supp/HER2_scRNA.pdf", p, height = 8, width = 8, units = "in")

temp2 <- temp %>% subset(b > 0 & b/d < 1)
p <- ggplot(temp2, aes(b/d)) + geom_histogram(fill="black") + theme_bw()
p <- p + ylab("Count") + xlab("Var reads / Depth") + ggtitle("HER2+ scRNA-seq histogram")
p <- PlotFunction(p)
ggsave(filename = "_figures/supp/HER2_scRNA_biallelic.pdf", p, height = 8, width = 8, units = "in")

# Supplementary Figure 7e.
write.csv(temp[,col_names], file = "data/NatComm/SupplementaryFigure7e.csv", row.names = F, quote = F)
# Supplementary Figure 7g.
write.csv(temp2[,col_names], file = "data/NatComm/SupplementaryFigure7g.csv", row.names = F, quote = F)
