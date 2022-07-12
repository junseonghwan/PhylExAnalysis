rm(list=ls())
library(biomaRt)
library(copykat)
library(dplyr)
library(infercnv)
library(ggplot2)
library(PhylExR)
library(Seurat)

setwd("~/PhylExAnalysis/")

# TSSB results:
reps <- 0:3
rep_count <- length(reps)
DATA_PATH <- "data/HER2_POS_SS3/"
rep_no <- FindBestRep(paste(DATA_PATH, "tssb", sep="/"), chains = 0:3)
rep_path <- paste(DATA_PATH, "/tssb/chain", rep_no, sep="")

#sc <- read.table(paste(DATA_PATH, "sc.txt", sep="/"), header=T, sep="\t")
datum2node_tssb <- read.table(paste(rep_path, "joint/tree0/datum2node.tsv", sep="/"), header=F, sep="\t", as.is = T)
names(datum2node_tssb) <- c("ID", "Node")
table(datum2node_tssb$Node)

# PhylEx results:
sc <- read.table("data/HER2_POS_SS3/sc.txt", header = T)
sc_hp <- read.table("data/HER2_POS_SS3/sc_hp.txt", header = T)
best_chain <- FindBestRep("data/HER2_POS_SS3/phylex/", chains = 0:3)
datum2node <- read.table(paste("data/HER2_POS_SS3/phylex/chain", best_chain, "/joint/tree0/datum2node.tsv", sep=""), header = F)
names(datum2node) <- c("ID", "Node")
bursty_hp <- list(alpha=0.01, beta=0.01)
cells.df <- AssignCellsBursty(sc, datum2node, bursty_hp = bursty_hp, biallelic_hp = sc_hp, include_normal_clone = F)
cells.df_w_normal <- AssignCellsBursty(sc, datum2node, bursty_hp = bursty_hp, biallelic_hp = sc_hp, include_normal_clone = T)
cells.df_tssb <- AssignCellsBursty(sc, datum2node_tssb, bursty_hp = bursty_hp, biallelic_hp = sc_hp, include_normal_clone = F)
table(datum2node$Node)
table(cells.df$Node)
table(cells.df_w_normal$Node)
table(cells.df_tssb$Node)

# Let's take a deeper look at clone 0_0_0_0_0_1 from TSSB.
branched_node <- "0_0_0_0_0_1"
branched_clone <- datum2node_tssb[datum2node_tssb$Node == branched_node,]
# Where are they in datum2node?
subset(datum2node, ID %in% branched_clone$ID)
# Let's check single cell data.
branched_cells <- subset(cells.df_tssb, Node %in% branched_node)
temp_sc <- subset(sc, ID %in% branched_clone$ID & Cell %in% branched_cells$Cell)
ret <- temp_sc %>% group_by(Cell) %>% summarise(n = sum(d - a > 0))
# All the cells have at least one of these mutations.
ret$n
# Cell count for each mutation.
ret <- temp_sc %>% group_by(ID) %>% summarise(n = sum(d - a > 0))
ret
# Do these branched cells have any other mutation?
temp_sc <- subset(sc, Cell %in% branched_cells$Cell)
ret <- temp_sc %>% group_by(Cell) %>% summarise(n = sum(d - a > 0))
# Yes.
ret$n
# Are any of the mutations from its sibling?
muts <- temp_sc %>% group_by(ID) %>% summarise(n = sum(d - a > 0))
sibling_node <- "0_0_0_0_0_0"
# Yes, 24 of them are actually in the sibling node.
sum(subset(datum2node_tssb, ID %in% muts$ID)$Node == sibling_node)
# Generate a heatmap figure focusing on the mutations found in branched_node and sibling_node.
ids <- datum2node_tssb[datum2node_tssb$Node %in% c(branched_node, sibling_node),]
sc_ <- subset(sc, ID %in% ids$ID)
sc_ <- left_join(sc_, cells.df_tssb)
sc_ <- subset(sc_, Node %in% c(branched_node, sibling_node))
sc_$b <- sc_$d - sc_$a
sc_ <- subset(sc_, b > 0)
ids <- unique(sc_$ID)
cells_clustered <- cells.df_tssb[order(cells.df_tssb$Node), "Cell"]
cells_clustered <- cells_clustered[cells_clustered %in% unique(sc_$Cell)]
sc_$Cell <- factor(sc_$Cell, levels = cells_clustered)
mut_ids_clustered <- datum2node_tssb[order(datum2node_tssb$Node), "ID"]
mut_ids_clustered <- mut_ids_clustered[mut_ids_clustered %in% ids]
sc_$ID <- factor(sc_$ID, levels = mut_ids_clustered)
base_size <- 11
p <- ggplot(subset(sc_, b > 0), aes(ID, Cell, fill = b)) + geom_tile(color = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cell")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x = element_text(size = base_size * 2))
p <- p + theme(axis.title.y = element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.text = element_text(size = 7, face = "bold"))
p <- p + guides(fill = guide_legend(title = "Variant reads", size = 7))
p <- p + annotate("rect", xmin = which(mut_ids_clustered == "s343") - 0.5, xmax = which(mut_ids_clustered == "s352") + 0.5, ymin = 0, ymax = length(cells_clustered) + 0.5, fill = "red", alpha = 0.2)
#p <- p + annotate("rect", xmin = which(mut_ids_clustered == "s15") - 0.5, xmax = which(mut_ids_clustered == "s343") - 0.5, ymin = 0, ymax = length(cells_clustered) + 0.5, fill = "gray", alpha = 0.2)
ggsave(p, filename = "~/phylo-express-paper/figures/revision/HER2_pos/TSSB_branching_sc_plot.pdf", width = 6, height = 11, units = "in")
ggsave(p, filename = "_figures/HER2_POS/TSSB_branching_sc_plot.pdf", width = 6, height = 11, units = "in")

# Load SS3 data.
fc <- read.table("data/HER2_POS_SS3/featureCounts.txt", header = T)
MIN_CELLS <- 5
seurat <- CreateSeuratObject(fc, min.cells = MIN_CELLS)
exp.rawdata <- as.matrix(seurat@assays$RNA@counts)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convert_to_ensembl_id <- function(gene_names, mart) {
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                   filters = "external_gene_name", values = gene_names,
                   mart = mart)
  return(results)
}
convert_to_gene_name <- function(ensembl_gene_ids, mart) {
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = "ensembl_gene_id", values = ensembl_gene_ids,
                   mart = mart)
  return(results)
}

# Prepare input files for InferCNV.
chrs <- c(1:22, "X", "Y")

ANNOTATION_FILE <- "data/HER2_POS_SS3/infercnv/annotation.tsv"
GENE_ORDER_FILE <- "data/HER2_POS_SS3/infercnv/gene_order.tsv"

ret <- convert_to_gene_name(rownames(exp.rawdata), mart)
ret <- subset(ret, chromosome_name %in% chrs)
ret$chromosome_name <- factor(ret$chromosome_name, levels = chrs)
ret_ <- ret[order(ret$chromosome_name, ret$start_position),]
unique(ret_$chromosome_name)
ret_$chromosome_name <- paste("chr", ret_$chromosome_name, sep="")
write.table(ret_[,c(1, 3, 4, 5)], file = GENE_ORDER_FILE, sep="\t", quote = F, row.names = F, col.names = F)
write.table(data.frame(Cell=colnames(exp.rawdata), Annotation="Cancer"), file = ANNOTATION_FILE, quote = F, col.names = F, row.names = F, sep="\t")

# Run InferCNV.
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=exp.rawdata,
                                     ref_group_names=NULL,
                                     annotations_file = ANNOTATION_FILE,
                                     gene_order_file=GENE_ORDER_FILE) 

infercnv_obj_hmm_i6 <- infercnv::run(infercnv_obj,
                                     cutoff=1, 
                                     out_dir="data/HER2_POS_SS3/infercnv/", 
                                     cluster_by_groups=TRUE, 
                                     denoise=TRUE,
                                     HMM=TRUE,
                                     HMM_type = "i6",)
saveRDS(infercnv_obj_hmm_i6, file = "data/HER2_POS_SS3/infercnv/infercnv_obj_hmm_i6.RDS")

plot(infercnv_obj_hmm_i6@tumor_subclusters$hc$Cancer, labels = FALSE)
infercnv_clones <- cutree(infercnv_obj_hmm_i6@tumor_subclusters$hc$Cancer, 5)

# Run CopyKAT and analyze the results.
dim(exp.rawdata)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="E", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4, cell.line = F)
pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)
saveRDS(copykat.test, file = "data/HER2_POS_SS3/copykat//copykat_object.RDS")

plot(copykat.test$hclustering, label = FALSE)
copykat_clones <- cutree(copykat.test$hclustering, k = 5)

# Let's check concordances between the two methods.
infercnv_clones_ <- infercnv_clones[names(infercnv_clones) %in% names(copykat_clones)]
length(infercnv_clones_) == length(copykat_clones)
mean(names(infercnv_clones_) == names(copykat_clones))
# Low concordance between the two as measured by adjusted Rand Index.
mclust::adjustedRandIndex(copykat_clones, infercnv_clones_)

# Compare to SNV clones.
sc_ <- left_join(sc, cells.df)
names(sc_)
snv_clones <- unique(sc_[,c("Node", "SampleName")])
snv_clones$SampleName <- paste(snv_clones$SampleName, "Aligned", sep="")
snv_clones <- snv_clones[snv_clones$SampleName %in% names(infercnv_clones_),]
snv_clones <- snv_clones[match(names(copykat_clones), snv_clones$SampleName),]
mclust::adjustedRandIndex(copykat_clones, snv_clones$Node)
mclust::adjustedRandIndex(infercnv_clones_, snv_clones$Node)

seurat <- seurat %>% NormalizeData() %>% ScaleData
seurat <- FindVariableFeatures(seurat)
seurat <- RunPCA(seurat)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.8)
seurat <- RunUMAP(seurat, dims = 1:30)
DimPlot(seurat, label = TRUE)
cell_names_to_analyze <- unique(paste(sc$SampleName, "Aligned", sep=""))
seurat$cell_names <- colnames(seurat)
seurat_ <- subset(seurat, subset = cell_names %in% cell_names_to_analyze)
ordering <- match(seurat_$cell_names, cell_names_to_analyze)
seurat_$cell_names == cell_names_to_analyze[ordering]
Idents(seurat_) <- cells.df$Node[ordering]
DimPlot(seurat_)

# Compute coverage for the genes used in Figure 4b.
genes <- convert_to_ensembl_id(gene_names = c("CDC6", "FN1", "WNT10A", "PITX2", "EZH2", "PRKACG", 
                                     "CACNG4", "NF1", "POLE2", "DKK2", "ETS2", "FGF14", 
                                     "ACVR1B", "PRKDC", "COL4A5", "VPS33B", "IL2RA", "MDC1", 
                                     "CACNA2D2", "PIK3R3", "TP53", "FOS", "DDX50", "MAP3K8"), mart)

genes <- genes[genes$chromosome_name %in% chrs,]
fc_ <- fc[rownames(fc) %in% genes$ensembl_gene_id,]
cell_count <- rowSums(fc_ > 0)
mean_reads <- rowMeans(fc_)
genes_ <- genes[match(rownames(fc_), genes$ensembl_gene_id),]
genes_$cell_count <- 0
genes_$mean_reads <- 0
dat.df <- data.frame(cell_count = cell_count, mean_reads = mean_reads, genes = genes_$external_gene_name)
dat.df_ <- dat.df[dat.df$cell_count > 0,]

idx <- match(dat.df_$genes, genes_$external_gene_name)
genes_[idx,]$cell_count <- temp$cell_count
genes_[idx,]$mean_reads <- temp$mean_reads

# Output as csv file.
genes_ <- genes_[,-c(1, 3)]
names(genes_) <- c("Gene", "Cell Count", "Mean Reads")
xlsx::write.xlsx(genes_, file = "~/Dropbox/seong/PhylExNatComm/SupplementTable5.xlsx", row.names = F)
xlsx::write.xlsx(genes_, file = "~/phylo-express-paper/tables/SupplementTable5.xlsx", row.names = F)


