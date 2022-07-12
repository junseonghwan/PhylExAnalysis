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

sc <- read.table("data/HER2_POS_SS3/sc.txt", header = T)
sc_hp <- read.table("data/HER2_POS_SS3/sc_hp.txt", header = T)
bursty_hp <- list(alpha=0.01, beta=0.01)


rep_no <- FindBestRep(paste(DATA_PATH, "tssb", sep="/"), chains = 0:3)
rep_path <- paste(DATA_PATH, "/tssb/chain", rep_no, sep="")
datum2node_tssb <- read.table(paste(rep_path, "joint/tree0/datum2node.tsv", sep="/"), header=F, sep="\t", as.is = T)
names(datum2node_tssb) <- c("ID", "Node")
table(datum2node_tssb$Node)

#sc <- read.table(paste(DATA_PATH, "sc.txt", sep="/"), header=T, sep="\t")
datum2node_tssb <- read.table(paste(rep_path, "joint/tree0/datum2node.tsv", sep="/"), header=F, sep="\t", as.is = T)
names(datum2node_tssb) <- c("ID", "Node")
table(datum2node_tssb$Node)

cell.df_tssb <- AssignCellsBursty(sc, datum2node_tssb, bursty_hp = bursty_hp, biallelic_hp = sc_hp, include_normal_clone = F)
table(cell.df_tssb$Node)

# Let's take a deeper look at clone 0_0_0_0_0_1 from TSSB.
branched_node <- "0_0_0_0_0_1"
branched_clone <- datum2node_tssb[datum2node_tssb$Node == branched_node,]
# Let's check single cell data.
branched_cells <- subset(cell.df_tssb, Node %in% branched_node)
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
sc_ <- left_join(sc_, cell.df_tssb)
sc_ <- subset(sc_, Node %in% c(branched_node, sibling_node))
sc_$b <- sc_$d - sc_$a
sc_ <- subset(sc_, b > 0)
ids <- unique(sc_$ID)
cells_clustered <- cell.df_tssb[order(cell.df_tssb$Node), "Cell"]
cells_clustered <- cells_clustered[cells_clustered %in% unique(sc_$Cell)]
sc_$Cell <- factor(sc_$Cell, levels = cells_clustered)
mut_ids_clustered <- datum2node_tssb[order(datum2node_tssb$Node), "ID"]
mut_ids_clustered <- mut_ids_clustered[mut_ids_clustered %in% ids]
sc_$ID <- factor(sc_$ID, levels = mut_ids_clustered)
base_size <- 11

sc_$Clone <- sc_$Node
sc_[sc_$Node %in% sibling_node,"Clone"] <- "Clone5"
sc_[sc_$Node %in% branched_clone$Node,"Clone"] <- "Clone6"
p <- ggplot(subset(sc_, b > 0), aes(ID, Cell, fill = Clone)) + geom_tile(color = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cell")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x = element_text(size = base_size * 2))
p <- p + theme(axis.title.y = element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.text = element_text(size = base_size * 1.5, face = "bold"))
p <- p + guides(fill = guide_legend(title = "TSSB (SNV) Clone", size = base_size))
p <- p + annotate("rect", xmin = which(mut_ids_clustered == "s343") - 0.5, xmax = which(mut_ids_clustered == "s352") + 0.5, ymin = 0, ymax = length(cells_clustered) + 0.5, fill = "red", alpha = 0.2)
#p <- p + annotate("rect", xmin = which(mut_ids_clustered == "s15") - 0.5, xmax = which(mut_ids_clustered == "s343") - 0.5, ymin = 0, ymax = length(cells_clustered) + 0.5, fill = "gray", alpha = 0.2)
p
ggsave(p, filename = "~/phylo-express-paper/figures/HER2_pos/TSSB_branching_sc_plot.pdf", width = 6, height = 11, units = "in")
ggsave(p, filename = "_figures/HER2_POS/TSSB_branching_sc_plot.pdf", width = 6, height = 11, units = "in")


# PhylEx
rep_no <- FindBestRep(paste(DATA_PATH, "phylex", sep="/"), chains = 0:3)
rep_path <- paste(DATA_PATH, "/phylex/chain", rep_no, sep="")

best_chain <- FindBestRep("data/HER2_POS_SS3/phylex/", chains = 0:3)
datum2node <- read.table(paste("data/HER2_POS_SS3/phylex/chain", best_chain, "/joint/tree0/datum2node.tsv", sep=""), header = F)
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)

datum2node <- CollapseClones(datum2node, MIN_SNV_COUNT = 1)
table(datum2node$Node)
table(cell.df$Node)
cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp = bursty_hp, biallelic_hp = sc_hp, include_normal_clone = F)
table(cell.df$Node)

ret <- CollapseClonesByCellCount(datum2node, cell.df)
datum2node <- ret$datum2node
cell.df <- ret$cell.df

datum2node$Node <- as.character(datum2node$Node)
datum2node.sorted <- datum2node[order(datum2node$Node),]
datum2node$Clone <- datum2node$Node

# Re-label the clone names.
nodes <- unique(datum2node$Node)
nodes_ordered <- nodes[order(nchar(nodes), nodes)]
clone_count <- length(nodes_ordered)
for (i in 1:clone_count) {
  idx <- datum2node$Node == nodes_ordered[i]
  datum2node$Clone[idx] <- i
  idx <- cell.df$Node == nodes_ordered[i]
  cell.df$Clone[idx] <- paste("Clone", i, sep="")
}

cell_order <- cell.df[order(cell.df$Node), "Cell"]
sc_join <- left_join(sc, cell.df, by = "Cell")
sc_join <- subset(sc_join, ID %in% datum2node.sorted$ID)
sc_join$Cell <- factor(sc_join$Cell, levels = cell_order)
sc_join$ID <- factor(sc_join$ID, levels = datum2node.sorted$ID)
sc_join$b <- sc_join$d - sc_join$a

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
p <- p + theme(legend.text = element_text(size=base_size)) + guides(fill=guide_legend(title="PhylEx (SNV) Clone", size = base_size))
p
ggsave(p, filename = "_figures/HER2_POS/Coclustering.pdf", width = 6, height = 11.5, units = "in")

# InferCNV

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

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=exp.rawdata,
                                     ref_group_names=NULL,
                                     annotations_file = ANNOTATION_FILE,
                                     gene_order_file=GENE_ORDER_FILE) 

infercnv_obj_hmm_i6 <- infercnv::run(infercnv_obj,
                                     cutoff=1, 
                                     out_dir="data/HER2_POS_SS3/infercnv/hmm_i6_subclusters/", 
                                     cluster_by_groups=FALSE, 
                                     denoise=TRUE,
                                     HMM=TRUE,
                                     HMM_type = "i6",
                                     analysis_mode = "subclusters")


infercnv_cnv_pred_genes <- read.table("data/HER2_POS_SS3/infercnv/hmm_i6_subclusters/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat", header = T)
infercnv_cnv_pred_regions <- read.table("data/HER2_POS_SS3/infercnv/hmm_i6_subclusters/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat", header = T)
infercnv_cell_groupings <- read.table("data/HER2_POS_SS3/infercnv/hmm_i6_subclusters/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings", header = T)

# Generate UMAP.
library(Seurat)
dim(exp.rawdata)
seurat <- CreateSeuratObject(counts = exp.rawdata)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)
DimPlot(seurat)

mean(cell.df$SampleName == colnames(seurat))
seurat_ <- seurat[,colnames(seurat) %in% paste(unique(sc_join$SampleName), "Aligned", sep="")]
Idents(seurat_) <- cell.df$Node
DimPlot(seurat_)

seurat_infercnv <- seurat[,colnames(seurat) %in% infercnv_cell_groupings$cell]
midx <- match(colnames(seurat_infercnv), infercnv_cell_groupings$cell)
mean(infercnv_cell_groupings$cell[midx] == colnames(seurat_infercnv))
Idents(seurat_infercnv) <- infercnv_cell_groupings$cell_group_name[midx]
DimPlot(seurat_infercnv)

# Generate heatmap.
# There should be clear separation in the SNVs.
sc_join$SampleName <- paste(sc_join$SampleName, "Aligned", sep="")
infercnv_cell_groupings_ <- infercnv_cell_groupings[infercnv_cell_groupings$cell %in% unique(sc_join$SampleName),]
table(infercnv_cell_groupings_$CNVClone)
names(infercnv_cell_groupings_) <- c("CNVClone", "SampleName")
sc_join_ <- left_join(sc_join, infercnv_cell_groupings_)
names(sc_join_)

sc_join_ <- sc_join_[order(sc_join_$CNVClone),]
cell_order <- unique(sc_join_$Cell)
sc_join_$Cell <- factor(sc_join_$Cell, levels = cell_order)
sc_join_$CNVClone <- gsub(pattern = "all_observations.all_observations.", replacement = "", sc_join_$CNVClone)
sc_join_$CNVClone <- gsub(pattern = "\\.", replacement = "_", sc_join_$CNVClone)

p <- ggplot(subset(sc_join_, b >0), aes(ID, Cell, fill=CNVClone)) + geom_tile(colour = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cell")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.text = element_text(size=base_size)) + guides(fill=guide_legend(title="CNV Clone", size = base_size))
p
ggsave(p, filename = "_figures/HER2_POS/Coclustering_InferCNV.pdf", width = 6, height = 11.5, units = "in")
