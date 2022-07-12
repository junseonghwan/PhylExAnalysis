# Replace HGSOC_SS3_DGE with this script.
rm(list=ls())
library(apeglm)
library(biomaRt)
library(BiocParallel)
library(copykat)
library(DESeq2)
library(edgeR)
library(gam)
library(ggplot2)
library(Glimma)
library(infercnv)
library(limma)
library(matrixStats)
library(magrittr)
library(mclust)
library(pheatmap)
library(Rcpp)
library(PhylExR)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(TSCAN)
library(tsne)
library(xtable)
library(zinbwave)

DATA_PATH <- "~/PhylExAnalysis/data/HGSOC_SS3/"
FEATURE_COUNTS_FILE <- "~/PhylExAnalysis/data/HGSOC_SS3/featureCounts.txt"
BULK_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "sc.txt", sep="/")
SC_HP_PATH <- paste(DATA_PATH, "sc_hp.txt", sep="/")
GT_PATH <- paste(DATA_PATH, "gt.txt", sep="/")
ANALYSIS_OUTPUT_PATH <- "~/PhylExAnalysis/_figures/HGSOC/"

topK <- 1000
base_size <- 12
MIN_CELLS <- 10
MIN_READS <- 5
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chrs <- c(1:22, "X", "Y")
#mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = "37")

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

plot_sig_genes <- function(dge_results, mart, gene_plot_count = 50) {
  # Plot the significantly differentially expressed genes along with their p-values.
  ordering <- order(dge_results$table$PValue)
  tbl <- dge_results$table[ordering,]
  tbl$ensembl_gene_id <- rownames(tbl)
  
  bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
                           'start_position', 'end_position', 'percentage_gene_gc_content'),
              filters = 'ensembl_gene_id',
              values = tbl$ensembl_gene_id,
              mart = mart, uniqueRows = T)
  
  ret <- dplyr::left_join(tbl, bm, by = "ensembl_gene_id")
  ret <- ret[!is.na(ret$hgnc_symbol),]
  ret.df <- data.frame(gene = ret$hgnc_symbol[1:gene_plot_count], logpvalue = -log(ret$PValue)[1:gene_plot_count], logFC=ret$logFC[1:gene_plot_count] > 0)
  ret.df$gene <- factor(ret.df$gene, levels = ret.df$gene)
  ret.df$logFC[ret.df$logFC == FALSE] <- "DOWN"
  ret.df$logFC[ret.df$logFC == TRUE] <- "UP"
  p <- ggplot(ret.df, aes(gene, logpvalue, col=logFC)) + geom_point() + theme_bw() + xlab("Gene") + ylab("-Log p-value")
  p <- p + scale_color_manual(values=c("DOWN"="blue", "UP"="red"))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title = element_blank())
  return(p)
}

# Get cell assignment.
# Find the best rep -- we will generate figures based on the best rep first.
reps <- 0:19
log_liks <- rep(0, length(reps))
for (rep_no in reps) {
  log_liks[rep_no + 1] <- read.table(paste(DATA_PATH, "/phylex/chain", rep_no, "/joint/tree0/log_lik.txt", sep=""), sep="", header=F)$V1
}
best_rep <- which.max(log_liks) - 1

bursty_hp <- list(alpha=0.01, beta=0.01)
dat <- read.table(paste(DATA_PATH, "bulk.txt", sep="/"), header=T, sep="\t")
#write.table(dat[,c(1,2,3)], paste(DATA_PATH, "/id2loc.tsv", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
sc <- read.table(SC_PATH, header=T, as.is = TRUE)
sc_hp <- read.table(paste(DATA_PATH, "sc_hp.txt", sep="/"), header=T, sep="\t")
datum2node <- read.table(paste(DATA_PATH, "/phylex/chain", best_rep, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t")
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)

# Get the tree
valid_clones <- c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "A_B", "C_D", "A", "B", "C", "D", "E_F_G_H_I", "E_F", "E", "F")
gt <- read.table(GT_PATH, header=T)
validation_idx <- which(gt$CloneName %in% valid_clones)
table(gt$CloneName[validation_idx])
unique(datum2node[validation_idx,"Node"])
temp <- cbind(datum2node[validation_idx,], gt[validation_idx,"CloneName"])
names(temp) <- c("ID", "Node", "CloneName")
tree <- temp[order(temp$Node),]
tree
dim(tree)

# Assign cell to nodes.
cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = FALSE)
#write.table(cell.df, paste(DATA_PATH, "/results/rep", best_rep, "/joint/tree0/cell2node.tsv", sep=""), quote = F, col.names = T, row.names = F, sep="\t")
cell.df$SampleName <- unique(sc$SampleName)
cell.df$Node <- as.character(cell.df$Node)
cell.df$Clone <- cell.df$Node
table(cell.df$Clone)
cell.df$Clone[cell.df$Clone == "0_0_0_0_0"] <- "C"
cell.df$Clone[cell.df$Clone == "0_0_0_0"] <- "CD"
# Although we have reason to believe that our method can detect AB clone, there are only 3 cells assigned to AB. 
# Also, the SNV that belongs in 0_0_0_1, s54 is mislabelled in the ground truth as A_B_I => applying the same criteria 
# to eliminate SNVs that are incorrectly labelled, we decide to omit this clone.
#cell.df$Clone[cell.df$Clone == "0_0_0_1"] <- "AB" 
cell.df$Clone[cell.df$Clone == "0_0_0" | cell.df$Clone == "0_0_0_1"] <- "ABCD"
cell.df$Clone[cell.df$Clone == "0_0_1"] <- "EF"
cell.df$Clone[cell.df$Clone == "0_0"] <- "Ancestral"
#valid_clones <- c("Ancestral", "ABCD", "AB", "CD", "C", "EF")
valid_clones <- c("Ancestral", "ABCD", "CD", "C", "EF")
cell.df$Clone <- factor(cell.df$Clone, levels = valid_clones)
table(cell.df$Clone)

cell_order <- cell.df[order(cell.df$Node), "Cell"]
nodes <- unique(cell.df$Node)
nodes <- nodes[order(nchar(nodes), nodes)]

# Retrieve the output of the feature counts.
fc <- read.table(FEATURE_COUNTS_FILE, header=T, row.names = 1, check.names = F)
feature_counts <- fc[,names(fc) %in% unique(sc$SampleName)]
dim(feature_counts)

# We will filter out some genes, e.g., remove MT genes, genes with unknown chromosome.
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', "chromosome_name"),
            filters = 'ensembl_gene_id',
            values = rownames(feature_counts),
            mart = mart, uniqueRows = T)
idx <- match(rownames(feature_counts), bm$ensembl_gene_id)
bm <- bm[idx,]
idx <- which(bm$chromosome_name %in% chrs)
feature_counts <- feature_counts[idx,]
bm <- bm[idx,]
mean(bm[,1] == rownames(feature_counts), na.rm = T)
mean(names(feature_counts) == cell.df$SampleName)
dim(feature_counts)

# First, perform t-SNE analysis.
sce <- SingleCellExperiment(assays = list(counts = as.matrix(feature_counts)),
                            rowData = data.frame(EnsemblID = rownames(feature_counts)),
                            colData = data.frame(CellName = cell.df$SampleName, Node=cell.df$Node, Clone=cell.df$Clone))
#logcounts(sce) <- log2(counts(sce) + 1)

filter <- rowSums(assay(sce) > MIN_READS) > MIN_CELLS
table(filter)
sce <- sce[filter,]

# Let's select top K genes.
vars <- assay(sce) %>% log1p %>% rowVars
length(vars)
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
head(vars)
sce <- sce[names(vars)[1:topK],]
sce

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
clone_colors <- c("Ancestral"="#999999", "ABCD"="#E69F00", "CD"="#56B4E9", "EF"="#009E73", "C"="#CC79A7")
cnv_clone_colors <- c("1_1_1_1"="#E69F00", "1_1_1_2"="#56B4E9", "1_1_2_1"="#009E73", "1_1_2_2"="#CC79A7",
                      "1_2_1_1"="#F0E442", "1_2_1_2"="#0072B2", "1_2_2"="#D55E00")
canopy_clone_colors <- c("clone2"="#999999", "clone3"="#E69F00", "clone4"="#56B4E9", 
                         "clone5"="#009E73", "clone6"="#CC79A7", "clone7"="#F0E442", 
                         "clone8"="#0072B2", "clone9"="#D55E00", "clone10"="#994F00")


# df: Data frame with columns W1, W2, Node, and clone.
ReducedDimensionPlots <- function(df, plot_title, base_size = 12, custom_cols = clone_colors) {
  #p <- ggplot(df, aes(W1, W2, colour=Clone, shape=Clone)) + geom_point(alpha=0.8) + theme_classic()
  p <- ggplot(df, aes(W1, W2, colour=Clone)) + geom_point(alpha=0.8) + theme_classic()
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
  p <- p + theme(legend.text = element_text(size=base_size*2), legend.title = element_text(size=base_size*2))
  p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  p <- p + xlab("Dim1") + ylab("Dim2")
  #p <- p + ggtitle(plot_title)
  p <- p + theme(title = element_text(size = base_size * 2))
  p <- p + scale_color_manual(values = custom_cols)
  return(p)
}

# 2D ZINB-WaVE.
zinb <- zinbwave(sce, K = 2, epsilon=1000)
W <- reducedDim(zinb)
df <- data.frame(W1=W[,1], W2=W[,2], Node=colData(sce)$Node, Clone=colData(sce)$Clone)
pl <- ReducedDimensionPlots(df, "ZINB-WaVE dimensions")
pl

# Plot just ABCD, CD, and EF clones.
cls <- c("ABCD", "CD", "EF")
cell_idx <- colData(sce)$Clone %in% cls
clones <- colData(sce)$Clone[cell_idx]
clones <- factor(clones, levels = cls)
df <- data.frame(W1=W[cell_idx,1], W2=W[cell_idx,2], Node=colData(sce)$Node[cell_idx], Clone=clones)
pl <- ReducedDimensionPlots(df, "ZINB-WaVE dimensions", custom_cols = clone_colors[cls])
pl2 <- pl + geom_point(aes(W1, W2, shape = Clone), size = 6, alpha = 0.8)
pl2 <- pl2 + scale_shape_manual(values = c(15, 17, 16))
pl2
ggsave(pl2, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_major_clones.pdf", sep="/"))

# 2D UMAP
logcounts(sce) <- log1p(counts(sce))
umap <- runUMAP(sce)
umapW <- reducedDim(umap)
df_umap <- data.frame(W1=umapW[,1], W2=umapW[,2], Node=colData(sce)$Node, Clone=colData(sce)$Clone)
pl <- ReducedDimensionPlots(df_umap, "UMAP dimensions")
pl

# Make a plot of each clone separately.
clones_to_plot <- c("Ancestral", "ABCD", "CD", "EF", "C")
pl <- ReducedDimensionPlots(df, paste("ZINB-WaVE dimensions. Clone", clone))
for (clone in clones_to_plot)
{
  pl2 <- pl + geom_point(data = subset(df, Clone == clone), size = 5)
  #pl2 <- pl2 + geom_point(data = subset(df, Clone != clone), size = 0.5)
  ggsave(pl2, filename = paste(ANALYSIS_OUTPUT_PATH, "/ZINBWaVE_Top1000Genes_Clone", clone, ".pdf", sep=""), width=8, height=8, units = "in")
}

# Enlarge the point sizes for trajectory analysis.
sce_zinb <- zinbwave(sce, K = 2, epsilon=1000, normalizedValues=TRUE, observationalWeights = TRUE)
# We will replace clone C by CD to reproduce Figure 4.
sce_zinb$Clone[sce_zinb$Clone == "C"] <- "CD"
rd <- reducedDims(sce_zinb)$zinbwave
cl <- Mclust(rd)$classification
colData(sce_zinb)$Cluster <- cl

sim <- slingshot(sce_zinb, clusterLabels = 'Cluster', reducedDim = 'zinbwave', reassign = T)
lineage <- SlingshotDataSet(sim)

W <- reducedDims(sim)$zinbwave
cluster_count <- length(unique(cl))
center <- matrix(0, ncol = 2, nrow = cluster_count)
for (i in 1:cluster_count) {
  center[i,] <- colMeans(W[sim$Cluster == i,])
}
rownames(center) <- paste("", 1:cluster_count, sep="")
center.df <- data.frame(W1=center[,1], W2=center[,2], Cluster=rownames(center), Clone=NA)

cluster_colors <- c("1"="#C3D7A4", "2"="#CC79A7", "3"="#0072B2", "4"="#D55E00")
df <- data.frame(W1 = W[,1], W2 = W[,2], Cluster=as.character(sim$Cluster), Clone=sim$Clone)
#p <- ggplot(df, aes(W1, W2, color = Cluster, shape = Clone)) + geom_point(size = 2.5) + theme_classic()
p <- ggplot(df, aes(W1, W2, color = Cluster)) + geom_point(size = 2.5) + theme_classic()
p <- p + theme(legend.text = element_text(size=base_size*2), legend.title = element_text(size=base_size*2))
p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
p <- p + scale_color_manual(values = cluster_colors)
p <- p + xlab("Dim1") + ylab("Dim2")
#p <- p + ggtitle("ZINB-WaVE using Mclust with pseudo-time analysis")
p <- p + theme(title = element_text(size = base_size * 2))
p <- p + annotate("point", x = center[,1], y = center[,2], size = 2)
lineages <- lineage@lineages
segment <- data.frame()
for (line in lineages) {
  line_num <- as.numeric(line)
  for (i in 1:(length(line_num)-1)) {
    idx1 <- line_num[i]
    idx2 <- line_num[i+1]
    segment_row <- data.frame(x = center.df[idx1,1], y = center.df[idx1,2], xend = center.df[idx2,1], yend = center.df[idx2,2])
    segment <- rbind(segment, segment_row)
  }
}
p <- p + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data = segment, arrow =  arrow(length = unit(0.05, "npc")), arrow.fill='blue', col='black', inherit.aes = FALSE)
p
ggsave(p, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_Top1000Genes_Trajectory_MClust_new.pdf", sep="/"), width = 8, height = 8, units = "in")

# Run Cardelino on Canopy tree.
load("data/HGSOC_SS3/canopy/canopy_postmcmc_image.rda")

sc$b <- sc$d - sc$a
ids <- rownames(output.tree$Z)
sc$ID <- factor(sc$ID, levels = ids)
head(sc)
cells <- unique(sc$Cell)
n_cells <- length(cells)
n_snvs <- length(ids)

# Prepare a matrix of variant reads and total reads for each cell.
# Matrix dimension is N x C, where N is the number of mutations and C is the number of cells.
B <- matrix(0, nrow = n_snvs, ncol = n_cells)
D <- matrix(0, nrow = n_snvs, ncol = n_cells)
rownames(B) <- ids
rownames(D) <- ids
colnames(B) <- cells
colnames(D) <- cells

# Do a simple loop since the number of cells isn't that large.
for (cell_idx in 1:n_cells) {
  temp <- subset(sc, Cell == cells[cell_idx])
  B[temp$ID,cell_idx] <- temp$b
  D[temp$ID,cell_idx] <- temp$d
}

cell <- "c191"
subset(sc, Cell == cell)
B[B[,cell] != 0,cell]
D[D[,cell] != 0,cell]

ape::write.tree(output.tree, "_figures/HGSOC/canopy.nwk")

# Run cardelino:
library(cardelino)
assignments <- clone_id(B, D, Config = output.tree$Z)
names(assignments)
prob_heatmap(assignments$prob)
cluster_idx <- assign_cells_to_clones(assignments$prob)
cluster_assignment <- cluster_idx$clone
table(cluster_assignment)
cell_cardelino <- data.frame(Cell = colnames(B), Node = cluster_assignment)

# Order the cells by clones.
cells <- cell_cardelino$Cell
nodes <- unique(cell_cardelino$Node)
nodes[order(nchar(nodes), nodes)]
nodes <- nodes[order(nchar(nodes), nodes)]
cell_ordering <- order(nchar(cell_cardelino$Node), cell_cardelino$Node)
cell_order <- cells[cell_ordering]
cell_cardelino[cell_ordering,]

# Order the ids:
ids <- tree$ID
tree$CloneName <- factor(tree$CloneName, levels = c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "C_D", "C", "E_F_G_H_I", "E_F"))
id_order <- ids[order(tree$CloneName)]

# Make a heatmap plot of the single cells and their mutation status.
sc_cardelino <- left_join(sc, cell_cardelino)
head(sc_cardelino)
sc_cardelino$Cell <- factor(sc_cardelino$Cell, levels = cell_order)
sc_cardelino <- subset(sc_cardelino, ID %in% ids)
sc_cardelino$ID <- factor(sc_cardelino$ID, levels = id_order)
sc_cardelino$Node <- factor(sc_cardelino$Node, levels = nodes)

snv_cell_coverage <- sc_cardelino %>% group_by(ID) %>% summarise(n = sum(b > 0))
snv_cell_coverage_ <- snv_cell_coverage[which(snv_cell_coverage$n > 0),]
sc_cardelino_ <- subset(sc_cardelino, ID %in% snv_cell_coverage_$ID & b >= 1)
dim(sc_cardelino)
dim(sc_cardelino_)

# Make a heatmap.
names(sc_cardelino_)
base_size <- 12
plot.df <- subset(sc_cardelino_, b > 0 & Node != "unassigned")
p <- ggplot(plot.df, aes(ID, Cell, fill=Node)) + geom_tile(colour = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cells")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + scale_fill_manual(values = canopy_clone_colors)
p

cell_count <- length(unique(plot.df$Cell))
# Ancestral
p_annot <- p + annotate("rect", xmin = which(id_order == "s0") - 0.5, xmax = which(id_order == "s19") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "gray", alpha = 0.2)
# ABCD
p_annot <- p_annot + annotate("rect", xmin = which(id_order == "s2") - 0.5, xmax = which(id_order == "s41") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "red", alpha = 0.2)
# CD
p_annot <- p_annot + annotate("rect", xmin = which(id_order == "s33") - 0.5, xmax = which(id_order == "s63") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "blue", alpha = 0.2)
# C
p_annot <- p_annot + annotate("rect", xmin = which(id_order == "s37") - 0.5, xmax = which(id_order == "s37") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "blue", alpha = 0.2)
# EFGHI
p_annot <- p_annot + annotate("rect", xmin = which(id_order == "s53") - 0.5, xmax = which(id_order == "s53") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "green", alpha = 0.2)
# EF
p_annot <- p_annot + annotate("rect", xmin = which(id_order == "s10") - 0.5, xmax = which(id_order == "s66") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "green", alpha = 0.2)
ggsave("_figures/HGSOC/canopy_cardelino_sc_coclustering_using_gt_ordering.pdf", p_annot, width = 6, height = 12, units = "in")

# Now plot using mutation order inferred from Canopy.
names(sc_cardelino_)
base_size <- 12
plot.df <- subset(sc_cardelino_, b > 0 & Node != "unassigned")
ids <- rownames(output.tree$Z)
values <- strtoi(apply(output.tree$Z, 1, function(row) { paste(row, collapse = "") }), base = 2)
ids <- ids[order(values, decreasing = T)]
ids <- ids[ids %in% plot.df$ID]
plot.df$ID <- factor(plot.df$ID, levels = ids)
p <- ggplot(plot.df, aes(ID, Cell, fill=Node)) + geom_tile(colour = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cells")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + scale_fill_manual(values = canopy_clone_colors)
p
ggsave("_figures/HGSOC/canopy_cardelino_sc_coclustering.pdf", p, width = 6, height = 12, units = "in")
output.tree$Z[ids,]

p_annot <- p + annotate("rect", xmin = which(ids == "s33") - 0.5, xmax = which(ids == "s53") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[7], alpha = 0.4)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s3") - 0.5, xmax = which(ids == "s24") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[5], alpha = 0.4)
ggsave("_figures/HGSOC/canopy_cardelino_sc_coclustering_annotated.pdf", p_annot, width = 6, height = 12, units = "in")

cell_count <- length(unique(plot.df$Cell))
p_annot <- p + annotate("rect", xmin = which(ids == "s8") - 0.5, xmax = which(ids == "s64") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "white", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s4") - 0.5, xmax = which(ids == "s51") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[1], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s2") - 0.5, xmax = which(ids == "s61") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[2], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s6") - 0.5, xmax = which(ids == "s66") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[3], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s0") - 0.5, xmax = which(ids == "s56") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[4], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s20") - 0.5, xmax = which(ids == "s32") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[6], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s16") - 0.5, xmax = which(ids == "s18") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[7], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s33") - 0.5, xmax = which(ids == "s53") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[8], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s10") - 0.5, xmax = which(ids == "s63") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = canopy_clone_colors[9], alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(ids == "s42") - 0.5, xmax = which(ids == "s57") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "green", alpha = 0.2)

# Use Cardelino using PhylEx tree.
output.tree$Z
phylex_config <- matrix(0, ncol = length(table(datum2node$Node)), nrow = dim(datum2node)[1])
phylex_config

# Run InferCNV.
ANNOTATION_FILE <- "data/HGSOC_SS3/annotation_denovo.txt"
GENE_ORDER_FILE <- "data/HGSOC_SS3/gene_order.txt"
MIN_CELLS <- 5
seurat <- CreateSeuratObject(counts = feature_counts, min.cells = MIN_CELLS)
exp.rawdata <- as.matrix(seurat@assays$RNA@counts)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=exp.rawdata,
                                     ref_group_names=NULL,
                                     annotations_file = ANNOTATION_FILE,
                                     gene_order_file=GENE_ORDER_FILE) 

infercnv_obj_hmm_i6_subclusters <- infercnv::run(infercnv_obj,
                                                 cutoff=1, 
                                                 out_dir="data/HGSOC_SS3/infercnv/hmm_i6_subclusters/", 
                                                 cluster_by_groups=FALSE, 
                                                 denoise=TRUE,
                                                 HMM=TRUE,
                                                 HMM_type = "i6",
                                                 analysis_mode = "subclusters")

infercnv_cnv_pred_genes <- read.table("data/HGSOC_SS3/infercnv/hmm_i6_subclusters/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat", header = T)
infercnv_cnv_pred_regions <- read.table("data/HGSOC_SS3/infercnv/hmm_i6_subclusters/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat", header = T)
infercnv_cell_groupings <- read.table("data/HGSOC_SS3/infercnv/hmm_i6_subclusters/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings", header = T)
infercnv_genes_used <- read.table("data/HGSOC_SS3/infercnv/hmm_i6_subclusters/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.genes_used.dat", header = T)
dim(infercnv_cnv_pred_genes)
dim(infercnv_cnv_pred_regions)
table(infercnv_cnv_pred$cell_group_name)
head(infercnv_cnv_pred_regions)
head(infercnv_cnv_pred_genes)
unique(infercnv_cnv_pred_regions$cell_group_name)

# Show that the clustering does not work?
library(Seurat)
seurat <- CreateSeuratObject(counts = exp.rawdata)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)

mean(cell.df$SampleName == colnames(seurat))
Idents(seurat) <- cell.df$Clone
DimPlot(seurat)

seurat$SampleName <- colnames(seurat)
seurat_ <- subset(seurat, SampleName %in% colnames(infercnv_obj_hmm_i6_subclusters@expr.data))
match_idx <- match(seurat_$SampleName, infercnv_cell_groupings$cell)
mean(seurat_$SampleName == infercnv_cell_groupings$cell[match_idx])
infercnv_cell_groupings_ <- infercnv_cell_groupings[match_idx,]
mean(infercnv_cell_groupings_$cell == seurat_$SampleName)
Idents(seurat_) <- infercnv_cell_groupings_$cell_group_name
DimPlot(seurat_)

# ZINB-WaVE plots.
zinb_ <- zinb[,zinb$CellName %in% colnames(infercnv_obj_hmm_i6_subclusters@expr.data)]
match_idx <- match(zinb_$CellName, infercnv_cell_groupings$cell)
mean(zinb_$CellName == infercnv_cell_groupings$cell[match_idx])
infercnv_cell_groupings_ <- infercnv_cell_groupings[match_idx,]
mean(infercnv_cell_groupings_$cell == zinb_$CellName)
W <- reducedDim(zinb_)
zinb_$InferCNVClone <- infercnv_cell_groupings_$cell_group_name
df <- data.frame(W1=W[,1], W2=W[,2], Clone=zinb_$Clone, InferCNVClone=zinb_$InferCNVClone)
p <- ggplot(df, aes(W1, W2, colour=InferCNVClone)) + geom_point(alpha=0.8) + theme_classic()
p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text = element_text(size=base_size*2), legend.title = element_text(size=base_size*2))
p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
p <- p + xlab("Dim1") + ylab("Dim2")
#p <- p + ggtitle(plot_title)
p <- p + theme(title = element_text(size = base_size * 2))
#p <- p + scale_color_manual(values = clone_colors)
p

# Let's look at the mutation plots for the CNVs.
sc_ <- subset(sc, SampleName %in% infercnv_cell_groupings$cell)
length(unique(sc_$SampleName))
names(infercnv_cell_groupings) <- c("CNVClone", "SampleName")
sc_ <- left_join(sc_, infercnv_cell_groupings)
sc_$CNVClone <- gsub(pattern = "all_observations.all_observations.", "", sc_$CNVClone)
sc_$CNVClone <- gsub(pattern = "\\.", replacement = "_", sc_$CNVClone)

sc_ <- left_join(sc_, cell.df)
sc_ <- subset(sc_, ID %in% tree$ID)
ids <- unique(sc_$ID)
length(ids)
cells_clustered <- sc_[order(sc_$CNVClone),"Cell"]
sc_$Cell <- factor(sc_$Cell, levels = unique(cells_clustered))
tree$CloneName <- factor(tree$CloneName, levels = c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "C_D", "C", "E_F_G_H_I", "E_F"))
mut_ids_clustered <- tree[order(tree$CloneName), "ID"]
mut_ids_clustered <- mut_ids_clustered[mut_ids_clustered %in% ids]
sc_$ID <- factor(sc_$ID, levels = mut_ids_clustered)
sc_$b <- sc_$d - sc_$a
base_size <- 8
p <- ggplot(subset(sc_, b > 0), aes(ID, Cell, fill = CNVClone)) + geom_tile(color = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cell")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x = element_text(size = base_size))
p <- p + theme(axis.title.y = element_text(size = base_size))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.text = element_text(size = 7, face = "bold"))
p <- p + guides(fill = guide_legend(title = "CNV Clone", size = 7))
p

sc_dat <- subset(sc_, b > 0)
cell_count <- length(unique(sc_dat$Cell))
tree[order(tree$CloneName),]
p_annot <- p + annotate("rect", xmin = which(mut_ids_clustered == "s0") - 0.5, xmax = which(mut_ids_clustered == "s19") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "gray", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(mut_ids_clustered == "s2") - 0.5, xmax = which(mut_ids_clustered == "s41") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "red", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(mut_ids_clustered == "s33") - 0.5, xmax = which(mut_ids_clustered == "s63") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "blue", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(mut_ids_clustered == "s37") - 0.5, xmax = which(mut_ids_clustered == "s37") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "blue", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(mut_ids_clustered == "s53") - 0.5, xmax = which(mut_ids_clustered == "s53") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "green", alpha = 0.2)
p_annot <- p_annot + annotate("rect", xmin = which(mut_ids_clustered == "s10") - 0.5, xmax = which(mut_ids_clustered == "s66") + 0.5, ymin = 0, ymax = cell_count + 0.5, fill = "green", alpha = 0.2)
p_annot
ggsave(p_annot, filename = paste(ANALYSIS_OUTPUT_PATH, "infer_cnv_clones_heatmap.pdf", sep="/"), width=8, height = 12, units = "in")

infercnv_obj_hmm_i6_cells <- infercnv::run(infercnv_obj,
                                                 cutoff=1, 
                                                 out_dir="data/HGSOC_SS3/infercnv/hmm_i6_cells/", 
                                                 cluster_by_groups=FALSE, 
                                                 denoise=TRUE,
                                                 HMM=TRUE,
                                                 HMM_type = "i6",
                                                 analysis_mode = "cells")

infercnv_obj_hmm_i6_samples <- infercnv::run(infercnv_obj,
                                           cutoff=1, 
                                           out_dir="data/HGSOC_SS3/infercnv/hmm_i6_samples/", 
                                           cluster_by_groups=FALSE, 
                                           denoise=TRUE,
                                           HMM=TRUE,
                                           HMM_type = "i6",
                                           analysis_mode = "samples")

# Read CNV files from Laks.
laks_cn <- read.table("~/data/cell-line/bulk/OV2295/ov2295_clone_cn.csv", header = T, sep=",")
laks_cn <- subset(laks_cn, clone_id %in% c("A", "B", "C", "D", "E", "F"))
head(laks_cn)
dim(infercnv_cnv_pred_regions)
table(infercnv_cell_groupings$CNVClone)
table(infercnv_cnv_pred_regions$cell_group_name)

# For each region showing up in infercnv_cnv_pred_regions, let's compare the copy numbers.
# We don't know which InferCNV clone corresponds to which Laks clone.
infer_cnv_state_to_cn <- c(0, 1, 2, 3, 4, 6)
cell_group_names <- unique(infercnv_cnv_pred_regions$cell_group_name)
laks_cn.gr <- ConstructGranges(paste("chr", laks_cn$chr, sep=""), start = laks_cn$start, width = laks_cn$end - laks_cn$start)

# Plot the CNV form Laks vs CNV states from InferCNV.
# We will plot Laks CN for each of the InferCNV clones.
library(ggpubr)
for (j in 1:length(cell_group_names))
{
  cell_group_name_ <- cell_group_names[j]
  temp <- subset(infercnv_cnv_pred_regions, cell_group_name == cell_group_name_)
  temp$cn <- infer_cnv_state_to_cn[temp$state]
  temp.gr <- ConstructGranges(temp$chr, start = temp$start, width = temp$end - temp$start)
  for (i in 1:dim(temp)[1]) {
    ret <- findOverlaps(temp.gr[i], laks_cn.gr)
    laks_temp <- laks_cn[ret@to,]
    y_ticks <- 0:max(c(laks_temp$total_cn, max(temp$cn)))
    pl1 <- ggplot(laks_temp) + geom_segment(aes(x = start, xend = end, y=total_cn, yend = total_cn)) + facet_grid(clone_id ~ .) + scale_y_continuous(breaks = y_ticks, labels = y_ticks, limits = c(0, max(y_ticks))) + theme_bw() + ylab("Total CN") + xlab("Position") + ggtitle("Ground truth") + xlim(c(temp[i,"start"], temp[i,"end"]))
    pl2 <- ggplot(temp[i,]) + geom_segment(aes(x = start, xend = end, y=cn, yend = cn)) + scale_y_continuous(breaks = y_ticks, labels = y_ticks, limits = c(0, max(y_ticks))) + theme_bw() + ggtitle("InferCNV") + ylab("Total CN") + xlab("")
    pl <- ggarrange(pl2, pl1, nrow = 2, heights = c(1, 3))
    ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "/infercnv/", cell_group_name_, "_", temp[i,"cnv_name"], ".pdf", sep=""), height = 11, units = "in")
  }
}

# Finally, we will demonstrate the performance of PhylEx when using different hyperparameters.
# This can be executed from Jupyter notebook (see notebook/).

