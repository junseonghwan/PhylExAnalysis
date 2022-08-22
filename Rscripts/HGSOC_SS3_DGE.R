# Replace HGSOC_SS3_DGE with this script.
rm(list=ls())
library(apeglm)
library(biomaRt)
library(BiocParallel)
library(DESeq2)
library(edgeR)
library(gam)
library(ggplot2)
library(Glimma)
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

DATA_PATH <- "data/HGSOC_SS3/"
FEATURE_COUNTS_FILE <- "data/HGSOC_SS3/featureCounts.txt"
BULK_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "sc.txt", sep="/")
SC_HP_PATH <- paste(DATA_PATH, "sc_hp.txt", sep="/")
GT_PATH <- paste(DATA_PATH, "gt.txt", sep="/")
ANALYSIS_OUTPUT_PATH <- "_figures/HGSOC/"

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
unique(datum2node[validation_idx,"Node"])
temp <- cbind(datum2node[validation_idx,], gt[validation_idx,"CloneName"])
names(temp) <- c("ID", "Node", "CloneName")
tree <- temp[order(temp$Node),]
tree

# Assign cell to nodes.
cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = FALSE)
#write.table(cell.df, paste(DATA_PATH, "/results/rep", best_rep, "/joint/tree0/cell2node.tsv", sep=""), quote = F, col.names = T, row.names = F, sep="\t")
cell.df$SampleName <- unique(sc$SampleName)
cell.df$Node <- as.character(cell.df$Node)
cell.df$Clone <- cell.df$Node
table(cell.df$Clone)
cell.df$Clone[startsWith(cell.df$Clone, "0_0_0_0")] <- "CD"
cell.df$Clone[startsWith(cell.df$Clone, "0_0_0")] <- "ABCD"
cell.df$Clone[cell.df$Clone == "0_0_1"] <- "EF"
cell.df$Clone[cell.df$Clone == "0_0"] <- "Ancestral"
cell.df$Clone <- factor(cell.df$Clone, levels = c("Ancestral", "ABCD", "CD", "EF"))
table(cell.df$Clone)

# Separate to 3 groups for better comparison.
cell.df$Group <- cell.df$Node
cell.df$Group[cell.df$Group == "0_0_1"] <- "EF"
cell.df$Group[startsWith(cell.df$Group, "0_0_0")] <- "ABCD"
cell.df$Group[cell.df$Group == "0_0"] <- "Ancestral"
cell.df$Group <- factor(cell.df$Group, levels = c("Ancestral", "ABCD", "EF"))

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
                            colData = data.frame(CellName = cell.df$SampleName, Node=cell.df$Node, Clone=cell.df$Clone, Group=cell.df$Group))
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

clone_colors <- c("Ancestral"="#999999", "ABCD"="#E69F00", "CD"="#56B4E9", "EF"="#009E73")

# df: Data frame with columns W1, W2, Node, and clone.
ReducedDimensionPlots <- function(df, plot_title, base_size = 12) {
  #p <- ggplot(df, aes(W1, W2, colour=Clone, shape=Clone)) + geom_point(alpha=0.8) + theme_classic()
  p <- ggplot(df, aes(W1, W2, colour=Clone)) + geom_point(alpha=0.8) + theme_classic()
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
  p <- p + theme(legend.text = element_text(size=base_size*2), legend.title = element_text(size=base_size*2))
  p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  p <- p + xlab("Dim1") + ylab("Dim2")
  #p <- p + ggtitle(plot_title)
  p <- p + theme(title = element_text(size = base_size * 2))
  p <- p + scale_color_manual(values = clone_colors)
  return(p)
}

# 2D ZINB-WaVE.
zinb <- zinbwave(sce, K = 2, epsilon=1000)
W <- reducedDim(zinb)
df <- data.frame(W1=W[,1], W2=W[,2], Node=colData(sce)$Node, Clone=colData(sce)$Clone)
zinb_df <- df
pl <- ReducedDimensionPlots(df, "ZINB-WaVE dimensions")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_Top1000Genes.pdf", sep="/"), width=8, height=8, units = "in")

# Save source data for Figure 3c, Supplementary Figure 3c, Supplementary Figure 4a-d.
write.csv(df, file = "data/NatComm/Figure3c_SupplementaryFigure3c_SupplementaryFigure4a-d.csv", quote = F, row.names = F)

# Remove the ancestral clone.
df_ <- subset(df, !(Clone %in% c("Ancestral") ))
pl_ <- ReducedDimensionPlots(df_, "ZINB-WaVE dimensions")
ggsave(pl_, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_Top1000Genes_WithoutAncestralClone.pdf", sep="/"), width=8, height=8, units = "in")

# t-SNE on topK genes.
set.seed(123455)
tsne_data <- Rtsne::Rtsne(t(counts(sce)), pca = TRUE, perplexity=30, max_iter=5000)
df <- data.frame(W1=tsne_data$Y[,1], W2=tsne_data$Y[,2], Clone=colData(sce)$Clone)
tsne_df <- df
pl <- ReducedDimensionPlots(df, "t-SNE dimensions")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "tSNE_Top1000Genes.pdf", sep="/"), width = 8, height = 8, units = "in")

# Save source data for Figure 3d, Supplementary Figure 3d
write.csv(df, file = "data/NatComm/SupplementaryFigure3d.csv", quote = F, row.names = F)

# Remove the ancestral clone.
df_ <- subset(df, !(Clone %in% c("Ancestral") ))
pl_ <- ReducedDimensionPlots(df_, "t-SNE dimensions")
ggsave(pl_, filename = paste(ANALYSIS_OUTPUT_PATH, "tSNE_Top1000Genes_WithoutAncestralClone.pdf", sep="/"), width = 8, height = 8, units = "in")

# Run with normalizing the values option set to TRUE.
sce_zinb <- zinbwave(sce, K = 2, epsilon=1000, normalizedValues=TRUE, observationalWeights = TRUE)
weights <- assay(sce_zinb, "weights")
W <- reducedDim(sce_zinb)
df <- data.frame(W1=W[,1], W2=W[,2], Clone=colData(sce)$Clone)
zinb_df2 <- df
pl <- ReducedDimensionPlots(df, "ZINB-WaVE dimensions")
# Not much of a difference -- if any.
pl


# Do trajectory analysis using slingshot on zinbwave dimensions.
sim <- slingshot(sce_zinb, clusterLabels = 'Clone', reducedDim = 'zinbwave', start.clus = "Ancestral", reassign = T)
lineage <- SlingshotDataSet(sim)

W <- reducedDims(sim)$zinbwave
center <- matrix(0, ncol = 2, nrow = 4)
center[1,] <- colMeans(W[sim$Clone == "Ancestral",])
center[2,] <- colMeans(W[sim$Clone == "ABCD",])
center[3,] <- colMeans(W[sim$Clone == "CD",])
center[4,] <- colMeans(W[sim$Clone == "EF",])
rownames(center) <- c("Ancestral", "ABCD", "CD", "EF")
center.df <- data.frame(W1=center[,1], W2=center[,2], Clone=rownames(center))

df <- data.frame(W1 = W[,1], W2 = W[,2], Clone=sim$Clone)
pl <- ReducedDimensionPlots(df, "ZINB-WaVE with pseudo-time analysis")
pl <- pl + annotate("point", x = center.df[,1], y = center.df[,2], fill=clone_colors, size = 3)
pl <- pl + geom_point(aes(W1, W2, color = Clone), data = center.df)
pl

# Draw trajectories.
lineages <- lineage@lineages
segment <- data.frame()
for (line in lineages) {
  for (i in 1:(length(line)-1)) {
    idx1 <- (rownames(center.df) == line[i])
    idx2 <- (rownames(center.df) == line[i+1])
    segment_row <- data.frame(x = center.df[idx1,1], y = center.df[idx1,2], xend = center.df[idx2,1], yend = center.df[idx2,2])
    segment <- rbind(segment, segment_row)
  }
}
pl <- pl + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data = segment, arrow =  arrow(length = unit(0.05, "npc")), inherit.aes = FALSE)
pl
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_Top1000Genes_Trajectory.pdf", sep="/"), width = 8, height = 8, units = "in")

# Do trajectory analysis using t-SNE dimensions rather than zinbwave dimensions.
# Update sce_zinb with the t-SNE results obtained from 1000 genes.
reducedDims(sce_zinb)$tsne <- tsne_data$Y
sim <- slingshot(sce_zinb, clusterLabels = 'Clone', reducedDim = 'tsne', start.clus = "Ancestral", reassign = T)
lineage <- SlingshotDataSet(sim)

W <- reducedDims(sim)$tsne
center <- matrix(0, ncol = 2, nrow = 4)
center[1,] <- colMeans(W[sim$Clone == "Ancestral",])
center[2,] <- colMeans(W[sim$Clone == "ABCD",])
center[3,] <- colMeans(W[sim$Clone == "CD",])
center[4,] <- colMeans(W[sim$Clone == "EF",])
rownames(center) <- c("Ancestral", "ABCD", "CD", "EF")
center.df <- data.frame(W1=center[,1], W2=center[,2], Clone=rownames(center))

df <- data.frame(W1 = W[,1], W2 = W[,2], Clone=sim$Clone)
pl <- ReducedDimensionPlots(df, "ZINB-WaVE with pseudo-time analysis")
pl <- pl + annotate("point", x = center.df[,1], y = center.df[,2], fill=clone_colors, size = 3)
pl <- pl + geom_point(aes(W1, W2, color = Clone), data = center.df)

# Draw trajectories.
lineages <- lineage@lineages
segment <- data.frame()
for (line in lineages) {
  for (i in 1:(length(line)-1)) {
    idx1 <- (rownames(center.df) == line[i])
    idx2 <- (rownames(center.df) == line[i+1])
    segment_row <- data.frame(x = center.df[idx1,1], y = center.df[idx1,2], xend = center.df[idx2,1], yend = center.df[idx2,2])
    segment <- rbind(segment, segment_row)
  }
}
pl <- pl + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data = segment, arrow =  arrow(length = unit(0.05, "npc")), inherit.aes = FALSE)
pl
ggsave(p, filename = paste(ANALYSIS_OUTPUT_PATH, "tSNE_Top1000Genes_Trajectory.pdf", sep="/"), width = 8, height = 8, units = "in")

# Do trajectory analysis using zinbwave dimensions and clustering based on zinbwave dimensions.
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
p <- ggplot(df, aes(W1, W2, color = Cluster, shape = Clone)) + geom_point() + theme_classic()
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
ggsave(p, filename = paste(ANALYSIS_OUTPUT_PATH, "ZINBWaVE_Top1000Genes_Trajectory_MClust.pdf", sep="/"), width = 8, height = 8, units = "in")

# Save source data for Figure 3d.
write.csv(df, file = "data/NatComm/Figure3d.csv", row.names = F, quote = F)

# We are going to do DGE using edgeR.
# We will use 3 groups rather than the clones.
dge <- DGEList(assay(sce_zinb), group = sce_zinb$Group, remove.zeros = TRUE)
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + Group, data = colData(sce_zinb))
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
colnames(fit$coefficients)

# We will generate plots to compare the 4 major clones: Ancestral vs ABCD vs EF vs CD.
FDR_THRESHOLD <- 0.1
FDR_THRESHOLD0 <- 0.01

cols = c("black", "red")

lrt <- glmWeightedF(fit, contrast = c(-1, 1, 0))
dge_results <- topTags(lrt, topK)
sig_genes <- rownames(dge_results$table[dge_results$table$FDR < FDR_THRESHOLD,])
volcano.df <- data.frame(logfc=dge_results$table$logFC, logpvalue=-log2(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)
p <- ggplot(volcano.df, aes(logfc, logpvalue, col = significant)) + geom_point() + theme_bw() + ylab("-Log p-value") + xlab("Log fold change") + labs(title="Ancestral vs ABCD")
p <- p + scale_color_manual(breaks=c("FALSE", "TRUE"), values=cols)
p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.position="bottom") + guides(color=guide_legend(title="Significantly differential?"))
p <- p + theme(legend.text = element_text(size=base_size), title = element_text(size=base_size))
p

p <- plot_sig_genes(dge_results, mart)
p <- p + labs(title="Ancestral vs ABCD")
p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text = element_text(size=base_size), title = element_text(size=base_size))
p
# Ancestral vs ABCD plots aren't very interesting -- we won't save these.

# Pathway analysis.
# We need to convert to Entrez gene IDs.
# We will use as many genes as possible rather than focusing on 1k genes for this analysis.
sce <- SingleCellExperiment(assays = list(counts = as.matrix(feature_counts)),
                            rowData = data.frame(EnsemblID = rownames(feature_counts)),
                            colData = data.frame(CellName = cell.df$SampleName, Node=cell.df$Node, Clone=cell.df$Clone, Group=cell.df$Group))
filter <- rowSums(assay(sce) > MIN_READS) > MIN_CELLS
sce <- sce[filter,]

# Going to use a lot more genes than 1,000 for pathway analysis for more detailed analysis.
sce_ <- zinbwave(sce, K = 2, epsilon=1000, normalizedValues=TRUE, observationalWeights = TRUE)

counts <- assay(sce_)
ensembl_gene_ids <- rownames(counts)
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
                         'start_position', 'end_position', 'percentage_gene_gc_content'),
            filters = 'ensembl_gene_id',
            values = ensembl_gene_ids,
            mart = mart, uniqueRows = T)
idx <- match(ensembl_gene_ids, bm$ensembl_gene_id)
bm <- bm[idx,]

not_na_idx <- !is.na(bm$entrezgene_id)
counts_ <- counts[not_na_idx,]
rownames(counts_) <- bm$entrezgene_id[not_na_idx]
dim(counts_)

dge <- DGEList(counts_, group = sce_$Group, remove.zeros = TRUE)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + Group, data = colData(sce_))
ww <- assay(sce_, "weights")
dge$weights <- ww[not_na_idx,]
dge <- estimateDisp(dge, design)

# We will do pathway analysis following https://www.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html.
# This article suggests to normalize the raw counts.
lcpm_ <- edgeR::cpm(dge, log=TRUE)
dim(lcpm_)
# DO NOT RUN
#glMDSPlot(lcpm_, labels = sce_$CellName, groups = sce_$Clone)

contr.matrix <- makeContrasts(AncestralvsABCD = GroupABCD - GroupAncestral,
                              AncestralvsEF = GroupEF - GroupAncestral,
                              ABCDvsEF = GroupEF - GroupABCD,
                              levels = colnames(design))
v <- voom(dge, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)

# Perform camera test on GO set.
# Download gene ontology set.
#download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p2.rdata", "human_c5_v5p2.rdata", mode = "wb")
load("data/human_c5_v5p2.rdata")
indices <- ids2indices(Hs.c5, rownames(dge))
test_ret <- camera(v, indices, contrast=contr.matrix[,1])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))

# Ancestral vs EF
test_ret <- camera(v, indices, contrast=contr.matrix[,2])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))

temp <- test_ret[test_ret$FDR < FDR_THRESHOLD,]
# we will need to keep track of down regulated genes for plotting volcano plots.
down <- subset(temp, Direction == "Down")
x <- data.frame("GeneOntology"=rownames(down), "P Value"=down$PValue, FDR=down$FDR)
x$GeneOntology <- gsub("GO_", "", x$GeneOntology)
x$GeneOntology <- gsub("_", " ", x$GeneOntology)
write.table(x, file = paste(ANALYSIS_OUTPUT_PATH, "Ancestral_EF_Down.csv", sep="/"), quote = F, row.names = F, col.names = T, sep=",")
print(xtable(x, type = "latex", display=c("d", "s", "e", "e")), file = paste(ANALYSIS_OUTPUT_PATH, "Ancestral_EF_Down.tex", sep="/"),
      include.rownames = FALSE, size=c("tiny"))

up <- subset(temp, Direction == "Up")
x <- data.frame("GeneOntology"=rownames(up), "P Value"=up$PValue, FDR=up$FDR)
x$GeneOntology <- gsub("GO_", "", x$GeneOntology)
x$GeneOntology <- gsub("_", " ", x$GeneOntology)
write.table(x, file = paste(ANALYSIS_OUTPUT_PATH, "Ancestral_EF_Up.csv", sep="/"), quote = F, row.names = F, col.names = T, sep=",")
print(xtable(x, type = "latex", display=c("d", "s", "e", "e")), file = paste(ANALYSIS_OUTPUT_PATH, "Ancestral_EF_Up.tex", sep="/"),
      include.rownames = FALSE, size=c("tiny"))

# Assume volcano.df is already ordered based on significance
MakeVolcanoPlot <- function(volcano.df, plot_title, num_genes_to_label = 10, base_size = 12) {
  p <- ggplot(volcano.df, aes(logfc, logpvalue, col = significant)) + geom_point() + theme_bw() + ylab("-Log p-value") + xlab("Log fold change")
  p <- p + labs(title=plot_title)
  p <- p + scale_color_manual(breaks=c("FALSE", "TRUE"), values=cols)
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
  p <- p + theme(legend.position="bottom") + guides(color=guide_legend(title="Significantly differential?"))
  p <- p + theme(legend.text = element_text(size=base_size*2), title = element_text(size=base_size*2))
  p <- p + theme(axis.text.x = element_text(size = base_size * 2), axis.text.y = element_text(size = base_size * 2))
  annotate.df <- volcano.df[1:num_genes_to_label,]
  temp.df <- volcano.df[order(volcano.df$logfc),]
  annotate.df2 <- temp.df[1:num_genes_to_label,]
  annotate.df <- rbind(annotate.df, annotate.df2)
  annotate.df <- unique(annotate.df)
  p <- p + ggrepel::geom_label_repel(data = annotate.df, aes(label=gene_name),
                                     size = 8,
                                     nudge_x = 0.2,
                                     col = "black",
                                     segment.size  = 1,
                                     segment.color = "grey50",
                                     direction     = "both", max.overlaps = 10)
  return(p)
}

# Let's make the volcano plots.
#temp <- Hs.c5[names(Hs.c5) %in% rownames(down)]
#entrez_ids <- unique(unlist(temp))

# bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
#                          'start_position', 'end_position', 'percentage_gene_gc_content'),
#             filters = 'entrezgene_id',
#             values = entrez_ids,
#             mart = mart, uniqueRows = T)

lrt <- glmWeightedF(fit, contrast = c(-1, 0, 1))
dge_results <- topTags(lrt, topK)
#sig_genes <- rownames(dge_results$table[dge_results$table$FDR < FDR_THRESHOLD,])
#gene_symbol <- convert_to_gene_name(sig_genes, mart)
#sig_genes_to_annotate <- intersect(unique(gene_symbol$external_gene_name), unique(bm$hgnc_symbol))
gene_names <- convert_to_gene_name(rownames(dge_results), mart)
gene_names <- gene_names[match(rownames(dge_results), gene_names$ensembl_gene_id),]
mean(gene_names$ensembl_gene_id == rownames(dge_results))

volcano.df <- data.frame(gene_name=as.character(gene_names$external_gene_name), logfc=dge_results$table$logFC, logpvalue=-log(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)
#annotate1 <- (volcano.df$gene_name %in% sig_genes_to_annotate)
#annotate2 <- (abs(volcano.df$logfc) > 3 | volcano.df$logpvalue > 40)
#volcano.df$annotate <- annotate1 & annotate2
volcano.df$annotate <- (abs(volcano.df$logfc) > 3 | volcano.df$logpvalue > 40)
subset(volcano.df, annotate == TRUE)
sum(volcano.df$annotate)

pl <- MakeVolcanoPlot(volcano.df, "Ancestral vs EF")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "Laks_Ancestral_vs_EF_Volcano.pdf", sep="/"), width=8, height=8, units = "in")

# Save source data for Figure 3e.
write.csv(volcano.df, file = "data/NatComm/Figure3e.csv", row.names = F, quote = F)

# ABCD vs EF GSEA table
test_ret <- camera(v, indices, contrast=contr.matrix[,3])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))
temp <- test_ret[test_ret$FDR < FDR_THRESHOLD0,]
# we will need to keep track of down regulated genes for plotting volcano plots.
down <- subset(temp, Direction == "Down")
x <- data.frame("GeneOntology"=rownames(down), "P Value"=down$PValue, FDR=down$FDR)
x$GeneOntology <- gsub("GO_", "", x$GeneOntology)
x$GeneOntology <- gsub("_", " ", x$GeneOntology)
write.table(x, file = paste(ANALYSIS_OUTPUT_PATH, "ABCD_EF_Down.csv", sep="/"), quote = F, row.names = F, col.names = T, sep=",")
print(xtable(x, type = "latex", display=c("d", "s", "e", "e")), file = paste(ANALYSIS_OUTPUT_PATH, "ABCD_EF_Down.tex", sep="/"),
      include.rownames = FALSE, size=c("tiny"))

up <- subset(temp, Direction == "Up")
x <- data.frame("GeneOntology"=rownames(up), "P Value"=up$PValue, FDR=up$FDR)
x$GeneOntology <- gsub("GO_", "", x$GeneOntology)
x$GeneOntology <- gsub("_", " ", x$GeneOntology)
head(x)
write.table(x, file = paste(ANALYSIS_OUTPUT_PATH, "ABCD_EF_Up.csv", sep="/"), quote = F, row.names = F, col.names = T, sep=",")
print(xtable(x, type = "latex", display=c("d", "s", "e", "e")), file = paste(ANALYSIS_OUTPUT_PATH, "ABCD_EF_Up.tex", sep="/"),
      include.rownames = FALSE, size=c("tiny"))

#temp <- Hs.c5[names(Hs.c5) %in% rownames(down)]
#entrez_ids <- unique(unlist(temp))

# bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
#                          'start_position', 'end_position', 'percentage_gene_gc_content'),
#             filters = 'entrezgene_id',
#             values = entrez_ids,
#             mart = mart, uniqueRows = T)

lrt <- glmWeightedF(fit, contrast = c(0, -1, 1))
dge_results <- topTags(lrt, topK)
#sig_genes <- rownames(dge_results$table[dge_results$table$FDR < FDR_THRESHOLD,])
#gene_symbol <- convert_to_gene_name(sig_genes, mart)
#sig_genes_to_annotate <- intersect(unique(gene_symbol$external_gene_name), unique(bm$hgnc_symbol))
gene_names <- convert_to_gene_name(rownames(dge_results), mart)
gene_names <- gene_names[match(rownames(dge_results), gene_names$ensembl_gene_id),]
mean(gene_names$ensembl_gene_id == rownames(dge_results))

volcano.df <- data.frame(gene_name=as.character(gene_names$external_gene_name), logfc=dge_results$table$logFC, logpvalue=-log(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)
#annotate1 <- (volcano.df$gene_name %in% sig_genes_to_annotate)
#annotate2 <- (abs(volcano.df$logfc) > 3 | volcano.df$logpvalue > 40)
#volcano.df$annotate <- annotate1 & annotate2

pl <- MakeVolcanoPlot(volcano.df, "ABCD vs EF")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "Laks_ABCD_vs_EF_Volcano.pdf", sep="/"), width=8, height=8, units = "in")

# Source data for Figure 3f.
write.csv(volcano.df, file = "data/NatComm/Figure3f.csv", row.names = F, quote = F)

# Ancestral vs ABCD
lrt <- glmWeightedF(fit, contrast = c(-1, 1, 0))
dge_results <- topTags(lrt, topK)
gene_names <- convert_to_gene_name(rownames(dge_results), mart)
gene_names <- gene_names[match(rownames(dge_results), gene_names$ensembl_gene_id),]
mean(gene_names$ensembl_gene_id == rownames(dge_results))

volcano.df <- data.frame(gene_name=as.character(gene_names$external_gene_name), logfc=dge_results$table$logFC, logpvalue=-log(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)

pl <- MakeVolcanoPlot(volcano.df, "Ancestral vs ABCD")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "Laks_Ancestral_vs_ABCD_Volcano.pdf", sep="/"), width=8, height=8, units = "in")
