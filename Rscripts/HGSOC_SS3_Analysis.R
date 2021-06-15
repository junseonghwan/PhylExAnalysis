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

DATA_PATH <- "~/PhylExAnalysis/data/"
#PHYLEX_OUTPUT_PATH <- "~/PhylExAnalysis/_output/HGSOC/phylex"
PHYLEX_OUTPUT_PATH <- "~/PhylExAnalysis/data/HGSOC_results/"
ANALYSIS_OUTPUT_PATH <- "~/PhylExAnalysis/_figures/HGSOC/"
FEATURE_COUNTS_PATH <- paste(DATA_PATH, "HGSOC_fc.txt", sep="/")
BULK_PATH <- paste(DATA_PATH, "HGSOC_bulk.txt", sep="/")
SC_PATH <- paste(DATA_PATH, "HGSOC_sc.txt", sep="/")
SC_HP_PATH <- paste(DATA_PATH, "HGSOC_sc_hp.txt", sep="/")
GT_PATH <- paste(DATA_PATH, "HGSOC_gt.txt", sep="/")

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
chains <- 0:3

# Load the data for evaluation.
dat <- read.table(BULK_PATH, header=T, sep="\t")
sc <- read.table(SC_PATH, header=T, as.is = TRUE)
sc_hp <- read.table(SC_HP_PATH, header=T)
fc <- read.table(FEATURE_COUNTS_PATH, header=T, row.names = 1, check.names = F)
valid_clones <- c("A_B_C_D_E_F_G_H_I", "A_B_C_D", "A_B", "C_D", "A", "B", "C", "D", "E_F_G_H_I", "E_F", "E", "F")
gt <- read.table(GT_PATH, header=T, as.is = T)
validation_idx <- which(gt$CloneName %in% valid_clones)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Find the best replicate to use for evaluation.
# best_chain <- FindBestRep(chain_paths = PHYLEX_OUTPUT_PATH, chains = chains)
# best_chain_path <- paste(PHYLEX_OUTPUT_PATH, "/chain", best_chain, sep="")
# datum2node <- read.table(paste(best_chain_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)

datum2node <- read.table(paste(PHYLEX_OUTPUT_PATH, "/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)

temp <- left_join(datum2node, gt)
names(temp) <- c("ID", "Node", "CloneName")
temp <- temp[validation_idx,]
head(temp)
tree <- temp[order(nchar(temp$Node), temp$Node),]
table(tree$CloneName)
table(tree$Node)
# We can map the inferred nodes to ground truth clones:
# 0_0: Ancestral
# 0_0_0: ABCD
# 0_0_1: EF/EFGHI
# 0_0_0_0: CD
# 0_0_0_0_0: C
# Note: node names may vary depending on the number of chains and length of MCMC iterations.
unique(tree[,-1])

cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = FALSE)
pl <- PlotTotalCounts(sc)
pl_ <- CoclusteringPlot(sc, cell.df)
ggsave(pl_, filename = paste(ANALYSIS_OUTPUT_PATH, "coclustering.pdf", sep="/"), height = 8, width = 3.5, units = "in")

# For the purposes of DGE, we will merge any node with MIN_SNV_COUNT or fewer into its parent.
datum2node_ <- CollapseClones(datum2node, MIN_SNV_COUNT)
table(datum2node_$Node)

# Assign cell to nodes.
cell.df <- AssignCellsBursty(sc, datum2node_, bursty_hp, sc_hp, include_normal_clone = FALSE)
cell.df$Node <- as.character(cell.df$Node)
cell.df$Clone <- cell.df$Node
table(cell.df$Node)

# Analyze the gene expression data
feature_counts <- fc[,names(fc) %in% unique(sc$Cell)]
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
dim(feature_counts)

sce <- SingleCellExperiment(assays = list(counts = as.matrix(feature_counts)),
                            rowData = data.frame(EnsemblID = rownames(feature_counts)),
                            colData = data.frame(CellName = cell.df$Cell, Node=cell.df$Node, Clone=cell.df$Clone))
filter <- rowSums(assay(sce) > MIN_READS) > MIN_CELLS
table(filter)
sce <- sce[filter,]
sce_topK <- SelectTopKGenes(sce, topK)

# Run ZINB-WaVE
sce_zinb <- zinbwave(sce_topK, K = 2, normalizedValues=TRUE, observationalWeights = TRUE)
W <- reducedDim(sce_zinb)
df <- data.frame(W1=W[,1], W2=W[,2], Node=colData(sce_topK)$Node, Clone=colData(sce_topK)$Clone)
pl <- ReducedDimensionPlots(df)

# Plot without the 0_0 node
pl_ <- ReducedDimensionPlots(subset(df, !(Clone %in% c("0_0") )))

# tSNE
set.seed(123455)
tsne_data <- Rtsne::Rtsne(t(counts(sce_topK)), pca = TRUE, perplexity=30, max_iter=5000)
df <- data.frame(W1=tsne_data$Y[,1], W2=tsne_data$Y[,2], Node=colData(sce_topK)$Node, Clone=colData(sce_topK)$Clone)
pl <- ReducedDimensionPlots(df)

# UMAP
logcounts(sce_topK) <- log1p(counts(sce_topK))
sce_umap <- runUMAP(sce_topK, name = "UMAP")
W <- reducedDim(sce_umap, "UMAP")
df <- data.frame(W1=W[,1], W2=W[,2], Node=colData(sce_topK)$Node, Clone=colData(sce_topK)$Clone)
pl <- ReducedDimensionPlots(df)

# Perform DGE: https://bioconductor.org/packages/devel/bioc/vignettes/zinbwave/inst/doc/intro.html#differential-expression.
cols = c("black", "red")

weights <- assay(sce_zinb, "weights")
dge <- DGEList(assay(sce_zinb), group = sce_zinb$Node, remove.zeros = TRUE)
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + Node, data = colData(sce_zinb))
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)

lrt <- glmWeightedF(fit, contrast = c(-1, 1, 0, 0))
dge_results <- topTags(lrt, topK)
sig_genes <- rownames(dge_results$table[dge_results$table$FDR < FDR_THRESHOLD,])
gene_names <- ConvertToGeneName(rownames(dge_results), mart)
gene_names <- gene_names[match(rownames(dge_results), gene_names$ensembl_gene_id),]
mean(gene_names$ensembl_gene_id == rownames(dge_results))
volcano.df <- data.frame(hgnc_symbol=as.character(gene_names$external_gene_name), logfc=dge_results$table$logFC, logpvalue=-log(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)
pl <- MakeVolcanoPlot(volcano.df, "Ancestral vs EF")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "/Ancestral_vs_EF_Volcano.pdf", sep=""), width = 8, height = 8, units = "in")

lrt <- glmWeightedF(fit, contrast = c(0, 1, -1, 0))
dge_results <- topTags(lrt, topK)
sig_genes <- rownames(dge_results$table[dge_results$table$FDR < FDR_THRESHOLD,])
gene_names <- ConvertToGeneName(rownames(dge_results), mart)
gene_names <- gene_names[match(rownames(dge_results), gene_names$ensembl_gene_id),]
mean(gene_names$ensembl_gene_id == rownames(dge_results))
volcano.df <- data.frame(hgnc_symbol=as.character(gene_names$external_gene_name), logfc=dge_results$table$logFC, logpvalue=-log(dge_results$table$PValue), significant = dge_results$table$FDR < FDR_THRESHOLD)
pl <- MakeVolcanoPlot(volcano.df, "ABCD vs EF")
ggsave(pl, filename = paste(ANALYSIS_OUTPUT_PATH, "/ABCD_vs_EF_Volcano.pdf", sep=""), width = 8, height = 8, units = "in")

# GSEA: https://www.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html.
# We will use as many genes as possible rather than top 1000 variable genes.
# This article suggests to normalize the raw counts.
counts <- assay(sce)

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

dge <- DGEList(counts_, group = sce$Node, remove.zeros = TRUE)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + Node, data = colData(sce))
colnames(design)
dge <- estimateDisp(dge, design)

lcpm_ <- edgeR::cpm(dge, log=TRUE)
dim(lcpm_)

contr.matrix <- makeContrasts(AncestralvsABCD = Node0_0_1 - Node0_0,
                              AncestralvsEF = Node0_0_0 - Node0_0,
                              ABCDvsEF = Node0_0_0 - Node0_0_1,
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
load(paste(DATA_PATH, "human_c5_v5p2.rdata", sep="/"))
indices <- ids2indices(Hs.c5, rownames(dge))

# Ancestral vs ABCD
test_ret <- camera(v, indices, contrast=contr.matrix[,1])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))

# Ancestral vs EF
test_ret <- camera(v, indices, contrast=contr.matrix[,2])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))

# ABCD vs EF
test_ret <- camera(v, indices, contrast=contr.matrix[,3])
cbind(head(rownames(test_ret), 20), head(test_ret$FDR, 20))

