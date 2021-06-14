rm(list=ls())
library(apeglm)
library(biomaRt)
library(BiocParallel)
library(DESeq2)
library(edgeR)
library(gam)
library(ggplot2)
library(ggrepel)
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
library(zinbwave)
#source("R/CellAssign.R")
#source("TNBC/Figure/Common.R")

FEATURE_COUNTS_FILE <- "../data/HER2_fc.txt"
BULK_DATA_FILE <- "../data/HER2_bulk.txt"
SC_DATA_FILE <- "../data/HER2_sc.txt"
SC_HP_FILE <- "../data/HER2_sc_hp.txt"
MIN_CELLS <- 2
MIN_READS <- 5

PHYLEX_OUT_PATH <- "../_output/HER2/"

reps <- 0:3
REP_COUNT <- length(reps)

REGION_COUNT <- 5

topK <- 1000
base_size <- 12

shapes <- c(15, 16, 17, 8, 3)
chrs <- c(1:22, "X", "Y")

REGIONS <- c("A", "B", "C", "D", "E")

# Get biomart.
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Perform DGE.
fc <- read.table(FEATURE_COUNTS_FILE, header=T, row.names = 1, check.names = F)
colnames(fc) <- gsub("Aligned", "", colnames(fc))

# Remove MT genes
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', "chromosome_name"),
            filters = 'ensembl_gene_id',
            values = rownames(fc),
            mart = mart, uniqueRows = T)
idx <- match(rownames(fc), bm$ensembl_gene_id)
bm <- bm[idx,]
idx <- which(bm$chromosome_name %in% chrs)
fc <- fc[idx,]
bm <- bm[idx,]
mean(bm[,1] == rownames(fc), na.rm = T)
dim(fc)

log_liks <- rep(0, REP_COUNT)
for (rep_no in reps) {
    rep_path <- paste(PHYLEX_OUT_PATH, "/phylex/chain", rep_no, sep="")
    log_liks[rep_no + 1] <- read.table(paste(rep_path, "/joint/tree0/log_lik.txt", sep=""), header=F)$V1
}
best_rep <- which.max(log_liks) - 1

rep_path <- paste(PHYLEX_OUT_PATH, "/phylex/chain", best_rep, sep="")
dat <- read.table(BULK_DATA_FILE, header=T, sep="\t")
sc <- read.table(SC_DATA_FILE, header=T, as.is = TRUE)
sc_hp <- read.table(SC_HP_FILE, header=T, sep="\t")

datum2node <- read.table(paste(rep_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(datum2node) <- c("ID", "Node")
table(datum2node$Node)

# 0_0_1 does not have any cell assigned to it. Let's see why.
bulk_temp <- datum2node[datum2node$Node == "0_0_1",]
sc_temp <- subset(sc, ID %in% bulk_temp$ID)
ret <- sc_temp %>% group_by(Cell) %>% summarise(n_b = sum(d - a > 0))
ret$n_b # All cells have at most one variant read for mutation in 0_0_1.

# Let's check 0_0_0_0_1. It gets 20 cells assigned but hard to know if there indeed is a clone with one mutation branching out from the main lineage.
bulk_temp <- datum2node[datum2node$Node == "0_0_0_0_1",]
sc_temp <- subset(sc, ID %in% bulk_temp$ID)
ret <- sc_temp %>% group_by(Cell) %>% summarise(n_b = sum(d - a > 0))
ret$n_b # Alle cells have at most one variant read for mutations in 0_0_0_0_1.

# Merge any singleton cluster into its parent
datum2node <- CollapseClones(datum2node)
table(datum2node$Node)
datum2node[datum2node$ID == "s343",]

# Assign cell to nodes.
bursty_hp <- list(alpha=0.01, beta=0.01)
cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = FALSE)
#cell.df$SampleName <- unique(sc$SampleName)
#cell.df$Node <- as.character(cell.df$Node)
#cell.df$Clone <- cell.df$Node
table(cell.df$Node)

# Collapse clones by number of cells.
ret <- CollapseClonesByCellCount(datum2node, cell.df)
datum2node <- ret$datum2node
cell.df <- ret$cell.df
table(datum2node$Node)
table(cell.df$Node)

# Re-assign cells using the collapsed clones.
#cell.df <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = F)
cell.df$SampleName <- unique(sc$SampleName)
cell.df$Node <- as.character(cell.df$Node)
cell.df$Clone <- cell.df$Node

# Set clone names for datum2node
datum2node$Clone <- datum2node$Node
nodes <- unique(datum2node$Node)
nodes_ordered <- nodes[order(nchar(nodes), nodes)]
clone_count <- length(nodes_ordered)
for (i in 1:clone_count) {
  idx <- datum2node$Node == nodes_ordered[i]
  datum2node$Clone[idx] <- i
}


# Since we have a lineage tree, let's re-label the clones from 1 to N.
nodes <- unique(cell.df$Node)
nodes_ordered <- nodes[order(nchar(nodes), nodes)]
clone_count <- length(nodes_ordered)
for (i in 1:clone_count) {
  idx <- cell.df$Node == nodes_ordered[i]
  cell.df$Clone[idx] <- i
}

b_alleles <- lapply(as.character(dat$b), function(row) {
    as.numeric(strsplit(row, split = ",")[[1]])
})
b_alleles <- matrix(unlist(b_alleles), byrow = T, nrow = dim(dat)[1], ncol = REGION_COUNT)
total_alleles <- lapply(as.character(dat$d), function(row) {
    as.numeric(strsplit(row, split = ",")[[1]])
})
total_alleles <- matrix(unlist(total_alleles), byrow = T, nrow = dim(dat)[1], ncol = REGION_COUNT)
vafs <- b_alleles/total_alleles
hist(rowMeans(vafs))

clone_depth <- lapply(cell.df$Clone, function(x) {
  length(as.numeric(strsplit(x, split = "_")[[1]]))
})
cell.df$Depth <- unlist(clone_depth)

fc_ <- fc[,colnames(fc) %in% cell.df$Cell]
dim(fc_)

# Sanity check: make sure we got the colnames of fc_ matching the SampleNames.
cell.df_idx <- match(colnames(fc_), cell.df$Cell)
mean(colnames(fc_) == cell.df$Cell[cell.df_idx])

sce <- SingleCellExperiment(assays = list(counts = as.matrix(fc_)),
                            rowData = data.frame(EnsemblID = rownames(fc_)),
                            colData = data.frame(CellName = cell.df$SampleName[cell.df_idx], Node=cell.df$Node[cell.df_idx], Clone=cell.df$Clone[cell.df_idx], Region=cell.df$Region[cell.df_idx],
                                                 Depth = cell.df$Depth[cell.df_idx]))

# Filter out cells without sufficient reads.
total_reads <- log1p(colSums(counts(sce)))
plot(total_reads[order(total_reads, decreasing = T)])
# Elbow at around 12
sum(total_reads >= 12)
sce <- sce[,total_reads >= 12]
log1p(colSums(counts(sce)))


# Let's select the genes in NanoString panel.
library(xlsx)
nano_string <- read.xlsx("../data/LBL-10025_nCounter_PanCancer_Human_Pathways_Panel_Gene_List.xlsx",
                         sheetIndex = 2, header = T, startRow = 2, endRow = 774)

bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', "chromosome_name", "external_gene_name"),
            filters = 'external_gene_name',
            values = nano_string$Official.Symbol,
            mart = mart, uniqueRows = T)
sum(is.na(bm$ensembl_gene_id))

nano_list_idx <- rownames(sce) %in% bm$ensembl_gene_id
sum(nano_list_idx)
sce_nano <- sce[nano_list_idx,]
row_filter <- rowSums(assay(sce_nano) > MIN_READS) > MIN_CELLS
table(row_filter)
sce_nano <- sce_nano[row_filter,]
col_filter <- colSums(assay(sce_nano)) > 0
table(col_filter)
sce_nano <- sce_nano[,col_filter]
logcounts(sce_nano) <- log1p(counts(sce_nano))


# We are going to do DGE analysis using edgeR using NanoString genes -- 602 genes.
dim(sce_nano)
table(sce_nano$Clone)
dge <- DGEList(assay(sce_nano), group = sce$Clone, remove.zeros = TRUE)
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + Clone, data = colData(sce_nano))
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
conditions <- colnames(fit$coefficients)
conditions
clones <- colnames(design)
clone_count <- length(clones)

# We will generate plots to compare the major clones.
FDR_THRESHOLD <- 0.1

cols = c("black", "red")

gene_count <- dim(sce_nano)[1]

# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# Section 2.11.

num_genes_to_label <- 10
# Remove duplicates from bm.
bm <- bm[!duplicated(paste(bm$ensembl_gene_id, bm$hgnc_symbol, sep="_")),]
dim(bm)

dest_dir <- "../_figures/HER2_pos"
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = T)
}

GetContrast <- function(idx, clone_count) {
  con <- rep(0, clone_count)
  con[idx[1]] <- -1
  con[idx[2]] <- 1
  return(con)
}

### Perform parent-child clone comparisons. ###
for (i in 1:(clone_count - 1)) {
  idx <- c(i, i+1)
  filename <- paste(clones[idx], collapse = "_vs_")
  title <- paste(clones[idx], collapse = " vs ")
  lrt <- glmLRT(fit, contrast = GetContrast(idx, clone_count))
  dge_results <- topTags(lrt, gene_count)
  
  volcano.df <- data.frame(logfc=dge_results$table$logFC,
                           logpvalue=-log2(dge_results$table$PValue),
                           significant = dge_results$table$FDR < FDR_THRESHOLD,
                           ensembl_gene_id=rownames(dge_results))
  volcano.df <- left_join(volcano.df, bm)
  names(volcano.df)[names(volcano.df) == "external_gene_name"] <- "gene_name"
  MakeVolcanoPlot(volcano.df, plot_title = title)

  #PlotDGE(dge_results, bm, dest_dir, filename, title)
}

### Pathway analysis. Use all genes, i.e., use sce not sce_nano. ###
rm(bm)

filter <- rowSums(assay(sce) > MIN_READS) > MIN_CELLS
table(filter)
sce_ <- sce[filter,]

counts_ <- assay(sce_)
ensembl_gene_ids <- rownames(counts_)
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
                         'start_position', 'end_position', 'percentage_gene_gc_content'),
            filters = 'ensembl_gene_id',
            values = ensembl_gene_ids,
            mart = mart, uniqueRows = T)
idx <- match(ensembl_gene_ids, bm$ensembl_gene_id)
bm <- bm[idx,]
dim(bm)

not_na_idx <- !is.na(bm$entrezgene_id)
counts_ <- counts_[not_na_idx,]
rownames(counts_) <- bm$entrezgene_id[not_na_idx]
dim(counts_)
head(counts_)

dge <- DGEList(counts_, group = sce_$Clone, remove.zeros = TRUE)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + Clone, data = colData(sce_))
dge <- estimateDisp(dge, design)

lcpm_ <- edgeR::cpm(dge, log=TRUE)
dim(lcpm_)
colnames(design)

clones
contr.matrix <- makeContrasts("Clone1 vs Clone2" = Clone2 - Clone1,
                              "Clone2 vs Clone3" = Clone3 - Clone2,
                              "Clone3 vs Clone4" = Clone4 - Clone3,
                              "Clone4 vs Clone5" = Clone5 - Clone4,
                              "Clone5 vs Clone6" = Clone6 - Clone5,
                              "Clone6 vs Clone7" = Clone7 - Clone6,
                              "Clone7 vs Clone8" = Clone8 - Clone7,
                              levels = colnames(design))

v <- voom(dge, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)

# Perform camera test on Hallmark.
#load("human_c5_v5p2.rdata")
#indices <- ids2indices(Hs.c5, rownames(dge))

# Perform camera test on Hallmark set.
load("../data/human_H_v5p2.rdata")
indices <- ids2indices(Hs.H, rownames(dge))

df <- data.frame()

for (i in 1:dim(contr.matrix)[2]) {
    test_ret <- camera(v, indices, contrast=contr.matrix[,i])
    cbind(rownames(test_ret), test_ret$FDR, test_ret$Direction)
    df <- rbind(df, data.frame(Comparison=colnames(contr.matrix)[i], Hallmark=rownames(test_ret), NegLogPvalue=-log2(test_ret$PValue), Direction=test_ret$Direction, FDR=test_ret$FDR))
}

df$PrettyName <- as.character(df$Hallmark)
for (i in 1:length(df$Hallmark)) {
  df$PrettyName[i] <- FormatHallmarkName(as.character(df$Hallmark[i]))
}

dim(df)
ret <- df %>% group_by(Hallmark) %>% summarise(n_sig = sum(FDR < FDR_THRESHOLD))
hallmarks <- subset(ret, n_sig > 0)$Hallmark
df_ <- subset(df, Hallmark %in% hallmarks)
head(df_)

hallmarks <- unlist(strsplit(as.character(df$Hallmark), split = "HALLMARK_"))
hallmarks <- hallmarks[seq(2, length(hallmarks), 2)]
df$Hallmark <- gsub(pattern = "_", replacement = " ", x = hallmarks)

ret <- df %>% group_by(Hallmark) %>% summarise(n_sig = sum(FDR < FDR_THRESHOLD))
#df$Hallmark <- factor(df$Hallmark, levels = ret$Hallmark[order(ret$n_sig)])

sign <- 0
sign[df$Direction == "Up"] <- 1
sign[df$Direction == "Down"] <- -1
df$SignedNegLogPvalue <- df$NegLogPvalue * sign
hallmarks <- unique(df$Hallmark)
df$Hallmark <- factor(df$Hallmark, levels = df$Hallmark[order(hallmarks, decreasing = F)])


p <- ggplot(df, aes(Comparison, PrettyName, fill = SignedNegLogPvalue)) + geom_tile(col="white") + theme_bw() + geom_point(data = subset(df, FDR < FDR_THRESHOLD), aes(shape=Direction))
p <- p + scale_shape_manual(values = c(6, 2))
p <- p + theme_bw() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
p <- p + ylab("") + xlab("") + labs(fill="(Signed)\n-Log p-value") + labs(shape=paste("FDR <", FDR_THRESHOLD))
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text.x = element_text(size = base_size * 1.3, angle = 300, vjust = 0.1, hjust = 0, colour = "black"))
p <- p + theme(axis.text.y = element_text(size = base_size * 1.5))
p <- p + theme(legend.text=element_text(size=base_size * 1.5), legend.title = element_text(size=base_size * 1.5))
p
ggsave(paste(dest_dir, "Nano_Hallmark.pdf", sep="/"), p, height = 13.5, width=7.5, units = "in")

# Let's make boxplots on clones at pathway level.
rm(bm)
# Filter genes.
counts_ <- assay(sce_)
dim(counts_)

ensembl_gene_ids <- rownames(counts_)
bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
                         'start_position', 'end_position', 'percentage_gene_gc_content'),
            filters = 'ensembl_gene_id',
            values = ensembl_gene_ids,
            mart = mart, uniqueRows = T)
idx <- match(ensembl_gene_ids, bm$ensembl_gene_id)
bm <- bm[idx,]

not_na_idx <- !is.na(bm$entrezgene_id)
counts_ <- counts_[not_na_idx,]
rownames(counts_) <- bm$entrezgene_id[not_na_idx]
dim(counts_)
head(counts_)

cell.df_ <- subset(cell.df, Cell %in% colnames(counts_))
dim(cell.df_)

# counts: matrix of counts for each cell with rownames given by EntrezIDs.
PathwayAverageExpression <- function(cell.df, pathway_entrez_ids, counts, avg = TRUE) {
  gene_idxs <- rownames(counts_) %in% pathway_entrez_ids
  temp <- as.data.frame(matrix(counts_[gene_idxs,], ncol = dim(counts)[2], nrow = sum(gene_idxs)))
  names(temp) <- colnames(counts)
  if (avg) {
    pathway_expr <- colMeans(temp)
  } else {
    pathway_expr <- colSums(temp)
  }
  ret <- melt(pathway_expr)
  ret.df <- data.frame(SampleName=rownames(ret), PathwayExpr=ret$value)
  ret.df <- left_join(ret.df, cell.df)
  # Sort ret.df by CellID.
  ret.df <- ret.df[order(ret.df$Cell),]
  return(ret.df)
}

PlotPathwayExprAverages <- function(pathway_name, df, log_scale=F, base_size = 12) {
  if (log_scale) {
    #p <- ggplot(df, aes(Clone, log1p(PathwayExpr), fill = Clone))
    p <- ggplot(df, aes(Clone, log1p(PathwayExpr)))
  } else {
    #p <- ggplot(df, aes(Clone, PathwayExpr, fill = Clone))
    p <- ggplot(df, aes(Clone, PathwayExpr))
  }
  p <- p + geom_boxplot(fill='gray') + theme_bw()
  #p <- p + scale_fill_brewer(type = "qual", palette = "Set1")
  p <- p + theme(legend.position = "none") + xlab("Clone")
  if (log_scale) {
    p <- p + ylab("Log pathway expression")
  } else {
    p <- p + ylab("Pathway expression")
  }
  p <- p + ggtitle(FormatHallmarkName(pathway_name)) + theme(plot.title = element_text(size = base_size * 2))
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
  p <- p + theme(axis.text.x = element_text(size = base_size * 2, angle = 270), axis.text.y = element_text(size = base_size * 2))
  return(p)
}

pathway_by_clone <- data.frame()
pathway_by_cell_list <- list()
pathway_names <- c()
pathway_by_cell <- matrix(0, nrow = dim(cell.df_)[1], ncol = length(Hs.H))
for (i in 1:length(Hs.H)) {
  pathway_entrez_ids <- Hs.H[[i]]
  pathway_name <- names(Hs.H)[i]
  pathway_names <- c(pathway_names, pathway_name)
  ret <- PathwayAverageExpression(cell.df_, pathway_entrez_ids, counts_, avg=F)

  p <- PlotPathwayExprAverages(pathway_name, ret, log_scale = T)
  ggsave(filename = paste("../_figures/HER2_pos/Nano_", pathway_name, "_boxplot.pdf", sep=""),
         p, height = 8, units = "in")

  ret2 <- ret %>% group_by(Clone) %>% summarise(avg = mean(PathwayExpr), std = sd(PathwayExpr))
  pathway_by_clone <- rbind(pathway_by_clone, data.frame(Clone=ret2$Clone, MeanExpr=ret2$avg, Pathway=pathway_name))

  pathway_by_cell[,i] <- ret$PathwayExpr
  pathway_by_cell_list[[i]] <- data.frame(SampleName=ret$SampleName,
                                          Clone=ret$Clone,
                                          Region=ret$Region,
                                          Pathway=pathway_name,
                                          PathwayExpr=ret$PathwayExpr)
}

names(pathway_by_clone)
unique(pathway_by_clone$Pathway)

# t-SNE on clones at pathway level
set.seed(1)
ret <- dcast(pathway_by_clone, formula = "Clone ~ Pathway", value.var = "MeanExpr")
tsne_ret <- Rtsne::Rtsne(log1p(ret[,-1]), k = 2, perplexity=0.5, max_iter=5000)
df <- data.frame(Dim1=tsne_ret$Y[,1], Dim2=tsne_ret$Y[,2], Clone=ret$Clone)
p <- ggplot(df, aes(Dim1, Dim2, colour=Clone)) + geom_point() + theme_classic()
p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2), legend.text = element_text(size=base_size))
p # Not very intersting? Separates clones 1-3 form 4-8.

rownames(pathway_by_cell) <- pathway_by_cell_list[[1]]$SampleName
colnames(pathway_by_cell) <- pathway_names
pathway_sce <- SingleCellExperiment(assays = list(counts = t(as.matrix(pathway_by_cell))),
                                    rowData = data.frame(Pathway=pathway_names),
                                    colData = data.frame(SampleName = rownames(pathway_by_cell),
                                                         Clone=pathway_by_cell_list[[1]]$Clone,
                                                         Region=pathway_by_cell_list[[1]]$Region))
logcounts(pathway_sce) <- log1p(counts(pathway_sce))
pathway_sce_umap <- runUMAP(pathway_sce, name="UMAP")
W <- reducedDim(pathway_sce_umap, "UMAP")
df <- data.frame(Dim1=W[,1], Dim2=W[,2], Clone=colData(pathway_sce)$Clone, colData(pathway_sce)$Region)
p <- ggplot(df, aes(Dim1, Dim2, col = Clone)) + geom_point() + theme_classic()
p <- p + scale_color_brewer(type = "qual", palette = "Set1")
p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2), legend.text = element_text(size=base_size))
p # Nothing interesting.


# Plot of U-Map dimension but using color coding by hallmarks.
# For each cell, for each pathway, accumulate the gene expression
set.seed(1234551)
logcounts(sce_) <- log1p(counts(sce_))
sce_umap_ <- runUMAP(sce_, name="UMAP")
W <- reducedDim(sce_umap_, "UMAP")
df <- data.frame(Dim1=W[,1], Dim2=W[,2], Clone=colData(sce_umap_)$Clone)
p <- ggplot(df, aes(Dim1, Dim2, colour=Clone)) + geom_point(size=base_size*0.2) + theme_classic()
p

#col_data <- colData(sce_umap_)

# DO NOT RUN: not very interesting plots.
#clone_shapes <- c(15:18)
# clone_shapes <- c(4, 16, 17, 18)
# for (i in 1:length(Hs.H)) {
#   pathway_entrez_ids <- Hs.H[[i]]
#   pathway_name <- names(Hs.H)[i]
#   ret <- PathwayAverageExpression(cell.df_, pathway_entrez_ids, counts_, avg=F)
#
#   df <- data.frame(Dim1=W[,1], Dim2=W[,2], PathwayExpr = ret$PathwayExpr, Clone=as.numeric(colData(sce_umap_)$Clone))
#   df$CloneShape <- df$Clone
#   df$CloneShape[df$Clone == 1] <- "1"
#   df$CloneShape[df$Clone == 2] <- "2"
#   df$CloneShape[df$Clone == 3] <- "3"
#   df$CloneShape[df$Clone >= 4] <- "4-8"
#   df$CloneShape <- factor(df$CloneShape, levels = c("1", "2", "3", "4-8"))
#
#   p <- ggplot(df, aes(Dim1, Dim2, colour=log1p(PathwayExpr), shape=CloneShape)) + geom_point(size=base_size*0.2) + theme_classic()
#   p <- p + scale_color_gradient2(
#     low = "blue",
#     mid = "yellow",
#     high = "red",
#     space = "Lab",
#     midpoint = max(log1p(df$PathwayExpr))/2,
#     na.value = "white",
#     guide = "colourbar",
#     aesthetics = "colour")
#   p <- p + scale_shape_manual(values=clone_shapes)
#   p <- p + labs(colour = "Log1p(expr)", shape = "Clone")
#   p <- p + theme(axis.title.x = element_text(size = base_size*2)) + theme(axis.title.y = element_text(size = base_size*2))
#   p <- p + theme(axis.text.x = element_text(size = base_size*1.5), axis.text.y = element_text(size = base_size*1.5))
#   p <- p + ggtitle(FormatHallmarkName(pathway_name)) + theme(plot.title = element_text(size = base_size * 2))
#   p <- p + theme(legend.text=element_text(size=base_size * 1.5), legend.title = element_text(size=base_size * 1.5))
#   p <- p + guides(colour = guide_colorbar(order = 0), shape = guide_legend(order = 1))
#   ggsave(filename = paste(dest_dir, "/Nano_", pathway_name, "_UMAP.pdf", sep=""), p)
# }


# Generate mutation plot.
# For each clone, generate a plot of mutation status for genes in the NanoString list.
# First, take intersection of SNVs used in the analysis with the genes in the NanoString list.
# Plot the mutation status for each clone.

exons <- read.table("/Users/seonghwanjun/data/references/exons.bed", header = F, as.is = T)
names(exons) <- c("CHR", "START", "END", "GENE")
exons.gr <- ConstructGranges(exons$CHR, start = exons$START, width = exons$END - exons$START)
dat.gr <- ConstructGranges(dat$CHR, dat$POS, width = 0)
ret <- findOverlaps(dat.gr, exons.gr)

df <- dat[ret@from,]
exons.df <- exons[ret@to,]
df$GENE <- exons.df$GENE

library(stringr)
all_symbols <- unlist(lapply(nano_string$Alias...Prev.Symbol, function(row) {
  ret <- strsplit(as.character(row), split = ",")[[1]]
  str_trim(ret)
}))
all_symbols <- c(as.character(nano_string$Official.Symbol), all_symbols)
all_symbols <- unique(all_symbols)

df.nano_pre <- subset(df, GENE %in% all_symbols)

# Remove rows that have duplicate genes and duplicate loc.
loc_gene <- paste(df.nano_pre$ID, df.nano_pre$GENE, sep=":")
df.nano_pre <- df.nano_pre[!duplicated(loc_gene),]

# Check for mutation status for each clone.
df.nano <- left_join(df.nano_pre, datum2node, by = "ID")
df.nano <- df.nano[order(nchar(df.nano$Node), df.nano$Node),]
df.nano[,c("Node", "Clone")]
head(df.nano)

# Do it again, use original datum2node
datum2node_ <- read.table(paste(rep_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(datum2node_) <- c("ID", "Node")
table(datum2node_$Node)
df.nano_ <- left_join(df.nano_pre, datum2node_, by = "ID")
df.nano_ <- df.nano_[order(nchar(df.nano_$Node), df.nano_$Node),]
head(df.nano_)
table(df.nano_$Node) # Nothing assigned to clones 0_0_1 or 0_0_0_0_1 -- it seems safe to ignore/merge those mutations.


df.nano$GENE <- factor(df.nano$GENE, levels = df.nano$GENE)
df.nano$MUT <- 1

df.nano.melted <- dcast(df.nano, formula = "GENE ~ Clone", value.var = "MUT")

# Can't think of any smarter way to do this... come back and make it better later.
for (i in 3:dim(df.nano.melted)[2]) {
  df.nano.melted[1:(which(df.nano.melted[,i] == 1)[1] - 1),i] <- 1
}
df.nano.melted[is.na(df.nano.melted)] <- 0
df.nano <- melt(df.nano.melted)
names(df.nano) <- c("GENE", "Clone", "Mutation")
#df.nano$Clone <- factor(df.nano$Clone, levels = rev(unique(df.nano$Clone)))

df.nano$Mutation[df.nano$Mutation == 0] <- "Absent"
df.nano$Mutation[df.nano$Mutation == 1] <- "Present"
cols <- c("grey50", "dark blue")
p <- ggplot(df.nano, aes(Clone, GENE)) + geom_tile(aes(fill = Mutation), col = "white")
p <- p + scale_fill_manual(values = cols)
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme_bw()
p <- p + theme(axis.ticks = element_blank()) + ylab("GENE")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text.x = element_text(size = base_size*1.5), axis.text.y = element_text(size = base_size*1.5))
p <- p + theme(axis.title.x = element_text(size=base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text=element_text(size=base_size * 2), legend.title = element_text(size=base_size * 2))
p <- p + theme(legend.title = element_blank(), legend.position = "top")
p
ggsave(filename = paste(dest_dir, "/Mutation_by_Clone.pdf", sep=""), p, width = 5.5, height = 11, units = "in")

#### Plot inferred cellular prevalences for each region using a barplot. ###
cell_prev <- read.table(paste(rep_path, "/joint/tree0/cellular_prev.tsv", sep=""), header=F, sep="\t", as.is = T)
cell_prevs <- matrix(as.numeric(unlist(lapply(cell_prev[,2], function(row) {
  strsplit(row, split = ",")[[1]]
}))), ncol = 5, byrow = T)
cell_prevs.df <- data.frame(ID = cell_prev[,1], CellularPrevalence = cell_prevs)
names(cell_prevs.df) <- c("ID", REGIONS)

unmodified_datum2node <- read.table(paste(rep_path, "/joint/tree0/datum2node.tsv", sep=""), header=F, sep="\t", as.is = T)
names(unmodified_datum2node) <- c("ID", "Node")

lineage_datum2node <- subset(unmodified_datum2node, Node %in% datum2node$Node)

lineage_datum2node$Clone <- lineage_datum2node$Node
nodes <- unique(lineage_datum2node$Node)
nodes_ordered <- nodes[order(nchar(nodes), nodes)]
clone_count <- length(nodes_ordered)
for (i in 1:clone_count) {
  idx <- lineage_datum2node$Node == nodes_ordered[i]
  lineage_datum2node$Clone[idx] <- i
}

# Remove SNVs that aren't assigned to the lineage.
cell_prevs.df <- subset(cell_prevs.df, ID %in% lineage_datum2node$ID)
cell_prevs.df <- left_join(cell_prevs.df, lineage_datum2node)
cell_prev.clone <- cell_prevs.df[!duplicated(cell_prevs.df$Clone),]

cell_prev.clone.melted <- melt(cell_prev.clone)
cell_prev.clone.melted$Clone <- as.numeric(cell_prev.clone.melted$Clone)
names(cell_prev.clone.melted) <- c("ID", "Node", "Clone", "Region", "CellPrev")
p <- ggplot(cell_prev.clone.melted, aes(x = Clone, y = CellPrev, col = Region)) + geom_line() + geom_point()
p <- p + theme_classic() + ylab("Cellular prevalences")
p <- p + scale_x_continuous(breaks = 1:clone_count)
p <- p + theme(axis.text.x = element_text(size = base_size*1.5), axis.text.y = element_text(size = base_size*1.5))
p <- p + theme(axis.title.x = element_text(size=base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text=element_text(size=base_size * 1.5), legend.title = element_text(size=base_size * 1.5))
p <- p + theme(legend.position = c(0.9, 0.8))
p
ggsave(p, filename = paste(dest_dir, "CellPrev.pdf", sep="/"))

# We will make a plot for clone fraction.
cell_prev.clone.melted$CellFrac <- cell_prev.clone.melted$CellPrev
# Let's order by clone ID then region.
cell_prev.clone.melted <- cell_prev.clone.melted[order(cell_prev.clone.melted$Clone, cell_prev.clone.melted$Region),]
for (clone in 1:(clone_count-1)) {
  parent_clone_idx <- (cell_prev.clone.melted$Clone == clone)
  child_clone_idx <- (cell_prev.clone.melted$Clone == (clone + 1))
  cell_prev.clone.melted[parent_clone_idx,"CellFrac"] <- cell_prev.clone.melted[parent_clone_idx,"CellPrev"] - cell_prev.clone.melted[child_clone_idx,"CellPrev"]
}

region_colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")
p <- ggplot(cell_prev.clone.melted, aes(x = Clone, y = CellFrac, fill=Region)) + geom_bar(stat="identity", position="dodge")
p <- p + theme_bw() + ylab("Clone fraction")
p <- p + scale_x_continuous(breaks = 1:clone_count)
p <- p + theme(axis.text.x = element_text(size = base_size*1.5), axis.text.y = element_text(size = base_size*1.5))
p <- p + theme(axis.title.x = element_text(size=base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text=element_text(size=base_size * 1.5), legend.title = element_text(size=base_size * 1.5))
p <- p + theme(legend.position = c(0.9, 0.8))
p <- p + scale_fill_manual(values = region_colors)
p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6), breaks = seq(0, 1, 0.1))
#p <- p + ggtitle("Cell fraction by clone and region")
p
ggsave(p, filename = paste(dest_dir, "CellFrac.pdf", sep="/"))

# Plot with Region on the x-axis and use color for the clones.
clone_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#0072B2", "#D55E00")

p <- ggplot(cell_prev.clone.melted, aes(x = Region, y = CellFrac, fill=as.factor(Clone))) + geom_bar(stat="identity", position="dodge")
p <- p + theme_bw() + ylab("Clone fraction")
p <- p + theme(axis.text.x = element_text(size = base_size*1.5), axis.text.y = element_text(size = base_size*1.5))
p <- p + theme(axis.title.x = element_text(size=base_size * 2), axis.title.y=element_text(size = base_size * 2))
p <- p + theme(legend.text=element_text(size=base_size * 1.5), legend.title = element_text(size=base_size * 1.5))
p <- p + scale_fill_manual(values = clone_colors)
p <- p + scale_y_continuous(expand = c(0,0), limits = c(0, 0.6), breaks = seq(0, 1, 0.1))
p <- p + labs(fill="Clone")
p
ggsave(p, filename = paste(dest_dir, "CellFrac.pdf", sep="/"))
