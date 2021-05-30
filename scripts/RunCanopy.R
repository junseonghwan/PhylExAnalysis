args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
DATA_PATH <- args[2]
K_BEGIN <- as.numeric(args[3])
K_END <- as.numeric(args[4])

#SEED <- 157
#DATA_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/"
#K_BEGIN <- 3
#K_END <- 12

SNV_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
OUTPUT_PATH <- paste(DATA_PATH, "canopy", sep="/")
FALCON_CNV_PATH <- paste(DATA_PATH, "falcon.txt", sep="/")

library(Canopy)
library(cardelino)
library(GenomicRanges)
library(PhylExR)

set.seed(SEED)

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH)
}
setwd(OUTPUT_PATH)

bulk <- read.table(SNV_PATH, header=T, sep="\t")
X <- as.matrix(bulk$d - bulk$b)
R <- as.matrix(bulk$b)
rownames(X) <- bulk$ID
rownames(R) <- bulk$ID
colnames(X) <- "OV2295"
colnames(R) <- "OV2295"

K <- K_BEGIN:K_END
numchains <- 4
min_iter <- 10000
max_iter <- 100000
projname <- "canopy"

snv_count <- dim(bulk)[1]

if (file.exists(FALCON_CNV_PATH)) {
    falcon <- read.table(FALCON_CNV_PATH, sep="\t", header=T, as.is = TRUE)
    cna_count <- dim(falcon)[1]

    temp <- as.matrix(falcon[,-(1:3)])
    rownames(temp) <- paste("c", 0:(cna_count-1), sep="")

    WM <- as.matrix(temp[,"Major_copy"])
    Wm <- as.matrix(temp[,"Minor_copy"])
    epsilonM <- as.matrix(temp[,"Major.sd"])
    epsilonm <- as.matrix(temp[,"Minor.sd"])

    falcon.gr <- ConstructGranges(falcon$CHR, falcon$st_bp, width = falcon$end_bp - falcon$st_bp)
    bulk.gr <- ConstructGranges(bulk$CHR, bulk$POS, width = 0)
    ret <- findOverlaps(bulk.gr, falcon.gr)

    Y <- matrix(0, nrow = snv_count, ncol = cna_count + 1)
    for (i in 1:length(ret)) {
        Y[ret@from[i],ret@to[i]+1] <- 1
    }
    # Check for any SNV that doesn't fall on any CNA.
    non_cna_snv_idxs <- which(!(1:snv_count %in% ret@from))
    Y[non_cna_snv_idxs,1] <- 1
    # Check that every row has exactly one entry with value of 1.
    print(paste("rowSums(Y) has exactly one entry with value of 1: ", mean(rowSums(Y) == 1) == 1))

    sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM,
                              epsilonm = epsilonm, Y = Y, K = K, numchain = numchains,
                              max.simrun = max_iter, min.simrun = min_iter,
                              writeskip = 200, projectname = projname, cell.line = TRUE,
                              plot.likelihood = TRUE)
} else {
    sampchain <- canopy.sample.nocna(R, X, K, numchains, projectname = projname,
                                     max.simrun = max_iter, min.simrun = min_iter,
                                     writeskip = 200, cell.line=TRUE, plot.likelihood=T)
}

burnin <- 10
thin <- 5
bic <- canopy.BIC(sampchain = sampchain, projectname = projname, K = K,
                  numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK <- K[which.max(bic)]
post <- canopy.post(sampchain = sampchain, projectname = projname, K = K,
                    numchain = numchain, burnin = burnin, thin = thin,
                    optK = optK, post.config.cutoff = 0.01)

samptreethin <- post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik <- post[[2]]   # likelihoods of trees in samptree
config <- post[[3]]
config.summary <- post[[4]]
print(config.summary)

config.i <- config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree <- canopy.output(post, config.i, C=NULL)
pdf.name <- paste(OUTPUT_PATH, 'highest_likelihood.pdf', sep='/')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), Cluster=clone_name)
write.table(predicted, paste(OUTPUT_PATH, "predicted.csv", sep="/"), sep=",", col.names=T, row.names=F, quote=F)

snv_count <- dim(output.tree$Z)[1]
A <- matrix(0, snv_count, snv_count)
for (i in 1:snv_count) {
    clone_i <- sort(which(output.tree$Z[i,] == 1))
    for (j in 1:snv_count) {
        clone_j <- sort(which(output.tree$Z[j,] == 1))
        if (!setequal(clone_i, clone_j)) {
            A[i,j] <- as.numeric(all(clone_j %in% clone_i))
        }
    }
}
A[2,]
output.tree$sna[1:4,]
output.tree$Z[1:4,]
write.table(A, paste(OUTPUT_PATH, "ancestral_matrix.csv", sep="/"), row.names = F, quote=F, col.names = F)

save.image(file = paste(OUTPUT_PATH, "/canopy_postmcmc_image.rda", sep=""), compress = "xz")

# Now let's run Cardelino.
# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), CloneID=clone_name)
write.table(predicted, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/Canopy/predicted.csv", sep=",", col.names=T, row.names=F, quote=F)

# Now, we are ready to run Cardelino.
# First, get the reads from single cell data.
sc_reads <- read.table("/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/sc.txt", header=T)
sc_reads$b <- sc_reads$d - sc_reads$a
ids <- rownames(output.tree$Z)
sc_reads$ID <- factor(sc_reads$ID, levels = ids)
head(sc_reads)
cells <- unique(sc_reads$Cell)
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
    temp <- subset(sc_reads, Cell == cells[cell_idx])
    B[temp$ID,cell_idx] <- temp$b
    D[temp$ID,cell_idx] <- temp$d
}

cell <- "c191"
subset(sc_reads, Cell == cell)
B[B.mat[,cell] != 0,cell]
D[D.mat[,cell] != 0,cell]

# Run cardelino:
assignments <- clone_id(B, D, Config = output.tree$Z)
names(assignments)
prob_heatmap(assignments$prob)
cluster_idx <- assign_cells_to_clones(assignments$prob)
cluster_assignment <- cluster_idx$clone
table(cluster_assignment)
cell_cardelino <- data.frame(Cell = colnames(B), Node = cluster_assignment)

# Make a heatmap plot: use PhylEx tree to assign
# Order the cells by clones.
cells <- cell_cardelino$Cell
nodes <- unique(cell_cardelino$Node)
nodes[order(nchar(nodes), nodes)]
cell_ordering <- order(nchar(cell_cardelino$Node), cell_cardelino$Node)
cell_order <- cells[cell_ordering]

# Order the ids:
ids <- rownames(output.tree$Z)
id_order <- ids[order(rowSums(output.tree$Z))]

# Make a heatmap plot of the single cells and their mutation status.
sc_cardelino <- left_join(sc_reads, cell_cardelino)
head(sc_cardelino)
sc_cardelino$Cell <- factor(sc_cardelino$Cell, levels = cell_order)
sc_cardelino$ID <- factor(sc_cardelino$ID, levels = id_order)

snv_cell_coverage <- sc_cardelino %>% group_by(ID) %>% summarise(n = sum(b > 0))
snv_cell_coverage_ <- snv_cell_coverage[which(snv_cell_coverage$n > 0),]
sc_cardelino_ <- subset(sc_cardelino, ID %in% snv_cell_coverage_$ID & b >= 1)
dim(sc_cardelino)
dim(sc_cardelino_)

# Make a heatmap.
names(sc_cardelino_)
base_size <- 12
p <- ggplot(subset(sc_cardelino_, b >0), aes(ID, Cell, fill=Node)) + geom_tile(colour = "white")
p <- p + theme_bw()
p <- p + xlab("Loci") + ylab("Cells")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
#ggsave(paste(rep_path, "joint/tree0/cardelino_sc_after.pdf", sep="/") , p, width = 8.5, height = 11, units = "in")
