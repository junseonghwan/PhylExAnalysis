
# Download the cherry data from Zenodo and unzip under data/simulation
library(ggplot2)
library(dplyr)
library(rjson)
library(PhylExR)
library(Rcpp)
# TODO: This code relies on older version of the code -- update required.
sourceCpp(file = "PrePhylExR/src/rcpp_hello_world.cpp")

# dropout_hp: named list containing 'alpha' and 'beta'.
# bursty_hp: named list containing 'alpha' and 'beta'.
# biallelic_hp: numeric matrix with the first column the alpha and the second column the beta
AssignCells <- function(cell_data,
                        datum2node,
                        dropout_hp,
                        bursty_hp,
                        biallelic_hp) {
  names(datum2node) <- c("ID", "Node")
  snv_count <- length(datum2node$ID)
  
  nodes <- as.character(unique(datum2node$Node))
  nodes_breadth_first <- nodes[order(nchar(nodes), nodes)]
  
  mut_ids <- unique(as.character(cell_data$ID))
  mut_ids <- mut_ids[order(nchar(mut_ids), mut_ids)]
  cells <- unique(as.character(cell_data$Cell))
  cells <- cells[order(nchar(cells), cells)]
  
  cell_data$ID <- factor(cell_data$ID, levels = mut_ids)
  cell_data$Cell <- factor(cell_data$Cell, levels = cells)
  if (!("b" %in% names(cell_data))) {
    cell_data$b <- cell_data$d - cell_data$a
  }
  var_reads <- reshape2::dcast(cell_data, Cell ~ ID, value.var = "b")
  var_reads[is.na(var_reads)] <- 0
  total_reads <- reshape2::dcast(cell_data, Cell ~ ID, value.var = "d")
  total_reads[is.na(total_reads)] <- 0
  
  dropout_hp_mat <- matrix(c(dropout_hp$alpha, dropout_hp$beta), nrow = snv_count, ncol=2, byrow=T)
  bursty_hp_mat <- matrix(c(bursty_hp$alpha, bursty_hp$beta), nrow = snv_count, ncol=2, byrow=T)
  
  log_unnorm_liks <- IdentifyCellMutationStatus(datum2node,
                                                nodes_breadth_first,
                                                mut_ids,
                                                as.matrix(var_reads[,-1]),
                                                as.matrix(total_reads[,-1]),
                                                as.matrix(dropout_hp_mat),
                                                as.matrix(bursty_hp_mat),
                                                as.matrix(biallelic_hp))
  cells_with_no_reads <- (rowSums(log_unnorm_liks) == 0)
  log_unnorm_liks2 <- log_unnorm_liks[!cells_with_no_reads,]
  #dim(log_unnorm_liks2)
  log_norms <- apply(log_unnorm_liks2, 1, logSumExp)
  cell_assign_probs <- exp(log_unnorm_liks2 - log_norms)
  err <- (rowSums(cell_assign_probs) - 1) > 1e-6
  if (sum(err) > 0) {
    print(paste("Row does not sum to 1 for ", which(err)))
  }
  cell_assignment <- nodes_breadth_first[apply(cell_assign_probs, 1, which.max)]
  return(data.frame(Cell = var_reads[,1], Node = cell_assignment))
}

dat <- read.table("data/simulation/cherry/genotype_ssm.txt", header=T, sep="\t")
sc <- read.table("data/simulation/cherry/simul_sc.txt", header=T, sep="\t")
sc_hp <- read.table("data/simulation/cherry/simul_sc_hp.txt", header=T, sep="\t")

snv_count <- dim(dat)[1]

# Plot cell mutation profile.
sc$b <- sc$d - sc$a

cells <- unique(sc$Cell)
ids <- dat$ID
sc$Cell <- factor(sc$Cell, levels = cells)
sc$ID <- factor(sc$ID, levels = ids)

base_size <- 11
p <- ggplot(sc, aes(ID, Cell, fill = b)) + geom_tile(colour = "white")
p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
p <- p + ylab("Cells") + xlab("Loci")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(legend.position = "none", axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("_figures/simulation/cherry/cherry_sc.pdf", p)

cluster_labels <- read.table("data/simulation/cherry/genotype/joint/tree0/cluster_labels.tsv", header=F, sep="\t")
names(cluster_labels) <- c("ID", "Cluster")
cluster_labels$Cluster <- as.character(cluster_labels$Cluster)

datum2node <- read.table("data/simulation/cherry/genotype/joint/tree0/datum2node.tsv", header=F, sep="\t")

# Assign cell to nodes.
dropout_hp <- list(alpha=0.01, beta=1)
bursty_hp <- list(alpha=1, beta=0.01)
#cell_assignment <- AssignCellsBursty(sc, datum2node, bursty_hp, sc_hp)
cell_assignment <- AssignCells(sc, datum2node, dropout_hp, bursty_hp, sc_hp[,2:3])
cell_assignment_truth <- read.table("data/simulation/cherry/cell2node.txt", header=F, sep="\t")
names(cell_assignment_truth) <- c("Cell", "Node")

cell_assignment$Node <- as.character(cell_assignment$Node)

ids <- cluster_labels[order(cluster_labels$Cluster),"ID"]
cells <- cell_assignment[order(nchar(cell_assignment$Node), cell_assignment$Node),"Cell"]
sc$Cell <- factor(sc$Cell, levels = cells)
sc$ID <- factor(sc$ID, levels = ids)

p <- ggplot(sc, aes(ID, Cell, fill = b)) + geom_tile(colour = "white")
p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
p <- p + ylab("Cells") + xlab("Loci")
p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
p <- p + theme(legend.position = "none", axis.ticks = element_blank())
p <- p + theme(axis.title.x =element_text(size = base_size * 2))
p <- p + theme(axis.title.y =element_text(size = base_size * 2))
p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
#ggsave("_figures/simulation/cherry/cherry_inferred.pdf", p, width = 6, height=6, units="in")
ggsave("_figures/simulation/cherry/cherry_inferred.pdf", p)

### Save the source data for Figure 1e-f ###
head(sc)
sc_ <- left_join(sc, cell_assignment)
sc_ <- left_join(sc_, cluster_labels)
names(sc_) <- c("ID", "Cell", "a", "d", "SampleName", "b", "CellCluster", "SNVCluster")
write.csv(x = sc_[,-c(5)], file = "data/NatComm/Figure1e-f.csv", quote = F, row.names = F)

