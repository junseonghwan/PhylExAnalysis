args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
DATA_PATH <- args[2]
OUTPUT_PATH <- args[3]
#SEED <- 1
#DATA_PATH <- "~/Dropbox/seong/PhylExNatComm/zenodo/binary_cn/"
#OUTPUT_PATH <- "~/PhylExAnalysis/_output/simul/binary_cn/"
set.seed(SEED)
if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH, recursive = TRUE)
}

library(cardelino)
library(Canopy)
library(PhylExR)

####################################################################
bursty_hp <- list(alpha=0.01, beta=0.01)

# gt: binary matrix of dimension n_cells x n_snvs with value of 1 at row c column n indicates presence of mutation n for cell c.
# probs: matrix of dimension n_cells x n_clones indicating probability of assigning a cell to a clone.
# Z: binary matrix of dimension n_snvs x n_clones indicating presence of mutation n in clone.
EvaluatePrediction <- function(gt, probs, Z)
{
  total_exp_loss <- 0
  n_cells <- dim(gt)[1]
  clones <- colnames(Z)
  for (cell in 1:n_cells)
  {
    exp_loss <- 0
    for (clone in clones)
    {
      exp_loss <- exp_loss + probs[cell,clone] * sum(abs(gt[cell,] - Z[,clone]))
    }
    total_exp_loss <- total_exp_loss + exp_loss
  }
  return(total_exp_loss/n_cells)
}
####################################################################

loss <- data.frame()
for (rep_no in 0:19) {
  REP_PATH <- paste(DATA_PATH, "/rep", rep_no, sep="")
  case_no <- 0
  CASE_PATH <- paste(REP_PATH, "/case", case_no, sep="")
  SNV_PATH <- paste(CASE_PATH, "genotype_ssm.txt", sep="/")
  SNV_GT_PATH <- paste(CASE_PATH, "datum2node.tsv", sep="/")
  CANOPY_PATH <- paste(CASE_PATH, "canopy", sep="/")
  
  ####################################################################
  # Read canopy tree.
  output.tree <- readRDS(paste(CANOPY_PATH, "canopy_output.RDS", sep="/"))
  
  # Read bulk data and ground truth assignment of bulk to clones.
  bulk <- read.table(SNV_PATH, header=T, sep="\t", as.is = T)
  snv_gt <- read.table(SNV_GT_PATH, header = F, sep="\t", as.is = T)
  
  # Get SNVs for each clone in the Canopy tree.
  nodes <- unique(snv_gt$V2)
  nodes_sorted <- nodes[order(nchar(nodes), nodes)]
  snv_list <- list()
  for (node in nodes_sorted)
  {
    parent_node <- PhylExR::get_parent_name(node)
    snv_list[[node]] <- c(snv_gt$V1[snv_gt$V2 == node], snv_list[[parent_node]])
    snv_list[[node]] <- snv_list[[node]][order(nchar(snv_list[[node]]), snv_list[[node]])]
  }
  
  ids <- rownames(output.tree$Z)
  n_snvs <- length(ids)
  
  ####################################################################
  # Evaluate mapping of scRNA-seq data on PhylEx tree vs Canopy tree.
  for (case_no in 1:3)
  {
    CASE_PATH <- paste(REP_PATH, "/case", case_no, sep="")
    SC_PATH <- paste(CASE_PATH, "simul_sc.txt", sep="/")
    SC_GT_PATH <- paste(CASE_PATH, "cell2node.txt", sep="/")
    SC_HP_PATH <- paste(CASE_PATH, "simul_sc_hp.txt", sep="/")
    
    sc <- read.table(SC_PATH, header=T, as.is = TRUE)
    sc_gt <- read.table(SC_GT_PATH, header = F, as.is = TRUE)
    sc_hp <- read.table(SC_HP_PATH, header=T)
    
    cells <- unique(sc$Cell)
    n_cells <- length(cells)
    
    cell_gt <- matrix(0, ncol = n_snvs, nrow = n_cells)
    colnames(cell_gt) <- ids
    rownames(cell_gt) <- paste("c", 0:(n_cells-1), sep="")
    for (cell in 1:n_cells)
    {
      cell_gt[cell, snv_list[[sc_gt[cell,"V2"]]]] <- 1
    }
    
    ####################################################################
    # Maps cells to Canopy tree using Cardelino.
    B <- matrix(0, nrow = n_snvs, ncol = n_cells)
    D <- matrix(0, nrow = n_snvs, ncol = n_cells)
    rownames(B) <- ids
    rownames(D) <- ids
    colnames(B) <- cells
    colnames(D) <- cells
    
    # Do a simple loop since the number of cells isn't that large.
    for (cell_idx in 1:n_cells) {
      temp <- subset(sc, Cell == cells[cell_idx])
      B[temp$ID,cell_idx] <- temp$d - temp$a
      D[temp$ID,cell_idx] <- temp$d
    }
    
    # Run Cardelino.
    #assignments <- clone_id(B, D, Config = output.tree$Z)
    assignments <- clone_id_EM(B, D, Config = output.tree$Z)
    cluster_idx <- assign_cells_to_clones(assignments$prob)
    cluster_assignment <- cluster_idx$clone
    
    ####################################################################
    # Get PhylEx results.
    best_chain <- PhylExR::FindBestRep(paste(CASE_PATH, "genotype", sep="/"), chains = 0:3)
    datum2node_path <- paste(CASE_PATH, "/genotype/chain", best_chain, "/joint/tree0/datum2node.tsv", sep="")
    datum2node <- read.table(datum2node_path, header = F)
    
    phylex_pred_prob <- ComputeCellAssignmentProbability(sc, datum2node, bursty_hp, sc_hp, include_normal_clone = T, epsilon = 0.001)
    phylex_clones <- colnames(phylex_pred_prob)
    phylex_cluster_assignment <- phylex_clones[apply(phylex_pred_prob, 1, which.max)]
    
    phylex_nodes <- unique(datum2node$V2)
    phylex_nodes_sorted <- phylex_nodes[order(nchar(phylex_nodes), phylex_nodes)]
    phylex_snv_list <- list()
    phylex_snvs_by_clone <- matrix(0, nrow = n_snvs, ncol = length(phylex_clones))
    rownames(phylex_snvs_by_clone) <- ids
    colnames(phylex_snvs_by_clone) <- phylex_clones
    for (node in phylex_nodes_sorted)
    {
      parent_node <- PhylExR::get_parent_name(node)
      while(!(parent_node %in% phylex_clones)) {
        parent_node <- PhylExR::get_parent_name(parent_node)
      }
      if (parent_node != "0") {
        phylex_snvs_by_clone[,node] <- phylex_snvs_by_clone[,parent_node]
      }
      phylex_snvs_by_clone[datum2node$V1[datum2node$V2 == node],node] <- 1
    }
    
    ####################################################################
    # Use Cardelino on PhylEx tree to show that the problem lies with clonal tree reconstruction. 
    assignments2 <- clone_id_EM(B, D, Config = phylex_snvs_by_clone)
    cluster_idx2 <- assign_cells_to_clones(assignments2$prob)
    cluster_assignment2 <- cluster_idx2$clone
    
    ####################################################################
    # Predict mutation status of cells using cells mapped on Canopy vs PhylEx clonal tree.
    cardelino_pred_err <- EvaluatePrediction(cell_gt, assignments$prob, output.tree$Z)
    cardelino_pred_err_on_phylex_tree <- EvaluatePrediction(cell_gt, assignments2$prob, phylex_snvs_by_clone)
    phylex_pred_err <- EvaluatePrediction(cell_gt, phylex_pred_prob, phylex_snvs_by_clone)
    
    loss <- rbind(loss, data.frame(values = c(cardelino_pred_err, cardelino_pred_err_on_phylex_tree, phylex_pred_err), type = c("CanopyCardelino", "PhylExCardelino", "PhylEx"), cells = n_cells, rep=rep_no))
  }
}
write.table(loss, file = paste(OUTPUT_PATH, "loss.tsv", sep="/"), col.names = T, quote = F, row.names = F, sep="\t")
