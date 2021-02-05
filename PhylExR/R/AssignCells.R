# Assign cells to nodes.
# 1. Read in cell reads and mutations harboured by each node.
# 2. For each cell, compute the likelihood that it should be assigned to a node.
# 3. Compute the marginalization to obtain the cell assignment probabilities.
# 4. Choose the one with maximal probability.


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

ComputeCellAssignmentProbability <- function(cell_data,
                                             datum2node,
                                             bursty_hp,
                                             biallelic_hp,
                                             include_normal_clone = TRUE) {
  names(datum2node) <- c("ID", "Node")
  snv_count <- length(datum2node$ID)
  
  nodes <- as.character(unique(datum2node$Node))
  nodes_breadth_first <- nodes[order(nchar(nodes), nodes)]
  if (!("0" %in% nodes_breadth_first) & include_normal_clone) {
    nodes_breadth_first <- c("0", nodes_breadth_first)
  }
  
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
  
  mutation_count <- length(mut_ids)
  bursty_hp_mat <- matrix(c(bursty_hp$alpha, bursty_hp$beta), nrow = mutation_count, ncol=2, byrow=T)
  
  # We need to pass in biallelic_hp in the same order as mut_ids.
  idxs <- which(biallelic_hp$ID %in% mut_ids)
  
  log_unnorm_liks <- IdentifyCellMutationStatusBursty(datum2node,
                                                      nodes_breadth_first,
                                                      mut_ids,
                                                      as.matrix(var_reads[,-1]),
                                                      as.matrix(total_reads[,-1]),
                                                      as.matrix(bursty_hp_mat),
                                                      as.matrix(biallelic_hp[idxs,2:3]))
  cells_with_no_reads <- (rowSums(log_unnorm_liks) == 0)
  log_unnorm_liks2 <- log_unnorm_liks[!cells_with_no_reads,]
  log_unnorm_liks2 <- matrix(log_unnorm_liks2, ncol = length(nodes_breadth_first))
  #dim(log_unnorm_liks2)
  log_norms <- apply(log_unnorm_liks2, 1, logSumExp)
  cell_assign_probs <- exp(log_unnorm_liks2 - log_norms)
  err <- (rowSums(cell_assign_probs) - 1) > 1e-6
  if (sum(err) > 0) {
    print(paste("Row does not sum to 1 for ", which(err)))
  }
  colnames(cell_assign_probs) <- nodes_breadth_first
  return(cell_assign_probs)
}

AssignCellsBursty <- function(cell_data,
                              datum2node,
                              bursty_hp,
                              biallelic_hp,
                              include_normal_clone = TRUE) {
  cell_assignment_probs <- ComputeCellAssignmentProbability(cell_data, datum2node, bursty_hp, biallelic_hp, include_normal_clone)
  nodes_breadth_first <- colnames(cell_assignment_probs)
  cell_assignment <- nodes_breadth_first[apply(cell_assignment_probs, 1, which.max)]
  return(data.frame(Cell = unique(cell_data$Cell), Node = cell_assignment))
}
