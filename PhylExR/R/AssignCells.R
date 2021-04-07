# Assign cells to nodes.
# 1. Read in cell reads and mutations harboured by each node.
# 2. For each cell, compute the likelihood that it should be assigned to a node.
# 3. Compute the marginalization to obtain the cell assignment probabilities.
# 4. Choose the one with maximal probability.
ScLikelihood <- function(has_snv,
                         var_read_count,
                         total_read_count,
                         biallelic_alpha,
                         biallelic_beta,
                         bursty_alpha,
                         bursty_beta,
                         seq_err = 0.001,
                         beta_binom_mixture_prob = 0.5) {
  log_lik <- 0
  
  # if sc has mutation s, then there are 2 cases
  # 1. Bi-allelic,
  # 2. Bursty distribution.
  if (total_read_count == 0) {
    return(0)
  }
  if (has_snv) {
    log_lik_biallelic <- TailRank::dbb(var_read_count, total_read_count, biallelic_alpha, biallelic_beta, log = TRUE)
    log_lik_biallelic <- log_lik_biallelic + log(1-beta_binom_mixture_prob)
    log_lik_bursty <- TailRank::dbb(var_read_count,  total_read_count, bursty_alpha, bursty_beta, log = TRUE)
    log_lik_bursty <- log_lik_bursty + log(beta_binom_mixture_prob)
    log_lik <- LogAdd(log_lik_biallelic, log_lik_bursty)
  } else {
    log_lik <- TailRank::dbb(var_read_count, total_read_count, seq_err, 1 - seq_err, log = TRUE)
  }
  
  return(log_lik);
}

# cell_data: single cell input data for PhylEx.
# bursty_hp: named list containing 'alpha' and 'beta'.
# biallelic_hp: numeric matrix with the first column the alpha and the second column the beta
ComputeCellAssignmentProbability <- function(cell_data,
                                              datum2node,
                                              bursty_hp,
                                              biallelic_hp,
                                              include_normal_clone = TRUE,
                                             epsilon=0.01) {
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
  
  node2snvs <- IdentifyNodeMutationStatus(datum2node,
                                          nodes_breadth_first,
                                          mut_ids)

  # node2snvs is a matrix of dimension node_count x snv_count.
  # Entry (i,j) equal to 1 indicates that node i has snv j.
  # For each cell, we will compute the probability of assigning to each node.
  node2snvs.df <- as.data.frame(node2snvs)
  rownames(node2snvs.df) <- nodes_breadth_first
  colnames(node2snvs.df) <- mut_ids

  cells <- unique(cell_data$Cell)
  cell_count <- length(cells)
  node_count <- length(nodes_breadth_first)
  cell_data$b <- cell_data$d - cell_data$a
  log_prob_matrix <- matrix(0, nrow = cell_count, ncol = node_count)
  for (i in 1:length(cells)) {
    cell <- cells[i]
    for (node in 1:length(nodes_breadth_first)) {
      sc_cell <- subset(cell_data, Cell == cell)
      node_mut_status <- mut_ids[which(node2snvs.df[node,] == 1)]
      has_snv <- sc_cell$ID %in% node_mut_status
      mut_ids_to_compute <- biallelic_hp$ID %in% sc_cell$ID
      biallelic_alpha <- biallelic_hp[mut_ids_to_compute,"alpha"]
      biallelic_beta <- biallelic_hp[mut_ids_to_compute,"beta"]
      df <- data.frame(has_snv=has_snv, b=sc_cell$b, d=sc_cell$d, alpha=biallelic_alpha, beta=biallelic_beta)
      log_liks <- apply(df, 1, function(row) {
        row <- as.data.frame(t(row))
        ScLikelihood(row$has_snv, row$b, row$d, row$alpha, row$beta, bursty_hp$alpha, bursty_hp$beta, seq_err = epsilon)
      })
      log_prob_matrix[i,node] <- sum(log_liks)
    }
  }
  
  cells_with_no_reads <- (rowSums(log_prob_matrix) == 0)
  log_unnorm_liks2 <- log_prob_matrix[!cells_with_no_reads,]
  log_unnorm_liks2 <- matrix(log_unnorm_liks2, ncol = length(nodes_breadth_first))
  #dim(log_unnorm_liks2)
  log_norms <- apply(log_unnorm_liks2, 1, logSumExp)
  cell_assign_probs <- exp(log_unnorm_liks2 - log_norms)
  err <- (rowSums(cell_assign_probs) - 1) > 1e-6
  
  if (sum(err) > 0) {
    print(paste("Row does not sum to 1 for ", which(err)))
  }
  colnames(cell_assign_probs) <- nodes_breadth_first
  rownames(cell_assign_probs) <- cells
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

