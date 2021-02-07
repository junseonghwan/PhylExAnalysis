#' @export
ConvertToEnsemblID <- function(gene_names, mart) {
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                   filters = "external_gene_name", values = gene_names,
                   mart = mart)
  return(results)
}

#' @export
ConvertToGeneName <- function(ensembl_gene_ids, mart) {
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = "ensembl_gene_id", values = ensembl_gene_ids,
                   mart = mart)
  return(results)
}

#' @export
FindBestRep <- function(chain_paths, chains) {
  log_liks <- rep(0, length(chains))
  for (chain in chains) {
    log_liks[chain + 1] <- read.table(paste(chain_paths, "/chain", chain, "/joint/tree0/log_lik.txt", sep=""), sep="", header=F)$V1
  }
  best_chain <- which.max(log_liks) - 1
  return(best_chain)
}

#' @export
EstimateHyperParameters <- function(mut_ids, sc, alpha0 = 0.05, beta0 = 0.05, alpha = 1, beta = 1, delta0 = 0.5) {
  snv_count <- length(mut_ids)
  hyper_params.df <- data.frame(ID = mut_ids, alpha = rep(1, snv_count), beta = rep(1, snv_count), delta0 = rep(0.5, snv_count))
  sc$b <- sc$d - sc$a
  temp <- subset(sc, d > 0)
  for (i in 1:snv_count) {
    id <- as.character(mut_ids[i])
    temp2 <- subset(temp, ID == id)
    if (dim(temp2)[1] > 0) {
      ww <- matrix(0, ncol = 2, nrow = dim(temp2)[1])
      ww[,1] <- TailRank::dbb(temp2$b, temp2$d, alpha0, beta0, log=TRUE)
      ww[,2] <- TailRank::dbb(temp2$b, temp2$d, alpha, beta, log=TRUE)
      log_norm <- apply(ww, 1, matrixStats::logSumExp)
      probs <- exp(ww - log_norm)
      idx <- (probs[,2] > (delta0 + 1e-16))
      if (sum(idx) > 0) {
        hyper_params.df[i,2:4] <- c(mean(temp2$b[idx]) + 1, mean(temp2$a[idx]) + 1, delta0)
      }
    }
  }
  return(hyper_params.df)
}

#' @export
SelectTopKGenes <- function(sce, topK) {
  vars <- assay(sce) %>% log1p %>% rowVars
  names(vars) <- rownames(sce)
  vars <- sort(vars, decreasing = TRUE)
  head(vars)
  sce <- sce[names(vars)[1:topK],]
  return(sce)
}

#' @export
CollapseClones <- function(datum2node, MIN_SNV_COUNT = 1) {
  # Merge any singleton cluster into its ancestor
  nodes <- as.character(unique(datum2node$Node))
  nodes <- nodes[order(nchar(nodes), nodes, decreasing = T)]
  for (node in nodes) {
    idx <- (datum2node$Node == node)
    mutation_count <- dim(datum2node[idx,])[1]
    if (mutation_count <= MIN_SNV_COUNT) {
      # Get the parent node
      new_node_name <- get_parent_name(node)
      datum2node[idx,"Node"] <- new_node_name
    }
  }
  return(datum2node)
}

#' @export
CollapseClonesByCellCount <- function(datum2node, cell.df, MIN_CELL_COUNT = 10) {
  # Merge clones with less than MIN_CELL_COUNT cells into its parent.
  nodes <- as.character(unique(datum2node$Node))
  nodes <- nodes[order(nchar(nodes), nodes, decreasing = T)]
  
  for (node in nodes) {
    cell_count <- dim(subset(cell.df, Node == node))[1]
    if (cell_count <= MIN_CELL_COUNT) {
      # Get the parent node
      new_node_name <- get_parent_name(node)
      idx <- (datum2node$Node == node)
      datum2node[idx,"Node"] <- new_node_name
      
      idx <- (cell.df$Node == node)
      cell.df[idx,"Node"] <- new_node_name
    }
  }
  return(list("datum2node"=datum2node, "cell.df"=cell.df))
}

# Hardcode pathway names.
#' @export
FormatHallmarkName <- function(pathway_name) {
  if (pathway_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") {
    return("TNFA signaling via NFKB")
  }
  if (pathway_name == "HALLMARK_HYPOXIA") {
    return("Hypoxia")
  }
  if (pathway_name == "HALLMARK_CHOLESTEROL_HOMEOSTASIS") {
    return("Cholesterol homeostasis")
  }
  if (pathway_name == "HALLMARK_MITOTIC_SPINDLE") {
    return("Mitotic spindle")
  }
  if (pathway_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING") {
    return("WNT beta catenin signaling")
  }
  if (pathway_name == "HALLMARK_TGF_BETA_SIGNALING") {
    return("TGF beta signaling")
  }
  if (pathway_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING") {
    return("IL16 JAK STAT3 signaling")
  }
  if (pathway_name == "HALLMARK_DNA_REPAIR") {
    return("DNA repair")
  }
  if (pathway_name == "HALLMARK_G2M_CHECKPOINT") {
    return("G2M checkpoint")
  }
  if (pathway_name == "HALLMARK_APOPTOSIS") {
    return("Apoptosis")
  }
  if (pathway_name == "HALLMARK_NOTCH_SIGNALING") {
    return("Notch signaling")
  }
  if (pathway_name == "HALLMARK_ADIPOGENESIS") {
    return("Adipogenesis")
  }
  if (pathway_name == "HALLMARK_ESTROGEN_RESPONSE_EARLY") {
    return("Estrogen response early")
  }
  if (pathway_name == "HALLMARK_ESTROGEN_RESPONSE_LATE") {
    return("Estrogen response late")
  }
  if (pathway_name == "HALLMARK_ANDROGEN_RESPONSE") {
    return("Androgen response")
  }
  if (pathway_name == "HALLMARK_MYOGENESIS") {
    return("Myogenesis")
  }
  if (pathway_name == "HALLMARK_PROTEIN_SECRETION") {
    return("Protein secretion")
  }
  if (pathway_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") {
    return("Interferon alpha response")
  }
  if (pathway_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") {
    return("Interferon gamma response")
  }
  if (pathway_name == "HALLMARK_APICAL_JUNCTION") {
    return("Apical junction")
  }
  if (pathway_name == "HALLMARK_APICAL_SURFACE") {
    return("Apical surface")
  }
  if (pathway_name == "HALLMARK_HEDGEHOG_SIGNALING") {
    return("Hedgehog signaling")
  }
  if (pathway_name == "HALLMARK_COMPLEMENT") {
    return("Complement")
  }
  if (pathway_name == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE") {
    return("Unfolded protein response")
  }
  if (pathway_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING") {
    return("PI3K AKT MTOR signaling")
  }
  if (pathway_name == "HALLMARK_MTORC1_SIGNALING") {
    return("MTORC1 signaling")
  }
  if (pathway_name == "HALLMARK_E2F_TARGETS") {
    return("E2F targets")
  }
  if (pathway_name == "HALLMARK_MYC_TARGETS_V1") {
    return("MYC targets V1")
  }
  if (pathway_name == "HALLMARK_MYC_TARGETS_V2") {
    return("MYC targets V2")
  }
  if (pathway_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") {
    return("Epithelial mesenchymal transition")
  }
  if (pathway_name == "HALLMARK_INFLAMMATORY_RESPONSE") {
    return("Inflammatory response")
  }
  if (pathway_name == "HALLMARK_XENOBIOTIC_METABOLISM") {
    return("Xenobiotic metabolism")
  }
  if (pathway_name == "HALLMARK_FATTY_ACID_METABOLISM") {
    return("Fatty acid metabolism")
  }
  if (pathway_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") {
    return("Oxidative phosphorylation")
  }
  if (pathway_name == "HALLMARK_GLYCOLYSIS") {
    return("Glycolysis")
  }
  if (pathway_name == "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY") {
    return("Reactive oxigen species pathway")
  }
  if (pathway_name == "HALLMARK_P53_PATHWAY") {
    return("P53 pathway")
  }
  if (pathway_name == "HALLMARK_UV_RESPONSE_UP") {
    return("UV response up")
  }
  if (pathway_name == "HALLMARK_UV_RESPONSE_DN") {
    return("UV response down")
  }
  if (pathway_name == "HALLMARK_ANGIOGENESIS") {
    return("Angiogenesis")
  }
  if (pathway_name == "HALLMARK_HEME_METABOLISM") {
    return("Heme metabolism")
  }
  if (pathway_name == "HALLMARK_COAGULATION") {
    return("Coagulation")
  }
  if (pathway_name == "HALLMARK_IL2_STAT5_SIGNALING") {
    return("IL2 STAT5 signaling")
  }
  if (pathway_name == "HALLMARK_BILE_ACID_METABOLISM") {
    return("Bile acide  metabolism")
  }
  if (pathway_name == "HALLMARK_PEROXISOME") {
    return("Peroxisome")
  }
  if (pathway_name == "HALLMARK_ALLOGRAFT_REJECTION") {
    return("Allograft rejection")
  }
  if (pathway_name == "HALLMARK_SPERMATOGENESIS") {
    return("Spermatogenesis")
  }
  if (pathway_name == "HALLMARK_KRAS_SIGNALING_UP") {
    return("KRAS signaling up")
  }
  if (pathway_name == "HALLMARK_KRAS_SIGNALING_DN") {
    return("KRAS signaling down")
  }
  if (pathway_name == "HALLMARK_PANCREAS_BETA_CELLS") {
    return("Pancreas beta cells")
  }
}
