#' @export
GetAllTranscripts <- function(ensembl_mart, chrs = c(1:22, "X", "Y")) {
  all.transcripts <- unique( getBM(attributes = c("external_gene_name", "chromosome_name", 
                                            "start_position", "end_position",
                                            "transcript_start", "transcript_end", "ensembl_transcript_id"),
                             filters="chromosome_name",
                             values=chrs, 
                             mart = ensembl_mart))
  return(all.transcripts)
}

#' @export
GetAllGenes <- function(ensembl_mart, chrs = c(1:22, "X", "Y")) {
  all.genes <- unique( getBM(attributes = c("external_gene_name", "chromosome_name", 
                                            "start_position", "end_position"),
                             filters="chromosome_name",
                             values=chrs, 
                             mart = ensembl_mart))
  return(all.genes)
}

#' @export
ProcessFixedComponentVCF <- function(fix, strip_chr_prefix=FALSE) {
  fix$CHROM <- as.character(fix$CHROM)
  fix$POS <- as.numeric(as.character(fix$POS))
  fix$ALT <- as.character(fix$ALT)
  fix$REF <- as.character(fix$REF)
  if (strip_chr_prefix) {
    chrs <- sapply(fix$CHROM, function(str) {
      ret <- strsplit(str, split="chr")[[1]][2]
      ret
    })
    fix$CHROM <- chrs
  }
  return(fix)
}

#' @export
ConstructGranges <- function(chr, start, offset=0, width=1)
{
  pos <- IRanges::IRanges(start=start-offset, end=start+width)
  rle_ <- S4Vectors::Rle(paste(chr, "", sep=""))
  granges_obj <- GenomicRanges::GRanges(seqnames = rle_, ranges = pos)
  return(granges_obj)
}

#' @export
CombineSingleCellReads <- function(sc_reads_path, sc_reads_pattern = "*.txt", file_separator="\t") {
  #files <- list.files(sc_reads_path, pattern = "*.txt")
  files <- list.files(sc_reads_path, pattern = sc_reads_pattern)
  n_cells <- length(files)
  if (n_cells > 0) {
    df <- read.table(paste(sc_reads_path, files[1], sep="/"), header=T, sep=file_separator)
    ids <- as.character(df$ID)
    n_snvs <- length(ids)
    ret <- data.frame()
    cell_id <- 1
    for (file in files) {
      df <- read.table(paste(sc_reads_path, file, sep="/"), header=T, sep=file_separator)
      if (dim(df)[1] != n_snvs) {
        print(paste("Error: The number of SNVs do not equal", n_snvs))
      }
      # Below if stmt will remove any cells that do not have any reads.
      if (sum(df$REF_COUNT) + sum(df$ALT_COUNT) > 0) {
        cell_name <- paste("c", cell_id, sep="")
        ret <- rbind(ret, data.frame(ID=ids, Cell=cell_name, a = df$REF_COUNT, d = df$ALT_COUNT + df$REF_COUNT, SampleName = file))
        cell_id <- cell_id + 1
      }
    }
    return(ret)
  } 
  stop("Error: No reads found.")
}

# The same function as CombineSingleCellReads but to handle the case where ID is not available.
#' @export
CombineSingleCellReads2 <- function(sc_reads_path, sc_reads_pattern = "*.txt", file_separator = "\t")
{
  files <- list.files(sc_reads_path, pattern = sc_reads_pattern)
  n_cells <- length(files)
  if (n_cells > 0) {
    df <- read.table(paste(sc_reads_path, files[1], sep = "/"),
                     header = T, sep = file_separator)
    n_snvs <- dim(df)[1]
    ids <- paste("s", 0:(n_snvs-1), sep="")
    ret <- data.frame()
    cell_id <- 1
    for (file in files) {
      df <- read.table(paste(sc_reads_path, file, sep = "/"),
                       header = T, sep = file_separator)
      if (dim(df)[1] != n_snvs) {
        print(paste("Error: The number of SNVs do not equal",
                    n_snvs))
      }
      if (sum(df$REF_COUNT) + sum(df$ALT_COUNT) > 0) {
        cell_name <- paste("c", cell_id, sep = "")
        ret <- rbind(ret, data.frame(ID = ids, Cell = cell_name,
                                     chr = df$CHROM, pos = df$POS,
                                     a = df$REF_COUNT, d = df$ALT_COUNT + df$REF_COUNT,
                                     SampleName = file))
        cell_id <- cell_id + 1
      }
    }
    return(ret)
  }
  stop("Error: No reads found.")
}

# The same function as CombineSingleCellReads but to handle the case where the files are spread across different paths.
#' @export
CombineSingleCellReadsMultiPaths <- function(sc_reads_paths, sc_reads_pattern = "*.txt", file_separator = "\t")
{
  files <- c()
  for (sc_reads_path in sc_reads_paths)
  {
    file_paths <- paste(sc_reads_path, list.files(sc_reads_path, pattern = sc_reads_pattern), sep="/")
    files <- c(files, file_paths)
  }
  n_cells <- length(files)
  if (n_cells > 0) {
    df <- read.table(files[1], header = T, sep = file_separator)
    n_snvs <- dim(df)[1]
    ids <- paste("s", 0:(n_snvs-1), sep="")
    ret <- data.frame()
    cell_id <- 1
    for (file in files) {
      df <- read.table(file, header = T, sep = file_separator)
      if (dim(df)[1] != n_snvs) {
        print(paste("Error: The number of SNVs do not equal",
                    n_snvs))
      }
      if (sum(df$REF_COUNT) + sum(df$ALT_COUNT) > 0) {
        cell_name <- paste("c", cell_id, sep = "")
        ret <- rbind(ret, data.frame(ID = ids, Cell = cell_name,
                                     chr = df$CHROM, pos = df$POS,
                                     a = df$REF_COUNT, d = df$ALT_COUNT + df$REF_COUNT,
                                     SampleName = file))
        cell_id <- cell_id + 1
      }
    }
    return(ret)
  }
  stop("Error: No reads found.")
}

#' @export
#' # Hardcode pathway names.
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