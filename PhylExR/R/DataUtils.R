
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
        ret <- rbind(ret, data.frame(ID=ids, Cell=cell_name, a = df$REF_COUNT, d = df$ALT_COUNT + df$REF_COUNT))
        cell_id <- cell_id + 1
      }
    }
    return(ret)
  } 
  stop("Error: No reads found.")
}
