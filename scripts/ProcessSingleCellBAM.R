args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args) != 3) {
  stop("Path to loci file, sorted single cell bam file with index file in bam.bai form, and output path are necessary to run the program.")
}

print("Begin extracting reads in R.")

library(dplyr)
library(GenomicAlignments)
library(GenomicRanges)
library(reshape2)
library(Rsamtools)
library(vcfR)

loci_file <- args[1]
bam_file <- args[2]
sc_output_file <- args[3]

#loci_file <- "~/PhylExAnalysis/_temp/loci.txt"
#bam_file <- "~/PhylExAnalysis/_temp/BCSA2_P8Aligned.sortedByCoord.out.bam"
#sc_output_file <- "~/PhylExAnalysis/_temp/BCSA2_P8.txt"

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A","C","G","T")

snv_loci <- read.table(loci_file, header=T, sep="\t")
snv_loci$REF_COUNT <- 0
snv_loci$ALT_COUNT <- 0
snv_loci$CHROM <- as.character(snv_loci$CHROM)

# We will strip off chr prefix if it is there.
if (startsWith(snv_loci$CHROM[1], "chr")) {
  snv_loci$CHROM <- unlist(lapply(as.character(snv_loci$CHROM), function(str) {
    substr(str, 4, nchar(str))
  }))
}

# The reason why width = 0 is used: pileup will count both the POS and POS+1 even though
# we only want it at POS.
somatic.gr <- ConstructGranges(snv_loci$CHROM, snv_loci$POS, width = 0)

# get counts at somatic SNV loci
sbp <- ScanBamParam(which = somatic.gr)
p_param <- PileupParam(distinguish_nucleotides = TRUE,
                       distinguish_strands = FALSE,
                       include_insertions = TRUE)
res_somatic <- pileup(file=bam_file,
                      scanBamParam = sbp,
                      pileupParam = p_param)
if (dim(res_somatic)[1] > 0) {
  df <- data.frame(seqnames = as.character(res_somatic[,1]),
                   pos = as.numeric(as.character(res_somatic[,2])),
                   nucleotide = as.character(res_somatic[,3]),
                   count = as.numeric(as.character(res_somatic[,4])),
                   stringsAsFactors = F)
  df$nucleotide <- factor(df$nucleotide, nucleotides)
  df_somatic <- reshape2::dcast(df, seqnames + pos ~ nucleotide, fun.aggregate = sum, drop = FALSE)
  
  if (sum(colnames(df_somatic) == "-") > 0) {
    col_to_remove<-which(colnames(df_somatic) == "-")
    df_somatic<-df_somatic[,-col_to_remove]
  }
  n <- dim(df_somatic)[2]
  df_somatic$Total <- rowSums(df_somatic[,3:n])
  
  nucleotide_counts <- rowSums(df_somatic[,3:n] > 0)
  print("How many biallelic sites?")
  print(paste(sum(nucleotide_counts > 1), sum(nucleotide_counts > 0), sep="/"))
  
  # get REF and VAR counts
  nn <- dim(df_somatic)[1]
  idx <- paste(snv_loci$CHROM, snv_loci$POS, sep=":") %in% paste(df_somatic$seqnames, df_somatic$pos, sep=":")
  for (i in 1:nn) {
    idx <- which(snv_loci$CHROM == df_somatic[i,"seqnames"] & snv_loci$POS == df_somatic[i,"pos"])
    if (length(idx) > 0) {
      ref <- as.character(snv_loci[idx,"REF"])
      alt <- as.character(snv_loci[idx,"ALT"])
      if ((ref %in% nucleotides) & (alt %in% nucleotides)) {
        alt_idx <- which(alt == nucleotides)
        ref_idx <- which(ref == nucleotides)
        snv_loci[idx,"REF_COUNT"] <- df_somatic[i,2+ref_idx]
        snv_loci[idx,"ALT_COUNT"] <- df_somatic[i,2+alt_idx]
      }
    }
  }
}
head(snv_loci)
write.table(snv_loci, file=sc_output_file, sep="\t", row.names=F, quote=F)
