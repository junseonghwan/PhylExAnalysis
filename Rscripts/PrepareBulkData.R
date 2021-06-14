args = commandArgs(trailingOnly=TRUE)
print(args)

VCF_PATH <- as.character(args[1])
CNV_PATH <- as.character(args[2])
OUTPUT_PATH <- as.character(args[3])

#VCF_PATH <- "~/data/TNBC/BCSA1/BCSA1_filtered_hg38.vcf"
#CNV_PATH <- "~/data/TNBC/BCSA1/"
#OUTPUT_PATH <- "~/PhylExAnalysis/_temp/"

# path to TitanCNA output file within each bulk region.
CNV_SUFFIX_PATH <- "cna/results/titan/hmm/"
OPTIMAL_CNV_SOLN_FILE <- paste(CNV_SUFFIX_PATH, "optimalClusterSolution.txt", sep="/")

GERMLINE_COLUMN_IDX <- 1
MIN_VAF <- 0.03
APPLY_PASS_FILTER <- TRUE # Set to TRUE if there is a filter to be used.
#FILTER_HLA_GENES <- TRUE

library(biomaRt)
library(PhylExR)
library(vcfR)

vcf <- read.vcfR(VCF_PATH)

fix <- data.frame(vcf@fix)
fix <- ProcessFixedComponentVCF(fix, strip_chr_prefix = T)

allelic_depth <- extract.gt(vcf, element = "AD")
allelic_depth <- data.frame(allelic_depth)

# Select SNVs
filter_idx2 <- (nchar(fix$ALT) == 1) & (nchar(fix$REF) == 1)

# Pass filters
if (APPLY_PASS_FILTER) {
  filter_idx1 <- (fix$FILTER == "PASS")
  filter_idx <- filter_idx1 & filter_idx2
} else {
  filter_idx <- filter_idx2
}

fix <- fix[filter_idx,]
allelic_depth <- allelic_depth[filter_idx,]

dim(fix)
dim(allelic_depth)

region_names <- names(allelic_depth)
ref_counts <- data.frame(matrix(0, ncol = length(region_names), nrow = sum(filter_idx)))
alt_counts <- data.frame(matrix(0, ncol = length(region_names), nrow = sum(filter_idx)))
for (i in 1:length(region_names)) {
  ref_alt <- sapply(as.character(allelic_depth[,region_names[i]]), function(ad) {
    ref <- strsplit(ad, split=",")[[1]][1]
    alt <- strsplit(ad, split=",")[[1]][2]
    return(matrix(c(ref, alt), ncol=2, byrow = T))
  })
  ref_alt <- t(ref_alt)
  ref_counts[,i] <- as.numeric(ref_alt[,1])
  alt_counts[,i] <- as.numeric(ref_alt[,2])
}

names(ref_counts) <- region_names
names(alt_counts) <- region_names
depth <- ref_counts + alt_counts

# Filter out by MIN_VAF
vaf <- alt_counts[,-GERMLINE_COLUMN_IDX]/depth[,-GERMLINE_COLUMN_IDX]
region_count <- dim(vaf)[2]
filter_idx <- (rowSums(vaf >= MIN_VAF) == region_count)

fix <- fix[filter_idx,]
depth <- depth[filter_idx,]
alt_counts <- alt_counts[filter_idx,]

d <- apply(depth[,-GERMLINE_COLUMN_IDX], 1, function(row) {
  ret <- paste(row, collapse = ",")
  ret
})
b <- apply(alt_counts[,-GERMLINE_COLUMN_IDX], 1, function(row) {
  ret <- paste(row, collapse = ",")
  ret
})

# Construct multi-region data frame.
snv_count <- dim(fix)[1]
mut_ids <- paste("s", 0:(snv_count-1), sep="")

# Add copy number information.
snv.gr <- ConstructGranges(chr = fix[,c("CHROM")], start = fix[,c("POS")], width = 0)
major_cn <- matrix(1, nrow = snv_count, ncol = region_count)
minor_cn <- matrix(1, nrow = snv_count, ncol = region_count)
colnames(major_cn) <- region_names[-GERMLINE_COLUMN_IDX]
colnames(minor_cn) <- region_names[-GERMLINE_COLUMN_IDX]
for (i in 1:region_count) {
  region <- region_names[-GERMLINE_COLUMN_IDX][i]
  REGION_PATH <- paste(CNV_PATH, region, sep="/")
  opt_soln <- read.table(paste(REGION_PATH, OPTIMAL_CNV_SOLN_FILE, sep="/"), header=F, skip = 1)
  CNA_PATH <- paste(REGION_PATH, "/", CNV_SUFFIX_PATH, "/optimalClusterSolution/", opt_soln$V2, ".segs.txt", sep="")

  cna <- read.table(CNA_PATH, header=T)
  # Strip `chr`
  cna_chrs <- gsub("chr", "", cna$Chromosome)
  cna.gr <- ConstructGranges(cna_chrs, cna$Start_Position.bp., width = cna$End_Position.bp. - cna$Start_Position.bp.)
  ret <- findOverlaps(snv.gr, cna.gr)
  major_cn[ret@from,i] <- cna[ret@to,"MajorCN"]
  minor_cn[ret@from,i] <- cna[ret@to,"MinorCN"]
}
major_cn <- apply(major_cn, 1, function(row) {
  paste(row, collapse = ",")
})
minor_cn <- apply(minor_cn, 1, function(row) {
  paste(row, collapse = ",")
})


bulk <- data.frame(ID = mut_ids,
                   b = b, d = d,
#                   major_cn=paste(rep(1, region_count), collapse = ","), minor_cn=paste(rep(1, region_count), collapse=","))
                  major_cn=major_cn, minor_cn=minor_cn)
if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH, recursive = T)
}
write.table(bulk, file = paste(OUTPUT_PATH, "bulk.txt", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

# We will output loci file needed to look up SNV to their CHR, POS, REF, ALT
loci <- data.frame(ID=mut_ids, fix[,c("CHROM", "POS", "REF", "ALT")])
write.table(loci, file = paste(OUTPUT_PATH, "loci.txt", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
