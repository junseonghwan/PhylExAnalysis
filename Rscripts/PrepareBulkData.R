args = commandArgs(trailingOnly=TRUE)
print(args)

VCF_PATH <- as.character(args[1])
CNV_PATH <- as.character(args[2])
OUTPUT_PATH <- as.character(args[3])
GERMLINE_COLUMN_IDX <- 1
MIN_VAF <- 0.03

#VCF_PATH <- "~/data/TNBC/BCSA2/BCSA2_filtered_liftoverhg19.vcf"
#OUTPUT_PATH <- "~/PhylExAnalysis/_temp/"

library(PhylExR)
library(vcfR)

vcf <- read.vcfR(VCF_PATH)

fix <- data.frame(vcf@fix)
fix <- ProcessFixedComponentVCF(fix, strip_chr_prefix = T)

allelic_depth <- extract.gt(vcf, element = "AD")
allelic_depth <- data.frame(allelic_depth)

# Pass filters
filter_idx1 <- (fix$FILTER == "PASS")

# Select SNVs
filter_idx2 <- (nchar(fix$ALT) == 1) & (nchar(fix$REF) == 1)

filter_idx <- filter_idx1 & filter_idx2

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

bulk <- data.frame(ID = mut_ids, 
                   b = b, d = d,
                   major_cn=paste(rep(1, region_count), collapse = ","), minor_cn=paste(rep(1, region_count), collapse=","))
if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH, recursive = T)
}
write.table(bulk, file = paste(OUTPUT_PATH, "bulk.txt", sep=""), sep="\t", col.names = T, row.names = F, quote = F)

# We will output loci file needed to look up SNV to their CHR, POS, REF, ALT
loci <- data.frame(ID=mut_ids, fix[,c("CHROM", "POS", "REF", "ALT")])
write.table(loci, file = paste(OUTPUT_PATH, "loci.txt", sep=""), sep="\t", col.names = T, row.names = F, quote = F)
