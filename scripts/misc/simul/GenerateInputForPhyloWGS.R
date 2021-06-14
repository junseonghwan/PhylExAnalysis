args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    print("Provide a config file delimeted by tab or space.")
    stop("Terminating...")
}

print(args)

library(GenomicRanges)
library(ScRNAClone)
library(vcfR)

library(dplyr)
library(reshape2)
library(Rsamtools)

REP_PATH <- as.character(args[1])
#REP_PATH <- "/Users/seonghwanjun/data/simul/quadternary_multiregion/rep0/case0"
#REP_PATH <- "/Users/seonghwanjun/data/simul/binary/rep0/case0"
PWGS_SNV_OUT <- paste(REP_PATH, "pwgs_snv.txt", sep="/")
PWGS_CNV_OUT <- paste(REP_PATH, "pwgs_cnv.txt", sep="/")
SNV_FILE <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
snvs <- read.table(SNV_FILE, header=T, sep="\t", as.is=T)

if (!is.numeric(snvs$b[0])) {
    alt_counts <- t(as.matrix(sapply(snvs$b, function(x) {
        as.numeric(strsplit(x, split = ",")[[1]])
    })))
    total_counts <- t(as.matrix(sapply(snvs$d, function(x) {
        as.numeric(strsplit(x, split = ",")[[1]])
    })))
    ref_counts <- total_counts - alt_counts
    snvs$a <- apply(ref_counts, 1, function(row) {
        paste(row, collapse = ",")
    })
} else {
    snvs$a <- snvs$d - snvs$b
}

pwgs.df <- data.frame(id=as.character(snvs$ID), gene=paste(snvs$CHR, snvs$POS, sep="_"), a=snvs$a, d=snvs$d, mu_r=0.999, mu_v=0.499)

write.table(pwgs.df, PWGS_SNV_OUT, quote=F, row.names = F, col.names = T, sep="\t")

# Write an empty CNV file.
pwgs_cnv <- data.frame(cnv=character(0), a=numeric(0), d=numeric(0), ssms=character(0), physical_cnvs=character(0))
write.table(pwgs_cnv, PWGS_CNV_OUT, col.names = T, row.names = F, quote=F, sep="\t")
