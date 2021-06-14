rm(list=ls())
library(MARATHON)
library(falcon)
library(VariantAnnotation)

chrs <- c(1:22, "X", "Y")
nucleotides <- c("A", "C", "G", "T")
READ_THRESHOLD <- 30
LENGTH_THRESHOLD <- 10^6  # Threshold for length of segments, in base pair.
DELTA_CN_THRESHOLD <- 0.3  # Threshold of absolute copy number difference between consecutive segments.


vcfFile <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/results/variants/somatic.snvs.vcf"
vcfObj <- readVcf(vcfFile)
names(geno(vcfObj))
names(info(vcfObj))
header(vcfObj)

# We have an array.
sapply(geno(vcfObj), class)
dim(geno(vcfObj)$AU)
# Convert it to data frame. This yields 4 columns. The first two columns are the read counts for tier 1 (see Strelka).
# The lasst two columns are the read counts for tier 2 (we don't use).
AU <- as.data.frame(geno(vcfObj)$AU)
CU <- as.data.frame(geno(vcfObj)$CU)
GU <- as.data.frame(geno(vcfObj)$GU)
TU <- as.data.frame(geno(vcfObj)$TU)

normal_sample <- data.frame(A = AU$NORMAL.1, C = CU$NORMAL.1, G = GU$NORMAL.1, "T" = TU$NORMAL.1)
tumor_sample <- data.frame(A = AU$TUMOR.1, C = CU$TUMOR.1, G = GU$TUMOR.1, "T" = TU$TUMOR.1)

REF <- as.data.frame(ref(vcfObj))
ALT <- as.data.frame(unlist(alt(vcfObj)))
dim(REF) == dim(ALT)
site_count <- dim(REF)[1]

loc_gt <- matrix(unlist(strsplit(rownames(AU), "_")), ncol = 2, byrow = T)
loc <- loc_gt[,1]
chr_pos <- data.frame(matrix(unlist(strsplit(loc, ":")), ncol = 2, byrow = T))
names(chr_pos) <- c("CHR", "POS")

# Prepare read matrix with 4 columns:
# 'AT','BT','AN','BN'
# Number of A and B alleles for tumor in the first 2 columns.
# Number of A and B alleles for normal in the last 2 columns.
read_matrix.df <- as.data.frame(matrix(0, nrow = site_count, ncol = 4))
names(read_matrix.df) <- c("AT", "BT", "AN", "BN")
read_matrix.df$CHR <- chr_pos$CHR
read_matrix.df$POS <- as.numeric(as.character(chr_pos$POS))
head(read_matrix.df)
for (i in 1:length(nucleotides)) {
    base <- nucleotides[i]
    # Retrieve read counts matching the reference allele from tumor_sample and normal sample.
    idx <- which(REF[,1] == base)
    read_matrix.df[idx,1] <- tumor_sample[idx,base]
    read_matrix.df[idx,3] <- normal_sample[idx,base]

    # Retrieve read counts matching the alternate allele from tumor_sample and normal sample.
    idx <- which(ALT[,1] == base)
    read_matrix.df[idx,2] <- tumor_sample[idx,base]
    read_matrix.df[idx,4] <- normal_sample[idx,base]
}

rdep <- sum(read_matrix.df$AT + read_matrix.df$BT)/sum(read_matrix.df$AN + read_matrix.df$BN)

# Run Falcon for each chromosome.
out.df <- data.frame()
for (chr in chrs) {
    temp <- subset(read_matrix.df, CHR == chr)
    temp <- temp[(rowSums(temp[,1:2]) >= READ_THRESHOLD & rowSums(temp[,3:4]) >= READ_THRESHOLD),]
    dim(temp)
    tauhat <- getChangepoints(temp[,1:4])
    cn <- getASCN(temp[,1:4], tauhat=tauhat, rdep = rdep, threshold = 0.3)
    if(length(tauhat)>0){
        falcon.qc.list <- falcon.qc(readMatrix = temp[,1:4],
                                    tauhat = tauhat,
                                    cn = cn,
                                    st_bp = temp$POS,
                                    end_bp = temp$POS,
                                    rdep = rdep,
                                    length.thres = LENGTH_THRESHOLD,
                                    delta.cn.thres = DELTA_CN_THRESHOLD)
        tauhat2 <- falcon.qc.list$tauhat
        cn2 <- falcon.qc.list$cn
    } else {
        tauhat2 <- tauhat
        cn2 <- cn
    }

    falconOutput <- falcon.output(readMatrix = temp[,1:4],
                                  tauhat = tauhat2,
                                  cn = cn2,
                                  st_bp = temp$POS,
                                  end_bp = temp$POS,
                                  nboot = 5000)
    out <- data.frame(CHR=chr, falconOutput[,-c(1,2)])
    out.df <- rbind(out.df, out)
}
write.table(out.df, "/Users/seonghwanjun/data/cell-line/bulk/OV2295/falcon.txt", sep="\t", quote=F, row.names = F, col.names = T)
