args = commandArgs(trailingOnly=TRUE)

print(args)

REP_PATH <- args[1]
#REP_PATH <- "/Users/seonghwanjun/data/temp/rep0/case0/"
SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")

ssm <- read.table(SNV_PATH, header=T, as.is = TRUE)

B <- matrix(unlist(strsplit(ssm$b, split = ",")), ncol = 3, byrow=T)
D <- matrix(unlist(strsplit(ssm$d, split = ",")), ncol = 3, byrow=T)
major_cn <- matrix(unlist(strsplit(ssm$major_cn, split = ",")), ncol = 3, byrow=T)
minor_cn <- matrix(unlist(strsplit(ssm$minor_cn, split = ",")), ncol = 3, byrow=T)

for (i in 1:dim(B)[2]) {
    single_region_ssm <- cbind(ID=ssm$ID, b=B[,i], d=D[,i], minor_cn=major_cn[,i], minor_cn=minor_cn[,i])
    write.table(single_region_ssm, paste(REP_PATH, paste("single_region_genotype_ssm", i, ".txt", sep=""), sep="/"), col.names = T, row.names=F, quote=F, sep="\t")
}
