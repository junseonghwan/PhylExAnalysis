args = commandArgs(trailingOnly=TRUE)

print(args)

SEED <- as.numeric(args[1])
DATA_PATH <- args[2]
SC_MUT_MATRIX_PATH <- args[3]
MCMC_ITER <- as.numeric(args[4])

SNV_PATH <- paste(DATA_PATH, "bulk.txt", sep="/")
OUTPUT_PATH <- paste(DATA_PATH, "ddClone", sep="/")
if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = T)
}
setwd(OUTPUT_PATH)

#SNV_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/bulk.txt"
#SC_MUT_MATRIX_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/mut_matrix.txt"
#OUTPUT_PATH <- "/Users/seonghwanjun/data/cell-line/bulk/OV2295/genotype/ddClone/"
#MCMC_ITER <- 200
#SEED <- 157

library(ddclone, quietly = T)
library(dplyr, quietly = T)
library(reshape2, quietly = T)

# Load the bulk data.
bulk <- read.table(SNV_PATH, header=T)

if (sum(names(bulk) == "major_cn") == 1 & sum(names(bulk) == "minor_cn")) {
    print("Using major and minor copy numbers.")
    bulkDat <- data.frame("mutation_id" = bulk$ID,
                          "ref_counts" = bulk$d - bulk$b,
                          "var_counts" = bulk$b,
                          "normal_cn" = 2,
                          "minor_cn" = bulk$major_cn,
                          "major_cn" = bulk$minor_cn)
} else {
    print("Copy number information is not specified -- setting major_cn = minor_cn = 1.")
    bulkDat <- data.frame("mutation_id" = bulk$ID,
                          "ref_counts" = bulk$d - bulk$b,
                          "var_counts" = bulk$b,
                          "normal_cn" = 2,
                          "minor_cn" = 1,
                          "major_cn" = 1)
}

sc_mut_matrix <- read.table(SC_MUT_MATRIX_PATH, header=T)
#head(sc_mut_matrix)
# Check that the mutations are correctly ordered.
mean(colnames(sc_mut_matrix) == as.character(bulkDat$mutation_id))

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH, recursive = TRUE)
}
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = sc_mut_matrix, outputPath = OUTPUT_PATH, nameTag = '')
start <- proc.time()
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
                      numOfIterations = MCMC_ITER, thinning = 10, burnIn = 0,
                      seed = SEED)
end <- proc.time()
diff <- end - start
elapsed_time_seconds <- diff[3]

# Output the results.
df <- ddCloneRes$df
output_file <- paste(OUTPUT_PATH, "results.txt", sep="/")
write.table(df, file = output_file, row.names = F, col.names = T, quote=F)

timing_file <- paste(OUTPUT_PATH, "timing.txt", sep="/")
write.table(paste(elapsed_time_seconds, "seconds"), file = timing_file, row.names = F, col.names = F, quote=F)

# We may remove cells that have mutation at exactly one loci as these would more than likely to confuse ddClone.
#sc_mut_matrix2 <- sc_mut_matrix[rowSums(sc_mut_matrix) > 1,]
#dim(sc_mut_matrix2)
#ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = sc_mut_matrix2, outputPath = OUTPUT_PATH, nameTag = '')
#ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
#                      outputPath = OUTPUT_PATH, tumourContent = 1.0,
#                      numOfIterations = MCMC_ITER, thinning = 10, burnIn = 0,
#                      seed = SEED)

# Output the results.
#df <- ddCloneRes$df
#output_file <- paste(OUTPUT_PATH, "results2.txt", sep="/")
#write.table(df, file = output_file, row.names = F, col.names = T, quote=F)
