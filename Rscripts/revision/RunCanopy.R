args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
REP_PATH <- args[2]
#SEED <- 1
#REP_PATH <- "~/Dropbox/seong/PhylExNatComm/zenodo/binary_cn/rep0/"
set.seed(SEED)

#library(cardelino)
library(Canopy)
#library(PhylExR)

# Run Canopy.

####################################################################
# Constants
K_BEGIN <- 3
K_END <- 12
K <- K_BEGIN:K_END
numchains <- 4
min_iter <- 10000
max_iter <- 100000
thinning <- 10
burnin <- 1000
projname <- "canopy"

####################################################################

####################################################################
# Run Canopy. 
case_no <- 0
CASE_PATH <- paste(REP_PATH, "/case", case_no, sep="")
SNV_PATH <- paste(CASE_PATH, "genotype_ssm.txt", sep="/")
SNV_GT_PATH <- paste(CASE_PATH, "datum2node.tsv", sep="/")
OUTPUT_PATH <- paste(CASE_PATH, "canopy", sep="/")

if (!dir.exists(OUTPUT_PATH)) {
  dir.create(OUTPUT_PATH)
}
setwd(OUTPUT_PATH)

bulk <- read.table(SNV_PATH, header=T, sep="\t", as.is = T)
snv_gt <- read.table(SNV_GT_PATH, header = F, sep="\t", as.is = T)
snv_count <- dim(bulk)[1]
cnv_count <- dim(bulk)[1]

if (!is.numeric(bulk$b) & !is.numeric(bulk$d)) {
  bulk$b <- as.character(bulk$b)
  bulk$d <- as.character(bulk$d)
  R <- matrix(as.numeric(unlist(strsplit(bulk$b, ","))), nrow = snv_count, byrow = T)
  D <- matrix(as.numeric(unlist(strsplit(bulk$d, ","))), nrow = snv_count, byrow = T)
  X <- D - R
  
  WM <- matrix(as.numeric(unlist(strsplit(bulk$major_cn, ","))), nrow = cnv_count, byrow = T)
  Wm <- matrix(as.numeric(unlist(strsplit(bulk$minor_cn, ","))), nrow = cnv_count, byrow = T)
} else {
  X <- as.matrix(bulk$d - bulk$b)
  R <- as.matrix(bulk$b)
  WM <- as.matrix(bulk$major_cn)
  Wm <- as.matrix(bulk$minor_cn)
}

region_count <- dim(X)[2]
rownames(X) <- bulk$ID
rownames(R) <- bulk$ID
colnames(X) <- paste("Sample", 1:region_count, sep="")
colnames(R) <- paste("Sample", 1:region_count, sep="")

rownames(WM) <- paste("c", 1:cnv_count, sep="")
rownames(Wm) <- paste("c", 1:cnv_count, sep="")
colnames(WM) <- paste("Sample", 1:region_count, sep="")
colnames(Wm) <- paste("Sample", 1:region_count, sep="")

no_cna <- mean(WM == 1) == 1 &  mean(Wm == 1) == 1
if (no_cna) {
  sampchain <- canopy.sample.nocna(R, X, K, numchains, projectname = projname,
                                   max.simrun = max_iter, min.simrun = min_iter,
                                   writeskip = thinning, cell.line=FALSE, plot.likelihood=F)
} else {
  Y <- matrix(0, cnv_count, ncol = cnv_count + 1)
  diag(Y[,-1]) <- 1
  sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = rep(0.01, cnv_count),
                            epsilonm = rep(0.01, cnv_count), Y = Y, K = K, numchain = numchains,
                            max.simrun = max_iter, min.simrun = min_iter,
                            writeskip = thinning, projectname = projname, cell.line = FALSE,
                            plot.likelihood = TRUE)
}

bic <- canopy.BIC(sampchain = sampchain, projectname = projname, K = K,
                  numchain = numchains, burnin = burnin, thin = thinning, pdf = TRUE)
optK <- K[which.max(bic)]
post <- canopy.post(sampchain = sampchain, projectname = projname, K = K,
                    numchain = numchains, burnin = burnin, thin = thinning,
                    optK = optK, post.config.cutoff = 0.01)
samptreethin <- post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik <- post[[2]]   # likelihoods of trees in samptree
config <- post[[3]]
config.summary <- post[[4]]
print(config.summary)

config.i <- config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree <- canopy.output(post, config.i, C=NULL)
saveRDS(output.tree, file = paste(OUTPUT_PATH, "canopy_output.RDS", sep="/"))
write.table(output.tree$Z, file = paste(OUTPUT_PATH, "Z.txt", sep="/"), row.names = T, col.names = T, quote = F)
write.table(output.tree$P, file = paste(OUTPUT_PATH, "P.txt", sep="/"), row.names = T, col.names = T, quote = F)
####################################################################
