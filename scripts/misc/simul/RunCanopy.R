args = commandArgs(trailingOnly=TRUE)
print(args)
SEED <- as.numeric(args[1])
REP_PATH <- args[2]
K_BEGIN <- as.numeric(args[3])
K_END <- as.numeric(args[4])

#SEED <- 157
#REP_PATH <- "/Users/seonghwanjun/data/simul/quadternary_cn/rep0/case0/"
#K_BEGIN <- 3
#K_END <- 12
#OUTPUT_PATH <- "/Users/seonghwanjun/data/simul/quadternary_cn/rep0/case0/canopy/"

SNV_PATH <- paste(REP_PATH, "genotype_ssm.txt", sep="/")
OUTPUT_PATH <- paste(REP_PATH, "canopy", sep="/")

library(Canopy)

set.seed(SEED)

if (!dir.exists(OUTPUT_PATH)) {
    dir.create(OUTPUT_PATH)
}
setwd(OUTPUT_PATH)

bulk <- read.table(SNV_PATH, header=T, sep="\t", as.is = T)
snv_count <- dim(bulk)[1]
cnv_count <- dim(bulk)[1]

# Check if the data is multi-regional
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

K <- K_BEGIN:K_END
numchains <- 4
min_iter <- 10000
max_iter <- 100000
thinning <- 100
projname <- "canopy"

if (mean(WM == 1) == 1 &  mean(Wm == 1) == 1) {
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
write.table(output.tree$Z, file = paste(OUTPUT_PATH, "Z.txt", sep="/"), row.names = T, col.names = T, quote = F)
write.table(output.tree$P, file = paste(OUTPUT_PATH, "P.txt", sep="/"), row.names = T, col.names = T, quote = F)
write.tree(output.tree, file = paste(OUTPUT_PATH, "tree.newick", sep="/"))
pdf.name <- paste(OUTPUT_PATH, 'highest_likelihood.pdf', sep='/')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, pdf = FALSE)

# Output clustering of mutations for clustering accuracy evaluation.
clone_name <- paste(output.tree$sna[,2], output.tree$sna[,3], sep="_")
predicted <- cbind(ID=rownames(output.tree$sna), Cluster=clone_name)
write.table(predicted, paste(OUTPUT_PATH, "predicted.csv", sep="/"), sep=",", col.names=T, row.names=F, quote=F)

snv_count <- dim(output.tree$Z)[1]
A <- matrix(0, snv_count, snv_count)
for (i in 1:snv_count) {
    clone_i <- sort(which(output.tree$Z[i,] == 1))
    for (j in 1:snv_count) {
        clone_j <- sort(which(output.tree$Z[j,] == 1))
        if (!setequal(clone_i, clone_j)) {
            A[i,j] <- as.numeric(all(clone_j %in% clone_i))
        }
    }
}

write.table(A, paste(OUTPUT_PATH, "ancestral_matrix.csv", sep="/"), row.names = F, quote=F, col.names = F)
