
args = commandArgs(trailingOnly=TRUE)
print(args)

sc_reads_path <- args[1]
loci_path <- args[2]
output_path <- args[3]

#sc_reads_path <- "~/PhylExAnalysis/_temp/scRNA/"
#loci_path <- "~/PhylExAnalysis/_temp/loci.txt"
#output_path <- "~/PhylExAnalysis/_temp/"

library(PhylExR)

loci <- read.table(loci_path, header = T)

sc <- CombineSingleCellReads(sc_reads_path, file_separator = "\t")
sc_ <- sc[!(sc$a == 0 & sc$d == 0),]

if(!dir.exists(output_path)) {
  dir.create(output_path)
}

write.table(sc_, file = paste(output_path, "sc.txt"), sep="\t", row.names = F, quote = F, col.names = T)

# Estimate hyperparameters.
sc_hp <- EstimateHyperParameters(loci$ID, sc)
write.table(sc_hp, file = paste(output_path, "sc_hp.txt"), sep="\t", row.names = F, quote = F, col.names = T)


