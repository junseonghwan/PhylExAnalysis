library(biomaRt)
library(dplyr)
library(GenomicRanges)
library(PhylExR)

all_loci <- read.table("data/HGSOC/ov2295_clone_snvs.csv", header=T, sep=",")
all_loci$loc <- paste(all_loci$chrom, all_loci$coord, sep=":")
all_loci <- all_loci[!duplicated(all_loci$loc),]
exon_loci <- read.table("data/HGSOC/ov2295_exonic_loci.txt", header=T)


grch37_exons <- read.table("~/Dropbox/seong/data/GRCh37_exons.bed", header = T)
exon_ranges <- PhylExR::ConstructGranges(chr = grch37_exons$chr, start = grch37_exons$start, width = grch37_exons$end - grch37_exons$start)

loci_ranges <- PhylExR::ConstructGranges(chr = all_loci$chrom, start = all_loci$coord, width = 0)
ret2 <- GenomicRanges::findOverlaps(loci_ranges, exon_ranges)
length(unique(ret2@from))

#all_loci[ret2@from,]
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", GRCh = 37)
bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
            filters = 'ensembl_gene_id',
            values = grch37_exons[ret2@to,"ensembl_gene_id"], 
            mart = ensembl)
write.table(x = sort(unique(bm$hgnc_symbol)), file = "temp.txt", quote = F, row.names = F, col.names = F)

subset(all_loci[ret2@from,], chrom == 7 & coord > 40174575 & coord < 40900108) 

xx <- grch37_exons[ret2@to,"ensembl_gene_id"]
xx[xx == "ENSG00000175600"]
bm[bm$ensembl_gene_id == "ENSG00000175600",]

phylex_snvs <- read.table("data/HGSOC_SS3/bulk.txt", header = T)
chain_no <- PhylExR::FindBestRep("data/HGSOC_SS3/phylex_1_02_05/phylex/", chains = 0:19)
chain_path <- paste("data/HGSOC_SS3/phylex_1_02_05/phylex/chain", chain_no, sep="")
datum2node <- read.table(paste(chain_path, "/joint/tree0/datum2node.tsv", sep=""), header = F)
names(datum2node) <- c("ID", "Node")

phylex_snvs <- read.table("data/HGSOC_SS3/phylex/chain", header = T)
phylex_snvs <- dplyr::left_join(phylex_snvs, datum2node, by = "ID")
phylex_snvs_ranges <- ConstructGranges(phylex_snvs$CHR, phylex_snvs$POS, width = 0)
ret3 <- findOverlaps(phylex_snvs_ranges, exon_ranges)
phylex_exons <- grch37_exons[ret3@to,]
phylex_exons$from_idx <- ret3@from
phylex_exons$to_idx <- ret3@to

head(listAttributes(ensembl), 20)
bm3 <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'external_gene_name'),
            filters = 'ensembl_gene_id',
            values = phylex_exons$ensembl_gene_id, 
            mart = ensembl)
phylex_exons_ <- dplyr::left_join(phylex_exons, bm3)
#phylex_exons_ <- phylex_exons_[phylex_exons_$hgnc_symbol != "",]
is_dup <- duplicated(paste(phylex_exons_$from_idx, phylex_exons_$ensembl_gene_id, sep=":"))
phylex_exons_ <- phylex_exons_[!is_dup,]
head(phylex_exons_)
phylex_exons_$ID <- phylex_snvs[phylex_exons_$from_idx,"ID"]

unique(paste(phylex_exons_$from_idx, phylex_exons_$external_gene_name, sep=":"))
sort(phylex_exons_$external_gene_name)

phylex_snvs_ <- left_join(phylex_snvs, phylex_exons_, by = "ID")
head(phylex_snvs_)
phylex_snvs_[order(phylex_snvs_$Node),c("ID", "Node", "external_gene_name")]

library(xlsx)
nano <- xlsx::read.xlsx(file = "data/LBL-10025_nCounter_PanCancer_Human_Pathways_Panel_Gene_List.xlsx", sheetName = "Gene to Pathway Mapping", startRow = 2)
bm[bm$hgnc_symbol %in% nano$Gene,]
