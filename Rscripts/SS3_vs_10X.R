library(PhylExR)

SC_READS_SS3 <- "~/data/cell-line/smart-seq3/Reads/"
SC_READS_10X <- "~/data/cell-line/10X/scRNA-by-region/"
OV2295 <- "CHIP0013"
OV2295R <- c("CHIP0002", "CHIP0012")
LOCI_FILE <- "data/HGSOC/ov2295_exonic_loci.txt"
MIN_READS <- 1
base_size <- 12

loci <- read.table(LOCI_FILE, header = T)
gt <- read.table("data/HGSOC/ov2295_clone_snvs.csv", header = T, sep=",")
gt$loc <- paste(gt$chrom, gt$coord, sep=":")
gt_labels <- gt %>% group_by(loc) %>% summarise(clone = paste(clone_id[which(is_present == 1)], collapse = "_"))

### Takes long time ###
### Skip below and load the data if it has already been processed. ###
ov2295_dat <- PhylExR::CombineSingleCellReads2(sc_reads_path = paste(SC_READS_10X, OV2295, sep="/"))
ov2295R_dat1 <- PhylExR::CombineSingleCellReads2(sc_reads_path = paste(SC_READS_10X, OV2295R[1], sep="/"))
ov2295R_dat2 <- PhylExR::CombineSingleCellReads2(sc_reads_path = paste(SC_READS_10X, OV2295R[2], sep="/"))
ov2295R_dat <- rbind(ov2295R_dat1, ov2295R_dat2)

saveRDS(ov2295_dat, file = "data/HGSOC/ov2295_10x_reads.RDS")
saveRDS(ov2295R_dat1, file = "data/HGSOC/ov2295R1_10x_reads.RDS")
saveRDS(ov2295R_dat2, file = "data/HGSOC/ov2295R2_10x_reads.RDS")
###

### load data ###
ov2295_dat <- readRDS("data/HGSOC/ov2295_10x_reads.RDS")
ov2295R_dat1 <- readRDS("data/HGSOC/ov2295R1_10x_reads.RDS")
ov2295R_dat2 <- readRDS("data/HGSOC/ov2295R2_10x_reads.RDS")
###

# Identify loci with coverage to use for analysis.
ov2295_dat$loc <- paste(ov2295_dat$chr, ov2295_dat$pos, sep=":")
ov2295R_dat1$loc <- paste(ov2295R_dat1$chr, ov2295R_dat1$pos, sep=":")
ov2295R_dat2$loc <- paste(ov2295R_dat2$chr, ov2295R_dat2$pos, sep=":")
ov2295R_dat <- rbind(ov2295R_dat1, ov2295R_dat2)

ov2295_coverage <- ov2295_dat %>% group_by(loc) %>% summarise(n_b = sum(d - a > 0))
ov2295R1_coverage <- ov2295R_dat1 %>% group_by(loc) %>% summarise(n_b = sum(d - a > 0))
ov2295R2_coverage <- ov2295R_dat2 %>% group_by(loc) %>% summarise(n_b = sum(d - a > 0))

coverage <- left_join(left_join(ov2295_coverage, ov2295R1_coverage, by = "loc"), ov2295R2_coverage, by = "loc")
coverage.df <- as.data.frame(coverage)
coverage.df$b <- rowSums(coverage.df[,-1])

coverage.df <- left_join(coverage.df, gt_labels)
loci_with_coverage <- subset(coverage.df, b >= 2)
dim(loci_with_coverage)
loci_with_coverage

# Check that the loci that was identified by combining all cells together match above loci.
loci_old <- read.table("data/HGSOC_10X/loci.txt", header=T)
dim(loci_old)
head(loci_old)
loci_old$loc <- paste(loci_old$CHROM, loci_old$POS, sep=":")
mean(loci_with_coverage$loc %in% loci_old$loc)
mean(loci_old$loc %in% loci_with_coverage$loc)

sc <- read.table("data/HGSOC_10X/sc.txt", header=T)
temp <- sc %>% group_by(Cell) %>% summarise(n_b = sum(d - a > 0))
table(temp$n_b)
temp <- sc %>% group_by(Cell) %>% summarise(n_b = sum(d - a > 1))
table(temp$n_b)

# Of the cells that have coverage on these loci, how many of them have coverage on more than one SNV?
ov2295_dat_ <- ov2295_dat[ov2295_dat$loc %in% loci_with_coverage$loc,]
ov2295_dat_ <- ov2295_dat_[ov2295_dat_$d > 0,]
length(unique(ov2295_dat_$SampleName))
ret <- ov2295_dat_ %>% group_by(SampleName) %>% summarise(n_b = sum(d-a > MIN_READS))
table(ret$n_b)

ret <- subset(ov2295R_dat1, loc %in% loci_with_coverage$loc) %>% group_by(SampleName) %>% summarise(n_b = sum(d-a > MIN_READS))
table(ret$n_b)

ret <- subset(ov2295R_dat2, loc %in% loci_with_coverage$loc) %>% group_by(SampleName) %>% summarise(n_b = sum(d-a > MIN_READS))
table(ret$n_b)

##### Compute basic statistics on coverage and depth. #####
colnames(ov2295_dat)
colnames(ov2295R_dat)
tenx <- rbind(ov2295_dat, ov2295R_dat)
ss3 <- PhylExR::CombineSingleCellReads2(SC_READS_SS3, sc_reads_pattern = "*.txt", file_separator = ",")
ss3$loc <- paste(ss3$chr, ss3$pos, sep=":")

tenx$b <- tenx$d - tenx$a
ss3$b <- ss3$d - ss3$a

# Some cells are overlapping across tumors? Same barcodes? Just coincidence?
sum(unique(ov2295_dat$SampleName) %in% unique(ov2295R_dat$SampleName))
unique(ov2295_dat$SampleName)[which(unique(ov2295_dat$SampleName) %in% unique(ov2295R_dat$SampleName))]
sum(unique(ov2295R_dat1$SampleName) %in% unique(ov2295R_dat2$SampleName))

length(unique(ss3$SampleName))
length(unique(ov2295_dat$SampleName))
length(unique(ov2295R_dat$SampleName))
length(unique(tenx$SampleName))


names(tenx)
names(ss3)
column_names <- c("ID", "Cell", "a", "d", "SampleName", "b")
dat <- rbind(data.frame(tenx[,column_names], type="10X"), data.frame(ss3[,column_names], type="Smart-Seq3"))

# Compare the total and variant depth at a loci given that it is expressed.
pl <- subset(dat, d > 0) %>% ggplot(aes(x = d, fill = type)) + geom_histogram(position = "dodge") + theme_bw() + scale_y_sqrt()
pl <- pl + xlab("Total depth") + ylab("Sqrt of frequencies") + theme(legend.title = element_blank(), legend.text = element_text(size = 2*base_size), legend.position = "top")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
ggsave(pl, filename = "~/phylo-express-paper/figures/Laks/SS3_10X_total_depth.pdf", height = 8, width=8, units = "in")

pl <- subset(dat, b > 0) %>% ggplot(aes(x = b, fill = type)) + geom_histogram(position = "dodge", bins = 30) + theme_bw() + scale_y_sqrt()
pl <- pl + xlab("Variant depth") + ylab("Sqrt of frequencies") + theme(legend.title = element_blank(), legend.text = element_text(size = 2*base_size), legend.position = "top")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
ggsave(pl, filename = "~/phylo-express-paper/figures/Laks/SS3_10X_variant_depth.pdf", height = 8, width=8, units = "in")

dat$SampleName <- gsub(pattern = ".reads.txt", replacement = "", x = dat$SampleName)

### Figure 8ab
write.table(dat, file = "data/NatComm/SupplementaryFigure8a-b.csv", quote = F, row.names = F, col.names = T)


# Compare coverage per cell.
coverage10x <- tenx %>% group_by(SampleName) %>% summarise(n_d = sum(d > 0), n_b = sum(b > 0))
coverageSS3 <- ss3 %>% group_by(SampleName) %>% summarise(n_d = sum(d > 0), n_b = sum(b > 0))

length(unique(tenx$SampleName))
summary(coverage10x$n_b)
summary(coverageSS3$n_b)

table(coverage10x$n_d)
table(coverage10x$n_b)

hist(coverage10x$n_b)
hist(coverageSS3$n_b)

sum(coverage10x$n_b >= 2)
sum(coverageSS3$n_b >= 2)

coverage.df <- rbind(data.frame(coverage10x, type="10X"), data.frame(coverageSS3, type="Smart-Seq3"))
coverage.df$SampleName <- gsub(pattern = ".reads.txt", replacement = "", x = coverage.df$SampleName)
head(coverage.df)
coverage.df %>% group_by(type) %>% summarise(mean=mean(n_b), med=median(n_b), q1=quantile(n_b, 0.25), q3=quantile(n_b, 0.75), max = max(n_b))

pl <- ggplot(coverage.df, aes(y=n_b, x = type)) + geom_boxplot() + theme_bw() + theme(legend.title = element_blank()) + ylab("Co-occurrence of variants (counts)") + xlab("")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(axis.text.x = element_text(size = 2 * base_size))
pl
ggsave(pl, filename = "~/phylo-express-paper/figures/Laks/SS3_10X_co-occurrence.pdf", height = 8, width = 8, units = "in")

### Table 2 and Figure 8c.
write.table(out.df, file = "data/NatComm/Table2SupplementaryFigure8c.csv", sep=",", row.names = F, col.names = T, quote = F)
###

sc_10x <- read.table("~/PhylExAnalysis/data/HGSOC_10X_sc.txt", header=T)
length(unique(sc_10x$SampleName))

# Let's identify loci using stringent filtering.
MIN_READS <- 1
ov2295_coverage <- ov2295_dat %>% group_by(loc) %>% summarise(n_b = sum(d - a > MIN_READS))
ov2295R1_coverage <- ov2295R_dat1 %>% group_by(loc) %>% summarise(n_b = sum(d - a > MIN_READS))
ov2295R2_coverage <- ov2295R_dat2 %>% group_by(loc) %>% summarise(n_b = sum(d - a > MIN_READS))

coverage <- left_join(left_join(ov2295_coverage, ov2295R1_coverage, by = "loc"), ov2295R2_coverage, by = "loc")
coverage.df <- as.data.frame(coverage)
coverage.df$b <- rowSums(coverage.df[,-1])

coverage.df <- left_join(coverage.df, gt_labels)
loci_with_coverage <- subset(coverage.df, b >= 2)
dim(loci_with_coverage)
loci_with_coverage

tenx$cell_id <- paste(tenx$Cell, tenx$SampleName, sep="_")
length(unique(tenx$cell_id))

MIN_READS <- 1
temp <- tenx %>% group_by(cell_id) %>% summarise(coverage = sum(b >= MIN_READS))
table(temp$coverage) # 551 cells if MIN_READS = 1

MIN_READS <- 2
cell_mut_coverage <- tenx %>% group_by(cell_id) %>% summarise(coverage = sum(b >= MIN_READS))
table(cell_mut_coverage$coverage) # 12 cells cover 2 loci if MIN_READS = 2...

MIN_COVERAGE <- 1
snv_with_minimum_coverage <- subset(tenx, cell_id %in% subset(cell_mut_coverage, coverage >= MIN_COVERAGE)$cell_id)
mut_coverage <- snv_with_minimum_coverage %>% group_by(ID) %>% summarise(n_cells = sum(b > 0))
table(mut_coverage$n_cells)

# Let's select loci that have at least two cells covering it.
ids <- mut_coverage[mut_coverage$n_cells >= 2,]$ID
tenx_ <- subset(tenx, ID %in% ids)
tenx_$loc <- paste(tenx_$chr, tenx_$pos, sep=":")
locs <- paste(loci$CHROM, loci$POS, sep=":")
loci_ <- loci[which(locs %in% unique(tenx_$loc)),]

head(loci_)
loci_count <- dim(loci_)[1]
write.table(data.frame(ID=paste("s", 0:(loci_count-1), sep=""), loci_), file = "~/data/cell-line/10X2/loci.txt", col.names = T, row.names = F, quote = F, sep="\t")
