library(dplyr)
library(ggplot2)

DEPTHS <- c(2, 5)
DROPOUT_LEVELS <- c(90, 95, 98)
CASES <- 0:3
CELL_COUNT <- c(0, 100, 200, 400)
base_size <- 12

# Read in the results
CombineResults <- function(data_path, with_cn = TRUE)
{
    results.df <- data.frame()
    if (with_cn) {
        CN_PATH <- "cn"
    } else {
        CN_PATH <- "wo_cn"
    }
    for (DEPTH in DEPTHS)
    {
        DEPTH_PATH <- paste(data_path, "/depth", DEPTH, "/", CN_PATH, "/", sep="")

        for (DROPOUT_LEVEL in DROPOUT_LEVELS)
        {
            results_path <- paste(DEPTH_PATH, "/dropout", DROPOUT_LEVEL, sep="")
            for (CASE in CASES)
            {
                case_file <- paste(results_path, "/case", CASE, ".csv", sep="")
                case_results <- read.csv(case_file, header = T)
                dat.df <- data.frame(depth=DEPTH, dropout=DROPOUT_LEVEL, cells=CELL_COUNT[CASE+1], case_results)
                results.df <- rbind(results.df, dat.df)
            }
        }
    }
    return(results.df)
}

data_path <- "~/ScRNACloneEvaluation/data/simul/"
results_cn.df <- CombineResults(data_path, with_cn = TRUE)
results_wo_cn.df <- CombineResults(data_path, with_cn = FALSE)
names(results_cn.df)

results_cn.df$Coverage <- paste("Coverage", 100 - results_cn.df$dropout, sep=" ")
results_cn.df$Coverage <- factor(results_cn.df$Coverage, c("Coverage 2", "Coverage 5", "Coverage 10"))

results_wo_cn.df$Coverage <- paste("Coverage", 100 - results_wo_cn.df$dropout, sep=" ")
results_wo_cn.df$Coverage <- factor(results_wo_cn.df$Coverage, c("Coverage 2", "Coverage 5", "Coverage 10"))

#depth_ <- 2
depth_ <- 5
temp <- subset(results_cn.df, depth == depth_)
pl <- ggplot(temp, aes(x = as.factor(cells), y = VMeasure, fill = as.factor(dropout))) + geom_boxplot() + facet_grid(~ Coverage) + theme_bw()
pl <- pl + theme(legend.position = "none") + xlab("Cells") + ylab("V-measure")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(strip.text = element_text(size=base_size*2))
pl
ggsave(pl, filename = paste("~/phylo-express-paper/figures/simul/coverage_depth", depth_, "_simul_vmeasure.pdf", sep=""), width = 10, height = 8, units = "in")

pl <- ggplot(temp, aes(x = as.factor(cells), y = AncestralMetric, fill = as.factor(dropout))) + geom_boxplot() + facet_grid(~ Coverage) + theme_bw()
pl <- pl + theme(legend.position = "none") + xlab("Cells") + ylab("Anc. recon mean absolute error")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(strip.text = element_text(size=base_size*2))
pl
ggsave(pl, filename = paste("~/phylo-express-paper/figures/simul/coverage_depth", depth_, "_simul_ancestral_reconstruction_error.pdf", sep=""), width = 10, height = 8, units = "in")

pl <- ggplot(temp, aes(x = as.factor(cells), y = AdjRandIndex, fill = as.factor(dropout))) + geom_boxplot() + facet_grid(~ Coverage) + theme_bw()
pl <- pl + theme(legend.position = "none") + xlab("Cells") + ylab("Adj. rand index")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(strip.text = element_text(size=base_size*2))
pl
ggsave(pl, filename = paste("~/phylo-express-paper/figures/simul/coverage_depth", depth_, "_simul_adj_rand_index.pdf", sep=""), width = 10, height = 8, units = "in")

pl <- ggplot(temp, aes(x = as.factor(cells), y = AdjMutualInformation, fill = as.factor(dropout))) + geom_boxplot() + facet_grid(~ Coverage) + theme_bw()
pl <- pl + theme(legend.position = "none") + xlab("Cells") + ylab("Adj. mutual information")
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(strip.text = element_text(size=base_size*2))
pl
ggsave(pl, filename = paste("~/phylo-express-paper/figures/simul/coverage_depth", depth_, "_simul_adj_mutual_information.pdf", sep=""), width = 10, height = 8, units = "in")

# Compute the variant coverage.
GetCoverageData <- function(SIMUL_DATA_PATH, CASE)
{
    variant_coverage <- data.frame()
    for (rep in 0:19) {
        REP_PATH <- paste(SIMUL_DATA_PATH, "/rep", rep, sep="")
        CASE_PATH <- paste(REP_PATH, "/case", CASE, sep="")
        sc_data <- read.table(paste(CASE_PATH, "simul_sc.txt", sep="/"), header=T)
        sc_data$b <- sc_data$d - sc_data$a
        variant_coverage_by_cell <- sc_data %>% group_by(Cell) %>% summarise(n_d = sum(d > 0), n_b = sum(b > 0))
        tbl.df <- data.frame(table(variant_coverage_by_cell$n_b))
        variant_coverage <- rbind(variant_coverage, tbl.df)
    }
    return(variant_coverage)
}

ret.df <- data.frame()
for (dropout in DROPOUT_LEVELS)
{
    for (case in 1:3) {
        SIMUL_DATA_PATH <- paste("~/data/simul/10X/depth2/cn/dropout", dropout, sep="")
        ret.df <- rbind(ret.df, data.frame(GetCoverageData(SIMUL_DATA_PATH, case), case=paste(CELL_COUNT[case+1], "cells"), coverage=paste("Coverage", 100 - dropout)))
    }
}
head(ret.df)

ret.df$coverage <- factor(ret.df$coverage, levels = c("Coverage 2", "Coverage 5", "Coverage 10"))
pl <- ggplot(ret.df, aes(Var1, Freq)) + geom_boxplot() + theme_bw() + facet_grid(case ~ coverage)
pl <- pl + theme(axis.text.x = element_text(size = 2*base_size), axis.text.y = element_text(size = 2*base_size), axis.title.x = element_text(size = 2*base_size), axis.title.y = element_text(size = 2*base_size))
pl <- pl + theme(strip.text = element_text(size=base_size*2)) + xlab("Count of mutation coverage") + ylab("Frequency")
pl
ggsave(pl, filename = "~/phylo-express-paper/figures/simul/depth_coverage.pdf", width = 8, height = 8, units = "in")
