library(ggplot2)
library(scales)

#method_colors <- c("PhylEx"="#e41a1c", "Canopy"="#377eb8", "B-SCITE"="#4daf4a", "ddClone"="#984ea3", "PhyloWGS"="#ff7f00")
# Extrac colors in the color blind palette.
#, , , "#0072B2", "#D55E00"
method_colors <- c("PhylEx"="#999999", "Canopy"="#E69F00", "B-SCITE"="#56B4E9", "ddClone"="#009E73", "PhyloWGS"="#0072B2", "TSSB"="#CC79A7")

base_size <- 12

MakeSCComparisonPlots <- function(data_path, output_path, output_prefix, plot_title,
                                  multiregion = FALSE, include_sc = TRUE,
                                  ww = 8, hh = 8) {
    # Figure 2: Comparison to existing joint analysis approaches.
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = T)
    }

    canopy <- read.csv(paste(data_path, "canopy.csv", sep="/") , header=T)
    canopy$Method <- "Canopy"
    canopy$CellCount <- 0

    phylowgs_path <- paste(data_path, "phylowgs.csv", sep="/")
    if (file.exists(phylowgs_path)) {
        phylowgs <- read.csv(phylowgs_path, header=T)
        phylowgs$Method <- "PhyloWGS"
    } else {
        phylowgs <- read.csv(paste(data_path, "case0.csv", sep="/"), header=T)
        phylowgs$Method <- "TSSB"
    }
    phylowgs$CellCount <- 0

    df1 <- read.csv(paste(data_path, "case1.csv", sep="/"), header=T)
    df2 <- read.csv(paste(data_path, "case2.csv", sep="/"), header=T)
    df3 <- read.csv(paste(data_path, "case3.csv", sep="/"), header=T)

    df1$Method <- "PhylEx"
    df2$Method <- "PhylEx"
    df3$Method <- "PhylEx"
    df1$CellCount <- 100
    df2$CellCount <- 200
    df3$CellCount <- 400

    bscite_df1 <- read.csv(paste(data_path, "bscite_case1.csv", sep="/"), header=T)
    bscite_df2 <- read.csv(paste(data_path, "bscite_case2.csv", sep="/"), header=T)
    bscite_df3 <- read.csv(paste(data_path, "bscite_case3.csv", sep="/"), header=T)

    bscite_df1$Method <- "B-SCITE"
    bscite_df2$Method <- "B-SCITE"
    bscite_df3$Method <- "B-SCITE"
    bscite_df1$CellCount <- 100
    bscite_df2$CellCount <- 200
    bscite_df3$CellCount <- 400

    bscite_df <- rbind(bscite_df1, bscite_df2, bscite_df3)
    our_df <- rbind(df1, df2, df3)

    if (include_sc) {
        if (!multiregion) {
            ddclone_df1 <- read.csv(paste(data_path, "ddClone_case1.csv", sep="/"), header=T)
            ddclone_df2 <- read.csv(paste(data_path, "ddClone_case2.csv", sep="/"), header=T)
            ddclone_df3 <- read.csv(paste(data_path, "ddClone_case3.csv", sep="/"), header=T)

            ddclone_df1$AncestralMetric <- NA
            ddclone_df2$AncestralMetric <- NA
            ddclone_df3$AncestralMetric <- NA

            ddclone_df1$CellCount <- 100
            ddclone_df2$CellCount <- 200
            ddclone_df3$CellCount <- 400

            ddclone_df1$Method <- "ddClone"
            ddclone_df2$Method <- "ddClone"
            ddclone_df3$Method <- "ddClone"

            ddclone_df <- rbind(ddclone_df1, ddclone_df2, ddclone_df3)
            df <- rbind(canopy, phylowgs, our_df, bscite_df, ddclone_df)
        } else {
            df <- rbind(canopy, phylowgs, our_df, bscite_df)
        }
    } else {
        df <- rbind(canopy, phylowgs, our_df)
    }

    df$Method <- factor(df$Method, levels = c("PhylEx", "Canopy", "B-SCITE", "ddClone", "PhyloWGS", "TSSB"))
    method_colors_ <- method_colors[names(method_colors) %in% unique(df$Method)]

    p <- ggplot(df, aes(as.factor(CellCount), VMeasure, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors_)
    p <- p + theme_bw() + xlab("Cell count") + ylab("V-measure")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_vmeasure.pdf", sep=""), p, width = ww, height = hh, units = "in")

    p <- ggplot(df, aes(as.factor(CellCount), AdjMutualInformation, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors_)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Adjusted mutual information")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_mutual_information.pdf", sep=""), p, width = ww, height = hh, units = "in")

    p <- ggplot(df, aes(as.factor(CellCount), AdjRandIndex, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors_)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Adjusted rand index")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_rand_index.pdf", sep=""), p, width = ww, height = hh, units = "in")

    df_a <- rbind(canopy, phylowgs, our_df, bscite_df)
    df_a$Method <- factor(df_a$Method, levels = c("PhylEx", "Canopy", "B-SCITE", "ddClone", "PhyloWGS", "TSSB"))
    method_colors_ <- method_colors[names(method_colors) %in% unique(df_a$Method)]
    p <- ggplot(df_a, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors_)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Ancestral reconstruction mean absolute error")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_ancestral.pdf", sep=""), p, width = ww, height = hh, units = "in")
    
    return(list(cluster.df=df, ancestral_error.df=df_a))
}

data_path <- "_output//simul/binary/snvs/"
output_prefix <- "binary_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
MakeSCComparisonPlots(data_path, output_path, output_prefix, "Binary tree")

data_path <- "_output/simul/binary_cn/snvs/"
output_prefix <- "binary_cn_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeSCComparisonPlots(data_path, output_path, output_prefix, "Binary tree with CNV")
# Output source data for Supplementary Figure 1a-b.
write.csv(ret$cluster.df, file = "data/NatComm/SupplementaryFigure1a-b.csv", quote = F, row.names = F)

data_path <- "_output/simul/quadternary_cn/snvs/"
output_prefix <- "quadternary_cn_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeSCComparisonPlots(data_path, output_path, output_prefix, "Multifurcating tree with CNV")
# Output source data for Figure 2a-b.
write.csv(ret$cluster.df, file = "data/NatComm/Figure2a-b.csv", quote = F, row.names = F)

data_path <- "_output/simul/quadternary_cn_multiregion/snvs/"
output_prefix <- "quadternary_cn_multiregion_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeSCComparisonPlots(data_path, output_path, output_prefix, "Multiregion and multifurcating tree with CNV", multiregion = TRUE)
# Output source data for Figure 2e-f.
write.csv(ret$cluster.df, file = "data/NatComm/Figure2e-f.csv", quote = F, row.names = F)

data_path <- "_output/simul/quadternary_multiregion/snvs/"
output_prefix <- "quadternary_multiregion_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeSCComparisonPlots(data_path, output_path, output_prefix, "Multiregion and multifurcating tree", multiregion = TRUE, include_sc = TRUE, ww = 8, hh = 8)
# Output source data for Supplementary Figure 1e-f.
write.csv(ret$cluster.df, file = "data/NatComm/SupplementaryFigure1e-f.csv", quote = F, row.names = F)


