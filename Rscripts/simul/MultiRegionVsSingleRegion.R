library(ggplot2)
base_size <- 12

color_scheme <- c("PhylEx \nsingle \nregion\n"="#999999", "Canopy \nmulti \nregion\n"="#E69F00", "TSSB \nmulti \nregion\n"="#CC79A7", "PhyloWGS \nmulti \nregion\n"="#0072B2")

MakeMultiRegionVsSingleRegionComparisonPlots <- function(data_path, output_path, output_prefix, plot_title, method_colors=color_scheme) {
    
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = T)
    }
    # Figure 2: Comparison to existing joint analysis approaches.
    canopy <- read.csv(paste(data_path, "canopy.csv", sep="/") , header=T)
    canopy$Method <- "Canopy \nmulti \nregion\n"
    canopy$CellCount <- 0

    phylowgs_path <- paste(data_path, "phylowgs.csv", sep="/")
    if (file.exists(phylowgs_path)) {
        phylowgs <- read.csv(phylowgs_path, header=T)
        phylowgs$Method <- "PhyloWGS \nmulti \nregion\n"
    } else {
        phylowgs <- read.csv(paste(data_path, "case0.csv", sep="/"), header=T)
        phylowgs$Method <- "TSSB \nmulti \nregion\n"
    }
    phylowgs$CellCount <- 0

    df1 <- read.csv(paste(data_path, "single_region_case1.csv", sep="/"), header=T)
    df2 <- read.csv(paste(data_path, "single_region_case2.csv", sep="/"), header=T)
    df3 <- read.csv(paste(data_path, "single_region_case3.csv", sep="/"), header=T)

    df1$CellCount <- 100
    df2$CellCount <- 200
    df3$CellCount <- 400

    our_df <- rbind(df1, df2, df3)
    our_df$Method <- "PhylEx \nsingle \nregion\n"

    df <- rbind(canopy, phylowgs, our_df)

    df$Method <- factor(df$Method, levels = c("PhylEx \nsingle \nregion\n", "Canopy \nmulti \nregion\n", "PhyloWGS \nmulti \nregion\n", "TSSB \nmulti \nregion\n"))
    method_colors <- method_colors[names(method_colors) %in% unique(df$Method)]

    p <- ggplot(df, aes(as.factor(CellCount), VMeasure, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Cell count") + ylab("V-measure")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*1.6), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_vmeasure.pdf", sep=""), p, width = 8, height = 8, units = "in")

    p <- ggplot(df, aes(as.factor(CellCount), AdjMutualInformation, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Adjusted mutual information")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*1.6), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_mutual_information.pdf", sep=""), p, width = 8, height = 8, units = "in")

    p <- ggplot(df, aes(as.factor(CellCount), AdjRandIndex, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Adjusted rand index")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*1.6), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_rand_index.pdf", sep=""), p, width = 8, height = 8, units = "in")

    df_a <- rbind(canopy, phylowgs, our_df)
    df_a$Method <- factor(df_a$Method, levels = c("PhylEx \nsingle \nregion\n", "Canopy \nmulti \nregion\n", "TSSB \nmulti \nregion\n", "PhyloWGS \nmulti \nregion\n"))
    p <- ggplot(df_a, aes(as.factor(CellCount), AncestralMetric, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Cell count") + ylab("Ancestral reconstruction mean absolute error")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*1.6), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    p <- p + guides(colour = guide_legend(nrow = 2))
    ggsave(paste(output_path, "/", output_prefix, "_ancestral.pdf", sep=""), p, width = 8, height = 8, units = "in")
    
    return(df)
}

data_path <- "_output/simul/quadternary_cn_multiregion/snvs/"
output_prefix <- "quadternary_cn_multiregion_vs_single_region"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeMultiRegionVsSingleRegionComparisonPlots(data_path, output_path, output_prefix, "Multiregion and multifurcating tree with CNV")
ret$Method <- gsub(pattern = "\n", replacement = "", ret$Method)
write.csv(ret, file = "data/NatComm/Figure2c-d.csv", row.names = F, quote = F)

data_path <- "_output/simul/quadternary_multiregion/snvs/"
output_prefix <- "quadternary_multiregion_vs_single_region"
output_path <- paste("_figures/simulation/", output_prefix, sep="")
ret <- MakeMultiRegionVsSingleRegionComparisonPlots(data_path, output_path, output_prefix, "Multiregion and multifurcating tree")
