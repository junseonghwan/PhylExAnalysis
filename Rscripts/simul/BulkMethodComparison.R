library(ggplot2)
library(scales)

# Extrac colors in the color blind palette.
method_colors <- c("PhylEx"="#999999", "Canopy"="#E69F00", "B-SCITE"="#56B4E9", "ddClone"="#009E73", "PhyloWGS"="#0072B2", "TSSB"="#CC79A7")

base_size <- 12

CompareBulkMethods <- function(data_path, output_path, output_prefix, plot_title,
                               ww = 8, hh = 8) {
    # Figure 2: Comparison to existing joint analysis approaches.
    if (!dir.exists(output_path))
    {
        dir.create(output_path, recursive = T)
    }

    canopy <- read.csv(paste(data_path, "canopy.csv", sep="/") , header=T)
    canopy$Method <- "Canopy"

    tssb <- read.csv(paste(data_path, "case0.csv", sep="/"), header=T)
    tssb$Method <- "TSSB"

    phylowgs_path <- paste(data_path, "phylowgs.csv", sep="/")
    if (file.exists(phylowgs_path)) {
        phylowgs <- read.csv(phylowgs_path, header=T)
        phylowgs$Method <- "PhyloWGS"
    }

    df <- rbind(canopy, phylowgs, tssb)
    df$Method <- factor(df$Method, levels = c("Canopy", "PhyloWGS", "TSSB"))
    method_colors <- method_colors[names(method_colors) %in% unique(df$Method)]

    p <- ggplot(df, aes(Method, VMeasure, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Method") + ylab("V-measure")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_vmeasure.pdf", sep=""), p, width = ww, height = hh, units = "in")

    p <- ggplot(df, aes(Method, AdjMutualInformation, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Method") + ylab("Adjusted mutual information")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_mutual_information.pdf", sep=""), p, width = ww, height = hh, units = "in")

    p <- ggplot(df, aes(Method, AdjRandIndex, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Method") + ylab("Adjusted rand index")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_rand_index.pdf", sep=""), p, width = ww, height = hh, units = "in")

    p <- ggplot(df, aes(Method, AncestralMetric, fill = Method)) + geom_boxplot()
    p <- p + scale_fill_manual(values=method_colors)
    p <- p + theme_bw() + xlab("Method") + ylab("Ancestral reconstruction mean absolute error")
    p <- p + theme(axis.title.x =element_text(size = base_size * 2), axis.text.x = element_text(size = base_size*2))
    p <- p + theme(axis.title.y =element_text(size = base_size * 2), axis.text.y = element_text(size = base_size*2))
    p <- p + theme(legend.text = element_text(size = base_size*2), legend.title = element_blank())
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(size = base_size*2))
    ggsave(paste(output_path, "/", output_prefix, "_ancestral.pdf", sep=""), p, width = ww, height = hh, units = "in")
    
    return(df)
}

data_path <- "_output/simul/binary/snvs/"
output_prefix <- "bulk_binary_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="/")
CompareBulkMethods(data_path, output_path, output_prefix, "Binary tree")

data_path <- "_output/simul/quadternary_multiregion/snvs/"
output_prefix <- "bulk_quadternary_multiregion_sc"
output_path <- paste("_figures/simulation/", output_prefix, sep="/")
ret <- CompareBulkMethods(data_path, output_path, output_prefix, "Multiregion and multifurcating tree")
# Output source data for Supplementary Figure 1c-d.
write.csv(ret, file = "data/NatComm/SupplementaryFigure1c-d.csv", quote = F, row.names = F)


