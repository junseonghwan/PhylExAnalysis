
# df: Data frame with columns W1, W2, Node, and Clone.
#' @import ggplot2
ReducedDimensionPlots <- function(df, base_size = 12) {
  p <- ggplot(df, aes(W1, W2, colour=Node)) + geom_point(alpha=0.8) + theme_classic()
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size*2))
  p <- p + theme(legend.text = element_text(size=base_size*2), legend.title = element_text(size=base_size*2))
  p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  p <- p + xlab("Dim1") + ylab("Dim2")
  return(p)
}

# Assume volcano.df is already ordered based on significance
#' @import ggplot2
MakeVolcanoPlot <- function(volcano.df, plot_title, num_genes_to_label = 10, base_size = 12) {
  p <- ggplot(volcano.df, aes(logfc, logpvalue, col = significant)) + geom_point() + theme_bw() + ylab("-Log p-value") + xlab("Log fold change")
  p <- p + labs(title=plot_title)
  p <- p + scale_color_manual(breaks=c("FALSE", "TRUE"), values=cols)
  p <- p + theme(axis.title.x=element_text(size = base_size * 2), axis.title.y=element_text(size = base_size * 2))
  p <- p + theme(legend.position="bottom") + guides(color=guide_legend(title="Significantly differential?"))
  p <- p + theme(legend.text = element_text(size=base_size*2), title = element_text(size=base_size*2))
  p <- p + theme(axis.text.x = element_text(size = base_size * 2), axis.text.y = element_text(size = base_size * 2))
  # Choose genes to label based on FDR.
  annotate.df <- volcano.df[1:num_genes_to_label,]
  temp.df <- volcano.df[order(volcano.df$logfc),]
  # Choose genes to label based on log fold change.
  annotate.df2 <- temp.df[1:num_genes_to_label,]
  annotate.df <- rbind(annotate.df, annotate.df2)
  annotate.df <- unique(annotate.df)
  p <- p + ggrepel::geom_label_repel(data = annotate.df, aes(label=gene_name),
                                     size = 8,
                                     nudge_x = 0.2,
                                     col = "black",
                                     segment.size  = 1,
                                     segment.color = "grey50",
                                     direction     = "both", max.overlaps = 10)
  return(p)
}

#' @import ggplot2
PlotSigGenes <- function(dge_results, mart, gene_plot_count = 50) {
  # Plot the significantly differentially expressed genes along with their p-values.
  ordering <- order(dge_results$table$PValue)
  tbl <- dge_results$table[ordering,]
  tbl$ensembl_gene_id <- rownames(tbl)
  
  bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id",
                           'start_position', 'end_position', 'percentage_gene_gc_content'),
              filters = 'ensembl_gene_id',
              values = tbl$ensembl_gene_id,
              mart = mart, uniqueRows = T)
  
  ret <- dplyr::left_join(tbl, bm, by = "ensembl_gene_id")
  ret <- ret[!is.na(ret$hgnc_symbol),]
  ret.df <- data.frame(gene = ret$hgnc_symbol[1:gene_plot_count], logpvalue = -log(ret$PValue)[1:gene_plot_count], logFC=ret$logFC[1:gene_plot_count] > 0)
  ret.df$gene <- factor(ret.df$gene, levels = ret.df$gene)
  ret.df$logFC[ret.df$logFC == FALSE] <- "DOWN"
  ret.df$logFC[ret.df$logFC == TRUE] <- "UP"
  p <- ggplot(ret.df, aes(gene, logpvalue, col=logFC)) + geom_point() + theme_bw() + xlab("Gene") + ylab("-Log p-value")
  p <- p + scale_color_manual(values=c("DOWN"="blue", "UP"="red"))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title = element_blank())
  return(p)
}

PlotTotalCounts <- function(sc,  base_size = 12) {
  sc$b <- sc$d - sc$a
  p <- ggplot(sc, aes(ID, Cell, fill = log(d + 1))) + geom_tile(colour = "white")
  p <- p + theme_bw() + scale_fill_gradient(low = "white", high = "red")
  p <- p + ylab("Cell") + xlab("Loci")
  p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
  p <- p + theme(axis.ticks = element_blank())
  p <- p + theme(axis.title.x =element_text(size = base_size * 2))
  p <- p + theme(axis.title.y =element_text(size = base_size * 2))
  p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

CoclusteringPlot <- function(sc, cell.df, base_size = 12) {
  sc$b <- sc$d - sc$a
  sc_join <- left_join(sc, cell.df, by = "Cell")
  cells_clustered <- cell.df[order(cell.df$Node), "Cell"]
  sc_join$Cell <- factor(sc_join$Cell, levels = cells_clustered)
  
  mut_ids_clustered <- datum2node_[order(datum2node_$Node),"ID"]
  sc_join$ID <- factor(sc_join$ID, levels = mut_ids_clustered)
  names(sc_join)
  
  p <- ggplot(subset(sc_join, b >0), aes(ID, Cell, fill=Node)) + geom_tile(colour = "white")
  p <- p + theme_bw()
  p <- p + xlab("Loci") + ylab("Cell")
  p <- p + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
  p <- p + theme(axis.ticks = element_blank())
  p <- p + theme(axis.title.x =element_text(size = base_size * 2))
  p <- p + theme(axis.title.y =element_text(size = base_size * 2))
  p <- p + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + theme(legend.text = element_text(size=7, face="bold")) + guides(fill=guide_legend(title="Clone", size = 7))
  return(p)
}
