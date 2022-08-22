library(ggplot2)
library(ggpubr)
library(dplyr)

binary_loss <- read.table("_output/simul/binary/loss.tsv", header = T)
binary_cn_loss <- read.table("_output/simul/binary_cn/loss.tsv", header = T)

names(binary_loss)
pl <- ggplot(binary_loss, aes(as.factor(cells), values, fill = type)) + geom_boxplot() + theme_bw()
pl <- pl + xlab("Cells") + ylab("Expected absolute loss") +  theme(legend.title=element_blank()) + ggtitle("Binary tree w/o copy number evolution")
pl <- pl + theme(legend.text = element_text(size = 16), axis.title = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) 
pl

pl2 <- ggplot(binary_cn_loss, aes(as.factor(cells), values, fill = type)) + geom_boxplot() + theme_bw()
pl2 <- pl2 + xlab("Cells") + ylab("Expected absolute loss") + theme(legend.title=element_blank()) + ggtitle("Binary tree with copy number evolution")
pl2 <- pl2 + theme(legend.text = element_text(size = 16), axis.title = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) 
pl2

pl_combined <- ggarrange(pl, pl2, labels = c("a", "b"), common.legend = TRUE, legend = "bottom", font.label=list(color="black",size=16))
ggsave(pl_combined, filename = "_figures/simulation/two_step_comparison.pdf")

# Save the source data.
write.csv(binary_loss, file = "data/NatComm/SupplementaryFigure2a.csv", quote = F, row.names = F)
write.csv(binary_cn_loss, file = "data/NatComm/SupplementaryFigure2b.csv", quote = F, row.names = F)
