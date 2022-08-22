library(ggplot2)
#setwd("~/phylo-express-paper/figures/")

epsilon <- 0.01
chi <- rbeta(n = 10000, shape1 = 100 * epsilon, shape2 = 100 * (1 - epsilon))
hist(chi)

PlotFunction <- function(p, base_size = 12) {
  p <- p + theme(axis.title = element_text(size = base_size * 2))
  p <- p + theme(axis.text.x = element_text(size = base_size * 2), axis.text.y = element_text(size = base_size * 2))
  p <- p + theme(plot.title = element_text(size = base_size * 2))
  p <- p + theme(legend.text = element_text(size = base_size * 2), legend.title = element_text(size = base_size*1.5))
  return(p)
}

PlotMixture <- function(alpha0, beta0, alpha_n, beta_n, file_path, base_size = 12) {
  x_mono <- rbeta(n, shape1 = alpha0, shape2 = beta0)
  x_bi <- rbeta(n, shape1 = alpha_n, shape2 = beta_n)

  x_mono.df <- data.frame(y=x_mono, type="Mono-allelic")
  x_bi.df <- data.frame(y=x_bi, type="Bi-allelic")
  df <- rbind(x_mono.df, x_bi.df)
  df$type <- factor(df$type, levels = c("Mono-allelic", "Bi-allelic"))
  p <- ggplot(df, aes(y, fill = type)) + geom_histogram(position = "dodge")
  p <- p + theme_bw()
  p <- p + scale_y_continuous(expand = c(0, 0), breaks = seq(0, 5000, 500), limits = c(0, 5200))
  p <- p + ylab("") + xlab("") + labs(fill = "Distribution")
  p <- p + theme(axis.text.x = element_text(size = base_size*2), axis.text.y = element_text(size = base_size * 2))
  p <- p + theme(legend.text = element_text(size = base_size * 2), legend.title = element_text(size = base_size*1.5))
  ggsave(filename = file_path, p, width = 8, height = 8, units = "in")
  return(list(pl=p, dat=df))
}

set.seed(123)

n <- 10000
alpha_n <- 5
beta_n <- 7

alpha0 <- 0.05
beta0 <- 0.05
ret <- PlotMixture(alpha0, beta0, alpha_n, beta_n, "_figures/supp/scRNA_05.pdf")
write.csv(ret$dat, file = "data/NatComm/SupplementaryFigure7a.csv", quote = F, row.names = F)

alpha0 <- 0.01
beta0 <- 0.01
ret <- PlotMixture(alpha0, beta0, alpha_n, beta_n, "_figures/supp/scRNA_01.pdf")
write.csv(ret$dat, file = "data/NatComm/SupplementaryFigure7b.csv", quote = F, row.names = F)

alpha0 <- 0.005
beta0 <- 0.005
pl <- PlotMixture(alpha0, beta0, alpha_n, beta_n, "_figures/supp/scRNA_005.pdf")

alpha0 <- 0.001
beta0 <- 0.001
pl <- PlotMixture(alpha0, beta0, alpha_n, beta_n, "_figures/supp/scRNA_001.pdf")

#alpha0 <- 0.01
#beta0 <- 0.01
#alpha_n <- 3
#beta_n <- 10
#pl <- PlotMixture(alpha0, beta0, alpha_n, beta_n, "_figures/supp/scRNA_01_.pdf")


df <- data.frame()
epsilons <- c(0.05, 0.01, 0.005)
for (epsilon in epsilons) {
  xx <- rbeta(n, shape1 = epsilon, shape2 = 1 - epsilon)
  df <- rbind(df, data.frame(y=xx, type=paste("", epsilon, sep="")))
}
df$type <- factor(df$type, levels = c("0.05", "0.01", "0.005"))

p <- ggplot(df, aes(y, fill = type)) + geom_histogram(position = "dodge")
p <- p + theme_bw()
p <- p + scale_y_continuous(expand = c(0, 0))
p <- p + xlab("") + ylab("") + labs(fill = "Sequencing error")
p <- PlotFunction(p)
ggsave(filename = "_figures/supp/scRNA_error.pdf", p, width = 8, height = 8, units = "in")
write.csv(df, file = "data/NatComm/SupplementaryFigure7c.csv", quote = F, row.names = F)
