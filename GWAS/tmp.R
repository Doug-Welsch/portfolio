# 3/14/2025
# Testing
# Doug Welsch

setwd("~/RProjects/portfolio/GWAS/")

# https://cloufield.github.io/gwaslab/tutorial/
# "The dataset we will use as an example is the sumstats of type 2 diabetes from BBJ (K. Suzuki et al., Nature Genetics. 51, 379–386 (2019).)"
snp.df <- read.table(gzfile("./data/t2d_bbj.txt.gz"), header = T)

# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("MungeSumstats")

# Plotting
library(qqman)
library(ggplot2)
set.seed(123)

generateColors <- function(n, seed = NULL) {
  # Generate evenly spaced HCL hues
  hues <- seq(0, 360, length.out = n + 1)[-1]
  colors <- hcl(h = hues, l = 65, c = 100)
  # Shuffle colors
  if (!is.null(seed)) {
    set.seed(seed)
    colors <- sample(colors, n)
  }
  return(colors)
}

# Pre-processing
qqman.df <- data.frame(
  SNP = snp.df$SNP,
  CHR = snp.df$CHR,
  BP = snp.df$POS,
  P = snp.df$P
)
# Filtering to autosomal chromosomes
qqman.df <- qqman.df[qqman.df$CHR != "X", ]
qqman.df$CHR <- as.numeric(qqman.df$CHR)
# From hardware limitations, subsetting to 50,000 SNPs
qqman.df <- qqman.df[sample(1:nrow(qqman.df), 50000, replace=F), ]

## Manhattan Plot
# TODO(dawelsch): Fix text size
# TODO(dawelsch): ??? Assumes chromosome lengths from Hsapiens?
pdf(file = "./outs/manhattan/Manhattan_sample.pdf", width = 24, height = 8)
manhattan(qqman.df,
          main = "Manhattan Plot",
          cex = 0.6,
          cex.axis = 0.7, # 0.9
          col = c("orange", "blue"),
          suggestiveline = F,
          genomewideline = F,
          annotatePval = 1e-10,
          annotateTop = F,
          chrlabs = as.character(1:22))
dev.off()

## Q-Q Plot
pdf(file = "./outs/QQ/QQ_sample.pdf", width = 8, height = 8)
qq(qqman.df$P,
   main = "Q-Q plot of GWAS p-values",
   xlim = c(0, 7),
   ylim = c(0, 7),
   pch = 18,
   col = "blue4",
   cex = 1.5,
   las = 1)
dev.off()

# Annotate












