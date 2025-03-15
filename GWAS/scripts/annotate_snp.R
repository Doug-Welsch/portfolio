# 3/15/2025
# Annotating signficant SNPs with effect
# Doug Welsch

# Environment Configuration

setwd("~/RProjects/portfolio/GWAS/")

library(stats)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(openxlsx)

# Pre-processing Data

# See: https://cloufield.github.io/gwaslab/tutorial/
# "The dataset we will use as an example is the sumstats of type 2 diabetes from BBJ (K. Suzuki et al., Nature Genetics. 51, 379–386 (2019).)"
# https://doi.org/10.1038/s41588-018-0332-4
snp.df <- read.table(gzfile("./data/t2d_bbj.txt.gz"), header = T)

snp.sig.df <- snp.df
snp.sig.df$P.Adj <- p.adjust(snp.sig.df$P, method = "bonferroni")
snp.sig.df <- snp.sig.df[snp.sig.df$P.Adj < 5e-8, ]

# Annotate

# Prepare regions
snp.sig.gr <- GRanges(
  seqnames = paste0("chr", snp.sig.df$CHR),
  ranges = IRanges(snp.sig.df$POS, snp.sig.df$POS, names = snp.sig.df$SNP))
# Select annotations for intersection with regions
annots <- c('hg19_basicgenes', 'hg19_genes_intergenic', 'hg19_genes_intronexonboundaries')
# Build the annotations (a single GRanges object)
annotations <- build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we with the annotations
annotation.gr <- annotate_regions(
  regions = snp.sig.gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# Append annotations
# TODO(dawelsch): Merging can likely be optimized using tidyverse
annotation.df <- data.frame(annotation.gr)
write.csv(annotation.df, "./outs/annotation/Annotation_GRanges.csv", row.names = F)
annotation.df$annot.KEY <- paste0(gsub("chr", "", annotation.df$seqnames), ":", annotation.df$start)
annotation.df <- annotation.df[, colnames(annotation.df)[startsWith(colnames(annotation.df), "annot.")]]
colnames(annotation.df) <- gsub("annot.", "Annotation: ", colnames(annotation.df))
colnames(annotation.df)[colnames(annotation.df) == "Annotation: seqnames"] <- "Annotation: CHR"
snp.sig.df[["Annotation: KEY"]] <- paste0(snp.sig.df$CHR, ":", snp.sig.df$POS)
annotation.df <- merge(snp.sig.df, annotation.df, by = "Annotation: KEY")
annotation.df[["Annotation: KEY"]] <- NULL

# Save Outputs

# CSV
write.csv(annotation.df, "./outs/annotation/Annotation_SigSNP.csv", row.names = F)
# XLSX
# TODO(dawelsch): Fix number of digits for select columns
wb <- createWorkbook()
addWorksheet(wb, "Annotations")
writeDataTable(wb, "Annotations", annotation.df)
saveWorkbook(wb, "./outs/annotation/Annotation_SigSNP.xlsx", overwrite = TRUE)











