# 3/15/2025
# Testing
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # 4.3.0
library(data.table)
library(tidyverse)
# library(scCustomize)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465
meta.df <- fread("./data/GSE127465_mouse_cell_metadata_15939x12.tsv.gz")

raw.files <- list.files("./data/GSE127465/", full.names = T)
# sobj.ls <- list()
cnts.ls <- list()
for (fn in raw.files) {
  sample.id <- substr(basename(fn), 18, 22)
  sample.id.vec <- c(sample.id.vec, sample.id)
  print(sample.id)
  cnts.df <- fread(fn)
  cnts.df$barcode <- paste0(sample.id, "-", cnts.df$barcode)
  cnts.ls[[sample.id]] <- cnts.df
}
cnts.df <- do.call("rbind", cnts.ls)
mm.sobj <- CreateSeuratObject(counts = cnts.df)

# foo <- Matrix::readMM("./data/GSE127465_mouse_counts_normalized_15939x28205.mtx.gz")
# bar <- as.data.frame(foo)
# baz <- CreateSeuratObject(counts = bar)























