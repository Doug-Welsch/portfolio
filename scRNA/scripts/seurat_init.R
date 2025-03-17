# 3/15/2025
# Testing
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0
library(data.table)
library(tidyverse)
# library(scCustomize)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465
meta.df <- fread("./data/GSE127465_mouse_cell_metadata_15939x12.tsv.gz")
norm.mat <- Matrix::readMM(gzfile("./data/GSE127465_mouse_counts_normalized_15939x28205.mtx.gz"))

# Parse raw count files
raw.files <- list.files("./data/GSE127465/", full.names = T)
cnts.ls <- list()
for (fn in raw.files) {
  sample.id <- substr(basename(fn), 18, 22)
  print(sample.id)
  cnts.df <- fread(fn)
  cnts.df$barcode <- paste0(sample.id, "-", cnts.df$barcode)
  cnts.df <- cnts.df %>% column_to_rownames("barcode")
  cnts.ls[[sample.id]] <- Matrix::Matrix(as.matrix(t(cnts.df)), sparse = T)
}
cnts.mat <- do.call("cbind", cnts.ls)
rm(cnts.ls, cnts.df); gc()
# Subset to processed cells with meta information
meta.df$KEY <- paste0(meta.df$Library, "-", meta.df$Barcode)
cnts.mat <- cnts.mat[, meta.df$KEY]
mm.sobj <- CreateSeuratObject(counts = cnts.mat)
# Add meta information
all.equal(Cells(mm.sobj), meta.df$KEY)
for (cn in colnames(meta.df)) {
  mm.sobj <- AddMetaData(
    object = mm.sobj,
    metadata = meta.df[[cn]],
    col.name = cn
  )
}
mm.sobj$KEY <- NULL
# Add processed counts
mm.sobj[["RNA"]]@data <- as(t(norm.mat), "CsparseMatrix")
rownames(mm.sobj@assays$RNA@data) <- rownames(mm.sobj@assays$RNA@counts)
colnames(mm.sobj@assays$RNA@data) <- colnames(mm.sobj@assays$RNA@counts)
# Add SPRING
spring.df <- mm.sobj@meta.data[, c("x", "y")]
colnames(spring.df) <- c("SPRING_1", "SPRING_2")
mm.sobj[["spring"]] <- CreateDimReducObject(embeddings = as.matrix(spring.df),
                                            key = "SPRING_",
                                            assay = "RNA")
# Save
saveRDS(mm.sobj, "./outs/SeuratObject/SeuratObj_init.rds")




















