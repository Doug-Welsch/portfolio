# 3/16/2025
# Computing Major Cell Type DEG, Tumor vs. Healthy
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0

mm.sobj <- readRDS("./outs/SeuratObject/SeuratObj_init.rds")
Idents(mm.sobj) <- mm.sobj$Major.cell.type
# Healthy cell Seurat object
h.sobj <- subset(mm.sobj, subset = (Tumor.or.healthy == "h"))
print(table(Idents(h.sobj), h.sobj$Library.prep.batch))

# Compute DEG
# TODO(dawelsch): Batching by library prep batch, using FindConservedMarkers(),
# can alternatively use FindMarkers() with test LR and latent variable of batch

categories <- unique(mm.sobj@meta.data$Major.cell.type)
for (cat in categories) {
  print(cat)
  
  # Major Cell Types
  # TODO(dawelsch): These should be consistent with other literature when
  # considering only healthy cells
  print("Major Cell Type")
  try({
    markers <- FindConservedMarkers(h.sobj, ident.1 = cat, grouping.var = "Library.prep.batch")
    markers <- cbind(gene = rownames(markers), markers)
    write.csv(markers, paste0("./outs/DEG/Markers_Healthy_MajorCellType_", cat, ".csv"), row.names = F)
  })
  
  # Tumor vs. Healthy
  print("Tumor vs. Healthy")
  try({
    cat.sobj <- subset(mm.sobj, subset = (Major.cell.type == cat))
    Idents(cat.sobj) <- cat.sobj@meta.data$Tumor.or.healthy
    markers <- FindConservedMarkers(cat.sobj, ident.1 = "t", ident.2 = "h", grouping.var = "Library.prep.batch")
    markers <- cbind(gene = rownames(markers), markers)
    write.csv(markers, paste0("./outs/DEG/Markers_MajorCellType_Tumor-Healthy_", cat, ".csv"), row.names = F)
  })
}

# Summarize DEG
# TODO(dawelsch): Typically, top genes are ordered by fold change

marker.files <- list.files("./outs/DEG/", pattern = "Markers_Healthy_MajorCellType_*", full.names = T)
markers.ls <- list()
for (fn in marker.files) {
  cat <- substring(basename(fn), 31, nchar(basename(fn)) - 4)
  print(cat)
  
  markers <- read.csv(fn)
  # Subset to positive (upregulated) across batches
  try({markers <- markers[markers$round1_20151128_avg_log2FC > 0, ]}, T)
  try({markers <- markers[markers$round2_20151217_avg_log2FC > 0, ]}, T)
  try({markers <- markers[markers$round3_20160313_avg_log2FC > 0, ]}, T)
  # Get top markers
  is.empty <- T
  if ("max_pval" %in% colnames(markers) &
      "minimump_p_val" %in% colnames(markers)) {
    # TODO(dawelsch): Tentative interpretation of FindConservedMarkers()
    markers <- markers[markers$max_pval < 0.05, ]
    if (nrow(markers) > 0) {
      markers <- markers[order(markers$minimump_p_val), ][1:min(10, nrow(markers)), ]
      is.empty <- F
    }
  } else {
    cn <- colnames(markers)[grepl('p_val_adj', colnames(markers), fixed=T)]
    markers <- markers[markers[[cn]] < 0.05, ]
    if (nrow(markers) > 0) {
      cn <- colnames(markers)[grepl('FC', colnames(markers), fixed=T)]
      markers <- markers[order(markers[[cn]]), ][1:min(10, nrow(markers)), ]
      is.empty <- F
    }
  }
  if (!is.empty) {
    markers$Major.cell.type <- cat
    markers.ls[[cat]] <- markers
  } else {
    print("No significant positive genes")
  }
}
markers.summary <- dplyr::bind_rows(markers.ls)
write.csv(markers.summary, "./outs/DEG/Summary_Healthy_MajorCellType.csv", row.names = F)


