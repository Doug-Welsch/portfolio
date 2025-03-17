# 3/16/2025
# Computing Major Cell Type DEG
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0

mm.sobj <- readRDS("./outs/SeuratObject/SeuratObj_init.rds")

# Compute DEG
# TODO(dawelsch): Batching by library prep batch, using FindConservedMarkers(),
# can alternatively use FindMarkers() with test LR and latent variable of batch

categories <- unique(mm.sobj@meta.data$Major.cell.type)
for (cat in categories) {
  print(cat)
  try({
    cat.sobj <- subset(mm.sobj, subset = (Major.cell.type == cat))
    Idents(cat.sobj) <- cat.sobj@meta.data$Tumor.or.healthy
    markers <- FindConservedMarkers(cat.sobj, ident.1 = "t", ident.2 = "h", grouping.var = "Library.prep.batch")
    markers <- cbind(gene = rownames(markers), markers)
    write.csv(markers, paste0("./outs/DEG/Markers_Tumor-Healthy_MajorCellType_", cat, ".csv"), row.names = F)
  })
}

# Summarize DEG
# marker.files <- list.files("./outs/DEG/", pattern = "Markers_Tumor-Healthy_MajorCellType_*")
# markers.ls <- list()
# for (fn in marker.files) {
#   
# }








