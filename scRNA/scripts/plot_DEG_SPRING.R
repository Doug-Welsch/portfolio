# 3/18/2025
# Plotting Healthy Major Cell Type Positive DEG
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0
library(ggplot2)

mm.sobj <- readRDS("./outs/SeuratObject/SeuratObj_init.rds")
markers.summary <- read.csv("./outs/DEG/Summary_Healthy_MajorCellType.csv")

# Individual Gene Plots

for (i in 1:nrow(markers.summary)) {
  print(i)
  gene <- markers.summary$gene[i]
  cat <- markers.summary$Major.cell.type[i]
  
  # Plot
  p <- FeaturePlot(mm.sobj, reduction = "spring", features = gene) +
    theme(aspect.ratio = 1,
          plot.title = element_text(size = 24))
  # Save
  ggsave(filename = paste0("./outs/DEG/SPRING_", cat, "_", gene, ".pdf"),
         width = 6, height = 6,
         plot = p)
}

# Module Score
# TODO(dawelsch): Modify to be variable lists of significant genes

## Compute

mod.ls <- list()
for (cat in unique(markers.summary$Major.cell.type)) {
  mod.ls[[cat]] <- markers.summary[markers.summary$Major.cell.type == cat, ]$gene
}
mm.sobj <- AddModuleScore(
  object = mm.sobj,
  features = mod.ls,
  ctrl = 5,
  name = "ModuleScore",
  seed = 123
)
colnames(mm.sobj@meta.data)[grepl("ModuleScore", colnames(mm.sobj@meta.data), fixed=T)] <- paste0(names(mod.ls), ".ModuleScore")

## Visualize

for (cat in unique(markers.summary$Major.cell.type)) {
  print(cat)
  mod.name <- paste0(cat, ".ModuleScore")
  
  # Plot
  p <- FeaturePlot(mm.sobj, reduction = "spring", features = mod.name) +
    ggtitle(gsub(".", " ", mod.name, fixed = T)) +
    theme(aspect.ratio = 1,
          plot.title = element_text(size = 24))
  # Save
  ggsave(filename = paste0("./outs/DEG/SPRING_ModuleScore_", cat, ".pdf"),
         width = 6, height = 6,
         plot = p)
}

















