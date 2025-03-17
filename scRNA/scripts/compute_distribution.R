# 3/17/2025
# Compute Cell Distributions
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0
library(openxlsx)
# library(speckle)
library(limma)
source("./packages/propeller.R")
source("./packages/getTransformedProps.R")
source("./packages/propeller.ttest.R")

mm.sobj <- readRDS("./outs/SeuratObject/SeuratObj_init.rds")

# Compute cell distribution

propeller.df <- propeller(clusters = mm.sobj@meta.data$Major.cell.type,
                          sample = mm.sobj@meta.data$Biological.replicate,
                          group = mm.sobj@meta.data$Tumor.or.healthy)

# Save Outputs

# CSV
write.csv(propeller.df, "./outs/CellDistributions/Propeller.csv", row.names = F)
# XLSX
# TODO(dawelsch): Fix number of digits for select columns
wb <- createWorkbook()
addWorksheet(wb, "Propeller")
writeDataTable(wb, "Propeller", propeller.df)
saveWorkbook(wb, "./outs/CellDistributions/Propeller.xlsx", overwrite = TRUE)











