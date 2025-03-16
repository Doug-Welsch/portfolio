# 3/16/2025
# Generating Custom SPRING Plots
# Doug Welsch

setwd("~/RProjects/portfolio/scRNA/")

library(Seurat) # TODO(dawelsch): 4.3.0
library(ggplot2)

generateColors <- function(n, seed = NULL) {
  # Generates a character vector of n evenly spaced HCL hues. If seed is
  # provided, shuffles colors
  
  # Select colors
  hues <- seq(0, 360, length.out = n + 1)[-1]
  colors <- hcl(h = hues, l = 65, c = 100)
  # Shuffle colors
  if (!is.null(seed)) {
    set.seed(seed)
    colors <- sample(colors, n)
  }
  return(colors)
}

plotMetaSPRING <- function(sobj,
                       fill = NULL,
                       highlight = NULL,
                       size = 2,
                       alpha = 1,
                       title = "",
                       seed = 123,
                       save = NULL,
                       width = 8,
                       height = 8) {
  # Generates categorical SPRING reduction plot
  
  set.seed(seed)
  spring.df <- as.data.frame(sobj@reductions$spring@cell.embeddings)
  
  spring.df$fill <- "Unassigned"
  cols.fill <- list("Unassigned" = "grey")
  if (!is.null(fill))  {
    spring.df$fill <- mm.sobj@meta.data[[fill]]
    # Assign colors
    fill.cats <- unique(spring.df$fill)
    tmp.cols <- generateColors(n = length(fill.cats))
    names(tmp.cols) <- fill.cats
    for (i in seq_along(tmp.cols)) {
      cols.fill[[names(tmp.cols)[i]]] <- tmp.cols[i]
    }
  }
  if (!is.null(highlight)) {
    for (cat in names(cols.fill)) {
      if (cat != highlight) {
        cols.fill[[cat]] <- "grey"
      }
    }
  }
  
  # Plot
  spring.df <- spring.df[sample(1:nrow(spring.df), nrow(spring.df), replace = F), ]
  p <- ggplot(spring.df, aes(x=SPRING_1, y=SPRING_2)) +
    geom_point(color="black", aes(fill=fill), pch=21, size=size, alpha=alpha) +
    ggtitle(title) + 
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          aspect.ratio = 1) +
    scale_fill_manual(values = cols.fill) +
    labs(fill=if(is.null(fill)) "" else gsub(".", " ", fill, fixed = T), color="") + 
    guides(fill = guide_legend(override.aes = list(size=5)))
  
  # Save
  if (!is.null(save)) {
    ggsave(filename = save,
           width = width, height = height,
           plot = p)
  }
  
  return (p)
}

mm.sobj <- readRDS("./outs/SeuratObject/SeuratObj_init.rds")

# Major Cell Types
plotMetaSPRING(mm.sobj,
               fill = "Major.cell.type",
               save = "./outs/SPRING/SPRING_MajorCellTypes.pdf")
cell.types <- unique(mm.sobj@meta.data$Major.cell.type)
for (ct in cell.types) {
  print(ct)
  plotMetaSPRING(mm.sobj,
                 fill = "Major.cell.type",
                 highlight = ct,
                 save = paste0("./outs/SPRING/SPRING_MajorCellTypes_", gsub(" ", "", ct), ".pdf"))
}


