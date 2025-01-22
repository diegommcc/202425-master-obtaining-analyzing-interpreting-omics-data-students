###############################################################################
################## Helper functions for single-cell analysis ##################
###############################################################################
# Author: Diego Mañanes Cayero
# Date: 21/02/24
# Description: set of functions used in single-cell RNA-seq analysis.
###############################################################################

## list of colors used for discrete data
color.list <- function() {
  color.list.2 <- c(
    RColorBrewer::brewer.pal(12, "Paired"), "#d45b91", "#374738",
    RColorBrewer::brewer.pal(8, "Pastel2"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    "#333333", "#5D5D5D",
    "#888888", "#B3B3B3"
  )
  color.list.2[11] <- "#e3dc5b"
  color.list.2[15] <- "#60c4b4"
  color.list.2[16] <- "#7e8046"
  # color.list.2[17] <- "#3526bd"
  return(color.list.2)
}

color.numeric <- colorRampPalette(
  colors = c("#cfcfcf", "#958da6", "#7452b3", "#070084")
)(100)

## plot PCA
plotPCA <- function(
    pcaObject, col.points, shape.points = NULL, palette,
    legend.col, point.size = 3, title = "", pcs = c(1, 2)
){
  variance <- round(factoextra::get_eigenvalue(pcaObject)[pcs, 2], 1)
  
  p <- ggplot(
    data.frame(pcaObject[["x"]]),
    aes(
      x = .data[[paste0("PC", pcs[1])]], 
      y = .data[[paste0("PC", pcs[2])]], 
      color = col.points, shape = shape.points
    )
  ) + geom_point(size = point.size) +
    scale_color_manual(name = legend.col, values = palette) +
    xlab(paste0("PC", pcs[1], " (", variance[1], "%)")) + 
    ylab(paste0("PC", pcs[2], " (", variance[2], "%)")) +
    geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = 0, linetype = "dashed")  +
    guides(fill = guide_legend(override.aes = list(shape=21))) +
    ggtitle(title) + theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}

