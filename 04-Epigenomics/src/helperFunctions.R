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
# 
# do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
#                        do_labels = TRUE, nice_names, 
#                        palette_use = color.list(),
#                        pt_size = 4, point_size = .5, base_size = 12, 
#                        do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
#   umap_use <- umap_use[, 1:2]
#   colnames(umap_use) <- c('X1', 'X2')
#   plt_df <- umap_use %>% data.frame() %>% 
#     cbind(meta_data) %>% 
#     dplyr::sample_frac(1L) 
#   plt_df$given_name <- plt_df[[label_name]]
#   
#   if (!missing(nice_names)) {
#     plt_df %<>%
#       dplyr::inner_join(nice_names, by = "given_name") %>% 
#       subset(nice_name != "" & !is.na(nice_name))
#     
#     plt_df[[label_name]] <- plt_df$nice_name        
#   }
#   
#   plt <- plt_df %>% 
#     ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
#     theme_test(base_size = base_size) + 
#     theme(panel.background = element_rect(fill = NA, color = "black")) + 
#     guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
#                                                     shape = 16, size = 4)), 
#            alpha = FALSE) +
#     scale_color_manual(values = palette_use) + 
#     scale_fill_manual(values = palette_use) +    
#     theme(plot.title = element_text(hjust = .5)) + 
#     labs(x = "PC 1", y = "PC 2") 
#   
#   if (do_points) 
#     plt <- plt + geom_point(shape = '.')
#   if (do_density) 
#     plt <- plt + geom_density_2d()    
#   
#   
#   if (no_guides)
#     plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
#   
#   if (do_labels) {
#     data_labels <- plt_df %>% 
#       dplyr::group_by_(label_name) %>% 
#       dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>% 
#       dplyr::ungroup()
#     
#     plt <- plt + geom_label(data = data_labels, label.size = NA,
#                             aes_string(label = label_name), 
#                             color = "white", size = pt_size, alpha = 1,
#                             segment.size = 0) +
#       guides(col = FALSE, fill = FALSE)
#   }
#   
#   return(plt)
# }
# 
# ## get the explained variance from the total
# getVarianceExplainedSeurat <- function(
#   object, assay = "RNA", slot = "scale.data", reduction = "pca"
# ) {
#   mat <- Seurat::GetAssayData(object, assay = assay, slot = slot)
#   pca <- object[[reduction]]
#   # Get the total variance:
#   total_variance <- sum(matrixStats::rowVars(mat))
#   eigValues <- (pca@stdev)^2  ## EigenValues
#   return(eigValues / total_variance * 100)
# }
# 
# 
# plotDensities2 <- function(
#   matrix, 
#   color.by = NULL,
#   title = "", 
#   xlab = "",
#   ylim = 0.27,
#   cols = NULL, 
#   cutoff = NULL
# ) {
#   nsamples <- ncol(matrix)
#   plot(density(matrix[, 1]), col = cols[1], 
#        lwd = 2, las = 1, ylim = c(0, ylim), main = "", xlab = "")
#   grid()
#   title(main = title, xlab = xlab)
#   if (!is.null(cutoff)) abline(v = cutoff, lty = 3)
#   for (i in 2:nsamples){
#     den <- density(matrix[, i])
#     lines(den$x, den$y, col = cols[i], lwd = 2)
#   }
# }
# 
# plotDensities3 <- function(
#   matrix, 
#   color.by = NULL,
#   title = "", 
#   xlab = "",
#   ylim = 0.27,
#   cols = NULL, 
#   cutoff = NULL
# ) {
#   densList <- apply(
#     X = matrix, 
#     MARGIN = 2, 
#     FUN = function(x) {
#       dens <- density(x)
#       return(data.frame(x = dens$x, y = dens$y))
#     }
#   )
#   densDF <- do.call(what = rbind, densList)
#   densDF$Sample <- gsub(
#     pattern = "\\.\\d*", replacement = "", x = rownames(densDF)
#   )
#   if (!is.null(group.by)) densDF$Group <- group.by
#   nsamples <- ncol(matrix)
#   p <- ggplot(densDF, aes(x = x, y = y)) + 
#     p + geom_smooth(method = "loess")
#   
# }
# 
# 
## plot PCA with ggplot2
plotPCA2 <- function(pcaObject, col.points, shape.points = NULL, palette,
                     legend.col, point.size = 3, title = "", pcs = c(1, 2)){
  variance <- round(factoextra::get_eigenvalue(pcaObject)[pcs, 2], 1)

  p <- ggplot(data.frame(pcaObject[["x"]]),
              aes(x = .data[[paste0("PC", pcs[1])]], y = .data[[paste0("PC", pcs[2])]], color = col.points, shape = shape.points)) +
    geom_point(size = point.size) +
    scale_color_manual(name = legend.col, values = palette) +
    xlab(paste0("PC", pcs[1], " (", variance[1], "%)")) + ylab(paste0("PC", pcs[2], " (", variance[2], "%)")) +
    geom_vline(xintercept = 0, linetype="dashed") + geom_hline(yintercept = 0, linetype="dashed")  +
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    ggtitle(title) + theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  return(p)
}


plotPvalDist <- function(top.table, title, pval = "adj.P.Val") {
  dif <- sum(top.table[[pval]] <= 0.05)
  ggplot(data = as.data.frame(top.table), aes(x = .data[[pval]])) + 
    geom_histogram(bins = 100) + ggtitle(paste0(title, " (", dif, " DEGs)")) + 
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size=1) +
    theme_classic() + theme(plot.title = element_text(face = "bold"))  
}

MAplot.new <- function(
    res, 
    comp, 
    th.pval = 0.05, 
    th.fold = 1.5,
    col.points = "#45824a",
    col.sig = "#e68337",
    cex = 0.2
) {
  res.sig <- res %>% dplyr::filter(padj <= th.pval, abs(log2FoldChange) >= th.fold) 
  res.not <- res[!rownames(res) %in% rownames(res.sig), ]
  idx <- nrow(res.sig)
  
  xmin <- floor(min(log2(res$baseMean)))
  xmax <- ceiling(max(log2(res$baseMean)))
  ymin <- floor(min(res$log2FoldChange))
  ymax <- ceiling(max(res$log2FoldChange))
  smoothScatter(
    log2(res.not$baseMean),
    res.not$log2FoldChange,
    pch = 20,
    cex = cex, 
    col = col.points,
    colramp = colorRampPalette(c("white", col.points)),                      
    xaxp = c(xmin, xmax, xmax - xmin),
    xlim = c(xmin, xmax),
    xlab = 'log concentration',
    yaxp = c(ymin, ymax, (ymax - ymin)), 
    ylim = c(ymin, ymax),
    ylab = sprintf('log fold change: ', comp),
    main = sprintf(
      'Binding Affinity: %s (%s peaks, padj <= %1.2f, abs(logFC) >= %1.2f)',
      comp, sum(idx), th.pval, th.fold
    )
  )
  points(
    log2(res.sig$baseMean), res.sig$log2FoldChange, pch = 20, cex = cex, col = col.sig
  )
  abline(h = 0, col = 'dodgerblue')
}



#  
# ## customize volcano plot
# volcanoplotCustom <- function(
#   tableRes, 
#   xLim, 
#   yLim, 
#   column.pvalue = "adj.P.Val",
#   style = "pval",
#   colors = c("#9e2121", "#2d7fcc", "#4a5054"),
#   legend = TRUE, 
#   pos.legend = "topleft", 
#   legend.text = NULL, 
#   line.color = "black",
#   highlight = TRUE, 
#   numLables = 5,
#   main = "Volcano plot"
# ){
#   if (style == "pval") {
#     tableRes$logPval <- -log10(
#       p.adjust(tableRes[[column.pvalue]], method = "BH", n = length(tableRes[[column.pvalue]]))
#     )  
#     ylab <- "-log10(p-value)"
#   } else if (style == "log(CPM)") {
#     tableRes$logPval <- tableRes$logCPM
#     ylab <- "log(CPM)"
#   }
#   
#   # tableRes$logPval <- tableRes$logCPM
#   upGenes <- subset(tableRes, logFC > xLim & logPval > yLim)
#   downGenes <- subset(tableRes, logFC < -xLim & logPval > yLim)
#   nonGenes <- tableRes[!(rownames(tableRes) %in% rownames(upGenes)) &
#                          !(rownames(tableRes) %in% rownames(downGenes)), ]
#   plot(x = tableRes$logFC, y = tableRes$logPval, type = "n",
#        xlab = "logFC", ylab = ylab, main = main)
#   points(x = upGenes$logFC, y = upGenes$logPval, col = colors[1], pch = 20, cex = 0.4)
#   points(x = downGenes$logFC, y = downGenes$logPval, col = colors[2], pch = 20, cex = 0.4)
#   points(x = nonGenes$logFC, y = nonGenes$logPval, col = colors[3], pch = 20, cex = 0.4)
#   abline(v = c(xLim, -xLim), h = yLim, col = line.color, lty = "longdash")
#   
#   if (legend) {
#     if (is.null(legend.text)){
#       legend.text = c(paste0("Up-regulated: ", nrow(upGenes)), 
#                       paste0("Down-regulated: ", nrow(downGenes)),
#                       paste0("Non-significant: ", nrow(nonGenes)))
#     }
#     legend(pos.legend, legend = legend.text, fill = colors, 
#            cex = 0.8, bg = "white")  
#   }
#   
#   if (highlight) {
#     labelsUp <- head(with(upGenes, upGenes[order(-logFC, -logPval),], n = numLabels))
#     geneNamesUp <- match(rownames(labelsUp), rownames(tableRes))
#     maptools::pointLabel(labelsUp$logFC, labelsUp$logPval, 
#                          rownames(tableRes)[geneNamesUp], cex = 0.7)
#     
#     labelsDown <- head(with(downGenes, downGenes[order(-logPval, -logFC),], n = numLabels))
#     geneNamesDown <- match(rownames(labelsDown), rownames(tableRes))
#     maptools::pointLabel(labelsDown$logFC, labelsDown$logPval, 
#                          rownames(tableRes)[geneNamesDown], cex = 0.7)
#   }
# }
# 
# ################################################################################
# ## new functions
# 
# ## volcano plot with ggplot2
# interactiveVolcanoPlot <- function(
#   table, x.lim, y.lim, size = 0.5, 
#   title = "Volcano plot", interactive = TRUE, annotation = NULL
# ) {
#   table$Gene <- rownames(table)
#   table$log10.adj.pval <- -log10(table[["adj.P.Val"]])
#   table$Significant <- "Non-significant"
#   table[table$logFC >= x.lim & table$log10.adj.pval >= y.lim, "Significant"] <- "Up-regulated"
#   table[table$logFC <= -x.lim & table$log10.adj.pval >= y.lim, "Significant"] <- "Down-regulated"
#   table$Significant <- factor(table$Significant, levels = c("Up-regulated", "Down-regulated", "Non-significant"))
#   
#   if (!is.null(annotation)) {
#     table.sel <- table[annotation, ]
#     # ann <- geom_label(
#     #   data=table.sel, 
#     #   aes(x = logFC, y = log10.adj.pval, color = Significant, label = Gene), 
#     #   color='black', fill = "white"
#     # ) 
#     ann <- ggrepel::geom_label_repel(
#       data = table.sel, 
#       aes(x = logFC, y = log10.adj.pval, color = Significant, label = Gene), 
#       color='black', fill = "white",
#       box.padding   = 0.35, 
#       point.padding = 0.5,
#       segment.color = 'black'
#     )
#     point <- geom_point(
#       data=table.sel, 
#       aes(x = logFC, y = log10.adj.pval, label = Gene), color = "#32a852", size = 1
#     )
#   } else {
#     ann <- NULL
#     point <- NULL
#   }
#   
#   p <- ggplot(table, aes(x = logFC, y = log10.adj.pval, color = Significant, label = Gene)) + 
#     geom_point(size = size) + ggtitle(title) + 
#     scale_color_manual(values = c("#9e2121", "#2d7fcc", "#4a5054")) + 
#     theme_bw() + 
#     theme(
#       plot.title = element_text(face = "bold"), legend.title = element_text(face = "bold")
#     ) + ann + point
#   if (interactive) ggplotly(p) 
#   else p
# }
# 
# ## boxplot with ggplot
# findGene <- function(matrix.counts, pattern.gene, fix = FALSE) {
#   grep(pattern = pattern.gene, x = rownames(matrix.counts), ignore.case = TRUE, fixed = fix)
# }
# boxplotGene <- function(
#     matrix.counts, 
#     samples.metadata, 
#     x.var,
#     color.var,
#     disc.var = "CellType",
#     disc.fact = "cDC1",
#     gene, 
#     fix = FALSE,
#     z.score = FALSE
# ) {
#   pos.gene <- findGene(matrix.counts, gene, fix = fix)
#   gene.expression <- matrix.counts[pos.gene, , drop = FALSE]
#   lapply(
#     X = seq(nrow(gene.expression)), FUN = function(x) {
#       dfPlot <- cbind(Gene = gene.expression[x, ], samples.metadata)
#       if (color.var == x.var) {
#         ## plot
#         # dfPlot <- dfPlot %>% filter(dfPlot[[disc.var]] == disc.fact)
#         if (z.score) dfPlot$Gene <- scale(dfPlot$Gene)[, 1]
#         ggplot(dfPlot, aes(x = .data[[x.var]], y = Gene, fill = .data[[color.var]])) + 
#           # ggdist::stat_halfeye(
#           #   adjust = .5, 
#           #   width = .6, 
#           #   .width = 0, 
#           #   justification = -.3, 
#           #   point_colour = NA) + 
#           geom_boxplot(
#             width = 0.5, 
#             outlier.shape = NA,
#             position=position_dodge(0.60)
#           ) + 
#           geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), stackratio=0.8, dotsize=0.5) + 
#           coord_cartesian(xlim = c(1.2, NA), clip = "off") + 
#           scale_fill_manual(values = colorsConditions[[color.var]]) + 
#           theme_bw() + ggtitle(paste0("Expression levels of ", rownames(gene.expression)[x])) + 
#           theme(plot.title = element_text(face = "bold")) + ylab("Expression level") 
#       } else {
#         if (z.score) dfPlot$Gene <- scale(dfPlot$Gene)[, 1]
#         ## plot
#         ggplot(dfPlot, aes(x = .data[[x.var]], y = Gene, fill = .data[[color.var]])) + 
#           # ggdist::stat_halfeye(
#           #   adjust = .5, 
#           #   width = .6, 
#           #   .width = 0, 
#           #   justification = -.3, 
#           #   point_colour = NA) + 
#           geom_boxplot(
#             width = 0.5, 
#             outlier.shape = NA,
#             position=position_dodge(0.60)
#           ) + 
#           geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),stackratio=0.8, dotsize=0.5) + 
#           coord_cartesian(xlim = c(1.2, NA), clip = "off") + 
#           scale_fill_manual(values = colorsConditions[[color.var]]) + 
#           theme_bw() + ggtitle(paste("Expression levels of ", rownames(gene.expression)[x])) + 
#           theme(plot.title = element_text(face = "bold")) + ylab("Expression level") 
#       }
#     }
#   )
# }
# 
# heatmapGenes <- function(
#   expr.data, genes, 
#   samples.metadata, color,
#   colnames = TRUE, fix.order = TRUE,
#   ...
# ) {
#   matrix.expr <- expr.data[genes, ]
#   matrix.expr.z <- t(scale(t(matrix.expr)))
#   colnames(matrix.expr.z) <- samples.metadata$title
#   if (fix.order) {
#     order.cols <- 1:12
#     cluster.cols <- FALSE
#   } else{
#     order.cols <- NULL
#     cluster.cols <- TRUE
#   } 
#   # ha <- HeatmapAnnotation(
#   #   bar = names(color),
#   #   col = list(bar = color)
#   # )
#   Heatmap(
#     matrix = matrix.expr.z,
#     row_names_gp = gpar(fontsize = 8),
#     column_dend_reorder = FALSE,
#     column_names_gp = gpar(fontsize = 8), 
#     name = "Z-score",
#     column_title_gp = gpar(fontface = "bold"), 
#     column_split = samples.metadata[["Condition"]],
#     column_title = "Heatmap", 
#     cluster_columns = cluster.cols, 
#     column_order = order.cols,
#     show_row_dend = FALSE, 
#     show_row_names = TRUE,
#     show_column_names = colnames,
#     rect_gp = gpar(col = "white", lwd = 1)
#     # top_annotation = ha
#   )
# }
# 
# ## heatmap for deg genes (use it carefully, it's not completely checked)
# heatmapDEG <- function(
#   matrix.expr,
#   table.deg,
#   samples.metadata,
#   logFC.cut = 2.5,
#   pval.cut = 0.05,
#   n.genes = 500,
#   n.genes.ann = 10,
#   title = "DEG",
#   split.var = "CellType",
#   logFC.th = 8,
#   rigth.ann = TRUE,
#   ...
# ) {
#   ## deg
#   deg <- table.deg %>% filter(abs(logFC) > logFC.cut, adj.P.Val <= pval.cut)
#   # deg <- deg[order(abs(deg$logFC), deg$adj.P.Val), ]
#   if (nrow(deg) > n.genes) {
#     up <- ceiling(n.genes / 2)
#     down <- n.genes - up
#     deg <- rbind(
#       arrange(filter(deg, logFC >= 0), adj.P.Val)[seq(up), ],
#       arrange(filter(deg, logFC <= 0), adj.P.Val)[seq(down), ]
#     )
#     # deg <- deg[seq(n.genes), ]
#   }
#   matrix.expr <- matrix.expr[rownames(deg), ]
#   matrix.expr.z <- t(scale(t(matrix.expr)))
#   rownames(matrix.expr.z) <- deg[["SYMBOL"]]
#   topGenes <- order(deg$logFC, decreasing = TRUE)[seq(n.genes.ann)]
#   downGenes <- order(deg$logFC, decreasing = FALSE)[seq(n.genes.ann)]
#   if (isTRUE(rigth.ann)) {
#     hh <- rowAnnotation(
#       link = anno_mark(
#         at = c(topGenes, downGenes),
#         labels = deg[c(topGenes, downGenes), "SYMBOL"],
#         labels_gp = gpar(fontsize = 10),
#         padding = unit(1, "mm")))
#   } else {
#     hh <- NULL
#   }
#   Heatmap(
#     matrix = matrix.expr.z,
#     row_names_gp = gpar(fontsize = 8),
#     column_names_gp = gpar(fontsize = 8), name = "Z-score",
#     column_title_gp = gpar(fontface = "bold"),
#     column_split = samples.metadata[[split.var]],
#     column_title = title,
#     right_annotation = hh,
#     cluster_columns = FALSE,
#     ...
#   )
# }
# 
# 
# ## heatmap for deg genes (use it carefully, it's not completely checked)
# heatmapDEG2 <- function(
#   matrix.expr,
#   table.deg,
#   samples.metadata,
#   logFC.cut = 2.5,
#   pval.cut = 0.05,
#   n.genes.ann = 10,
#   title = "DEG",
#   split.var = "CellType",
#   logFC.th = 8,
#   rigth.ann = TRUE,
#   ...
# ) {
#   ## deg
#   deg <- table.deg %>% filter(abs(logFC) > logFC.cut, adj.P.Val <= pval.cut)
#   # deg <- deg[order(abs(deg$logFC), deg$adj.P.Val), ]
#   # if (nrow(deg) > n.genes) {
#   #   up <- ceiling(n.genes / 2)
#   #   down <- n.genes - up
#   #   deg <- rbind(
#   #     arrange(filter(deg, logFC >= 0), adj.P.Val)[seq(up), ],
#   #     arrange(filter(deg, logFC <= 0), adj.P.Val)[seq(down), ]
#   #   )
#   #   # deg <- deg[seq(n.genes), ]
#   # }
#   matrix.expr <- matrix.expr[rownames(deg), ]
#   matrix.expr.z <- t(scale(t(matrix.expr)))
#   rownames(matrix.expr.z) <- deg[["SYMBOL"]]
#   topGenes <- order(deg$logFC, decreasing = TRUE)[seq(n.genes.ann)]
#   downGenes <- order(deg$logFC, decreasing = FALSE)[seq(n.genes.ann)]
#   if (isTRUE(rigth.ann)) {
#     hh <- rowAnnotation(
#       link = anno_mark(
#         at = c(topGenes, downGenes),
#         labels = deg[c(topGenes, downGenes), "SYMBOL"],
#         labels_gp = gpar(fontsize = 10),
#         padding = unit(1, "mm")))
#   } else {
#     hh <- NULL
#   }
#   Heatmap(
#     matrix = matrix.expr.z,
#     row_names_gp = gpar(fontsize = 8),
#     column_names_gp = gpar(fontsize = 8),
#     name = "Z-score",
#     column_title_gp = gpar(fontface = "bold"),
#     column_split = samples.metadata[[split.var]],
#     column_title = title,
#     right_annotation = hh,
#     ...
#   )
# }
# 
# ## fgsea atuomated. see code just in case
# fgseaAutomated <- function(table, pathways, criterion = "logFC") {
#   sortedTable <- table
#   if (criterion == "logFC") {
#     sortedTable$criterion <- with(sortedTable, logFC)  
#   } else if (criterion == "logFC.adj.pval") {
#     sortedTable$criterion <- with(sortedTable, logFC * -log10(adj.P.Val))  
#   }
#   sortedTable <- sortedTable[order(sortedTable$criterion, decreasing = TRUE), ]
#   genesSorted <- sortedTable$criterion
#   names(genesSorted) <- rownames(sortedTable)
#   ## fgsea
#   fgseaRes <- fgseaMultilevel(pathways, stats = genesSorted, maxSize = 500)
#   topUp <- fgseaRes[NES > 0][head(order(padj), n = 10), pathway]
#   topDown <- fgseaRes[NES < 0][head(order(padj), n = 10), pathway]
#   topAll <- c(topUp, rev(topDown))
#   grid.newpage()
#   grDevices::dev.interactive()
#   plotGseaTable(pathways[topAll], genesSorted, fgseaRes, gseaParam = 0.5)  
#   return(fgseaRes)
# }
# 
# ## overrepresentation analysis with kegg using clusterProfiler
# keggORA <- function(
#   genes, changeGenes = TRUE, fromType = "SYMBOL", title = "KEGG results"
# ) {
#   if (changeGenes) {
#     genesEntrez <- bitr(
#       genes, fromType = "SYMBOL",
#       toType = c("ENTREZID"),
#       OrgDb = org.Mm.eg.db
#     )
#     genesEntrez <- genesEntrez$ENTREZID
#   } else {
#     genesEntrez <- genes
#   }
#   universeEntrez <- bitr(
#     genesMetadata$ENSEMBL, fromType = "ENSEMBL",
#     toType = "ENTREZID",
#     OrgDb = org.Mm.eg.db
#   )
#   
#   res <- clusterProfiler::enrichKEGG(
#     gene = genesEntrez, 
#     universe = universeEntrez,
#     organism = "mmu"
#   )
#   ## plot
#   print(
#     barplot(
#       res,
#       drop = TRUE,
#       showCategory = 15,
#       title = title,
#       font.size = 8
#     ) + theme(plot.title = element_text(face = "bold"))
#   )
#   return(res)
# }
# 
# ## overrepresentation analysis with GO using clusterProfiler
# goORA <- function(genes, changeGenes = TRUE, fromType = "SYMBOL", title = "GO results") {
#   if (changeGenes) {
#     genesEnsembl <- bitr(
#       genes, fromType = "SYMBOL",
#       toType = c("ENSEMBL"),
#       OrgDb = org.Mm.eg.db
#     )
#     genesEnsembl <- genesEnsembl$ENSEMBL
#   } else {
#     genesEnsembl <- genes
#   }
#   universe <- genesMetadata$ENSEMBL
#   res <- clusterProfiler::enrichGO(
#     gene = genesEnsembl,
#     universe = universe,
#     OrgDb = "org.Mm.eg.db", 
#     keyType = 'ENSEMBL',
#     readable = T,
#     ont = "ALL",
#     pvalueCutoff = 0.05, 
#     qvalueCutoff = 0.10
#   )
#   print(
#     barplot(
#       res, 
#       drop = TRUE, 
#       showCategory = 15, 
#       title = title,
#       font.size = 8
#     ) + theme(plot.title = element_text(face = "bold"))
#   )
#   return(res)
# }
# 
# 
