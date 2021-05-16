#' Simple PCA analysis
#'
#' Performs a principal component analysis of a
#' normalized expression/abundance matrix and plot the samples
#' in the bidimensional space formed by the two first principal
#' components. Samples are colored according to the variable of
#' interest in the sample information data frame.
#'
#' @param mat The input matrix containing normalized data or counts
#' @param sampInfo The sample information data frame
#' @param groupCol The column name of the grouping variable
#'
#' @return A ggplot2 object containing the resulting plot
#' @export
#'
#' @importFrom ggplot2 autoplot theme_bw
#' @importFrom stats prcomp
#' @importFrom dplyr %>%
#' @import ggfortify
#'
pcaPlot <- function(mat, sampInfo, groupCol) {

  pcRes <- mat %>%
    t() %>%
    stats::prcomp()

  rownames(sampInfo) <- sampInfo[,1]

  p <- ggplot2::autoplot(object = pcRes, data = sampInfo, colour = groupCol) +
    ggplot2::theme_bw()

  return(p)

}

#' Simple Volcano Plot
#'
#' Creates a simple volcano plot summarizing the number of differentially altered
#' features by comparison.
#'
#' @param intTable Data frame with differential analysis results
#' @param logFcCutoff Fold change cutoff
#' @param pCutoff Adjusted P value cutoff
#'
#' @return A ggplot2 object containing the resulting plot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
volcanoPlot <- function(intTable, logFcCutoff = log2(2), pCutoff = 0.05) {

  intTable  <- intTable %>%
    mutate(status = case_when(logFc >= logFcCutoff & pAdj <= pCutoff ~ "Up",
                              logFc <= -logFcCutoff & pAdj <= pCutoff ~ "Down",
                              TRUE ~ "Other")) %>%
    mutate(comparison = factor(comparison, levels = unique(comparison)))

  allFeatTable <-  intTable %>%
    group_by(comparison) %>%
    summarise(total_n = n())

  annotDf <- intTable %>%
    group_by(comparison, status) %>%
    summarise(n = n()) %>%
    inner_join(x = ., y = allFeatTable, by = c("comparison")) %>%
    mutate(pct = (n / total_n) * 100) %>%
    mutate(xpos = case_when(status == "Up" ~ Inf,
                            status == "Down" ~ -Inf,
                            TRUE ~ 0),
           ypos = Inf,
           hjust = case_when(status == "Up" ~ 1,
                             status == "Down" ~ 0,
                             TRUE ~ 0.5),
           vjust = 1,
           label = paste0("N = ", n, "\n", round(pct, digits = 2), "%"))

  outP <- intTable %>%
    ggplot(aes(x = logFc, y = -log10(pAdj), label = feature, color = status)) +
    geom_vline(xintercept = c(logFcCutoff, -logFcCutoff), lty = 2) +
    geom_hline(yintercept = -log10(pCutoff), lty = 2) +
    geom_label(data = annotDf, aes(x=xpos,y=ypos,hjust=hjust,vjust=vjust,label=label)) +
    geom_point(alpha = 0.3) +
    xlab("Log2 Fold Change") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Other" = "black")) +
    facet_wrap(facets = vars(comparison)) +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black", fill ="white"),
          strip.placement = "inside",
          panel.spacing = unit(0.2, "lines"),
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(colour="black"))

  return(outP)

}

#' Heatmap plot
#'
#' Creates a heatmap representation of the more variable features in a matrix
#'
#' @param mat The matrix to be plotted
#' @param sampInfo The sample information data frame
#' @param groupCol The column name of the grouping variable
#' @param scaleBy character indicating if the values should be centered and scaled in either the
#' row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param nTop Number of features to be plotted after ranking them by the standard deviation
#' @param ... Arguments to be passed to \link[pheatmap]{pheatmap}
#'
#' @return The phetamap plot
#' @export
#'
#' @import dplyr
#' @import pheatmap
#'
heatmapPlot <- function(mat, sampInfo, groupCol, scaleBy = "row", nTop = 100, ...) {

  rownames(sampInfo) <- sampInfo[, 1]
  sampInfo <- dplyr::select(sampInfo, !!sym(groupCol))

  intFeatures <- apply(mat, 1, sd) %>%
    sort(. , decreasing = TRUE) %>%
    .[1:nTop] %>%
    names()
  mat <- mat[intFeatures,]

  outHm <- pheatmap(mat, scale = scaleBy, annotation_col = sampInfo, ...)

  return(outHm)

}

#' Violin plot
#'
#' Creates a per-column violin plot to explore data density before and after normalization
#'
#' @param mat The matrix to be plotted
#'
#' @return A ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#'
violinPlot <- function(mat) {

  p <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(-feature, names_to = "sample") %>%
    ggplot(aes(x = sample, y = value, fill = sample)) +
    geom_violin() +
    guides(x = guide_axis(angle = 60)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), legend.position = "none")
  return(p)

}

#' GSEA plot
#'
#' Plots the results of a gsea analysis
#'
#' @param gseaResDf The output data frame from \link[biokit]{gseaFromStats}
#' @param pCutoff Adjusted P value cutoff to filter features before plotting
#'
#' @return A ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
gseaPlot <- function(gseaResDf, pCutoff = 0.05) {

  outP <- gseaResDf %>%
    mutate(status = ifelse(NES > 0, "Up", "Down")) %>%
    subset(padj <= pCutoff) %>%
    ggplot(aes(x = comparison, y = pathway, fill = NES, size = -log10(padj), shape = status)) +
    geom_point() +
    scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
    scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0) +
    guides(x = guide_axis(angle = 60)) +
    theme_bw() +
    ggtitle(paste0("MSigDb Hallmarks Results (FDR <= ", pCutoff, ")")) +
    theme(text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5))

  return(outP)

}


