#' Default volcano plot
#'
#' Creates a default volcano plot from comparison data frame.
#'
#' @param x Data frame with data to plot.
#' @param fcCol Fold change column.
#' @param pCol P value column
#' @param statusCol Status column.
#' @param compCol Comparison column.
#' @param featCol Feature column.
#' @param upColor Up reuglated features color.
#' @param noChangeColor Not chaning features color.
#' @param downColor Down regulated features color.
#' @param upLabel Up regulated features label.
#' @param noChangeLabel Not chaning features label.
#' @param downLabel Down regulated features label.
#'
#' @return A ggplot containing the default volcano plots.
#' @export
#'
#' @import ggplot2
#'
defaultVolcano <- function(x, fcCol = "logFc", pCol = "pAdj", statusCol = "status", compCol = NULL, featCol = "feature", upColor = "red", noChangeColor = "black", downColor = "blue", 
    upLabel = "Up", noChangeLabel = "No change", downLabel = "Down") {
    
    # prepare color values
    colVal <- c(upColor, noChangeColor, downColor)
    names(colVal) <- c(upLabel, noChangeLabel, downLabel)
    # tidy eval for column names (https://tidyeval.tidyverse.org/introduction.html)
    p <- ggplot(x, aes(x = !!sym(fcCol), y = -log10(!!sym(pCol)), color = !!sym(statusCol), feature = !!sym(featCol))) + geom_point(alpha = 0.5) + scale_color_manual(values = colVal) + 
        theme_bw()
    # if comparison column specified, facet wrap using it
    if (!is.null(compCol)) 
        p <- p + facet_wrap(facets = vars(!!sym(compCol)))
    return(p)
    
}


#' Default Principal Component Analysis Plot
#'
#' Creates a simple representation of the \link[stats]{prcomp} results formatted with the
#' \link[biokit]{pcaToList} function.
#'
#' @param x Result obtained with \link[biokit]{pcaToList}.
#' @param sampInfoDf Data frame containing sample information. Rownames should match \link[stats]{prcomp} matrix rownames.
#' @param groupCol Column of the data frame that indicates group of samples.
#' @param showLabel Show row labels?
#' @param encircle Encircle points by group?
#'
#' @return A ggplot containing the default pca plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggalt geom_encircle
#' @importFrom tibble rownames_to_column
#'
defaultPcaPlot <- function(x, sampInfoDf, groupCol, showLabel = TRUE, encircle = TRUE) {
    
    # check sampInfo rownames
    if (!all(rownames(sampInfoDf) == rownames(x$result$x))) 
        stop("Rownames of the sample information data frame could not be mapped to pca result.")
    # create df to plot
    toPlotDf <- cbind(tibble::rownames_to_column(sampInfoDf, var = "id"), x$result$x)
    # create plot
    p <- ggplot(toPlotDf, mapping = aes(x = !!sym("PC1"), y = !!sym("PC2"), color = !!sym(groupCol), fill = !!sym(groupCol), label = !!sym("id"))) + geom_point() + xlab(x$pcts[1]) + 
        ylab(x$pcts[2]) + theme_bw()
    # add label and encircle if required
    if (encircle) 
        p <- p + ggalt::geom_encircle(expand = 0, alpha = 0.3)
    if (showLabel) 
        p <- p + geom_label(fill = "white")
    return(p)
    
}
