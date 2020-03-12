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
#' @return A ggplot containing the default volcano plots
#' @export
#'
#' @import ggplot2
#'
defaultVolcano <- function(x, fcCol = "logFc", pCol = "pAdj", statusCol = "status", compCol = NULL, featCol = "feature", upColor = "red",
    noChangeColor = "black", downColor = "blue", upLabel = "Up", noChangeLabel = "No change", downLabel = "Down") {

    # prepare color values
    colVal <- c(upColor, noChangeColor, downColor)
    names(colVal) <- c(upLabel, noChangeLabel, downLabel)
    # tidy eval for column names (https://tidyeval.tidyverse.org/introduction.html)
    p <- ggplot(x, aes(x = !!sym(fcCol), y = -log10(!!sym(pCol)), color = !!sym(statusCol), feature = !!sym(featCol))) + geom_point(alpha = 0.5) +
        scale_color_manual(values = colVal) + theme_bw()
    # if comparison column specified, facet wrap using it
    if (!is.null(compCol))
        p <- p + facet_wrap(facets = vars(!!sym(compCol)))
    return(p)

}
