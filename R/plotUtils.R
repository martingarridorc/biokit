#' Default volcano plot
#'
#' Creates a default volcano plot from comparison data frame.
#'
#' @param x Data frame with data to plot.
#' @param fcCol Fold change column.
#' @param pCol P value column.
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
defaultVolcano <- function(x, fcCol = "logFc", pCol = "pAdj", statusCol = "status", compCol = NULL, featCol = "feature",
    upColor = "red", noChangeColor = "black", downColor = "blue", upLabel = "Up", noChangeLabel = "No change", downLabel = "Down") {

    # prepare color values
    colVal <- c(upColor, noChangeColor, downColor)
    names(colVal) <- c(upLabel, noChangeLabel, downLabel)
    # tidy eval for column names (https://tidyeval.tidyverse.org/introduction.html)
    p <- ggplot(x, aes(x = !!sym(fcCol), y = -log10(!!sym(pCol)), color = !!sym(statusCol), feature = !!sym(featCol))) +
        geom_point(alpha = 0.5) + scale_color_manual(values = colVal) + ggtitle("Volcano Plot") + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
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
#'
#' @return A ggplot containing the default pca plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#'
defaultPcaPlot <- function(x, sampInfoDf, groupCol, showLabel = TRUE) {

    # check sampInfo rownames
    if (!all(rownames(sampInfoDf) == rownames(x$result$x)))
        stop("Rownames of the sample information data frame could not be mapped to pca result.")
    # create df to plot
    toPlotDf <- cbind(tibble::rownames_to_column(sampInfoDf, var = "id"), x$result$x)
    # create plot
    p <- ggplot(toPlotDf, mapping = aes(x = !!sym("PC1"), y = !!sym("PC2"), color = !!sym(groupCol), fill = !!sym(groupCol),
        label = !!sym("id"))) + geom_point() + xlab(x$pcts[1]) + ylab(x$pcts[2]) + ggtitle("Principal Component Analysis") + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    # add label if required
    if (showLabel)
        p <- p + geom_label(fill = "white")
    return(p)

}

#' Default Over-representation analysis plot
#'
#' Creates a simple representation of the results of \link[clusterProfiler]{enricher} formatted
#' with \link[biokit]{cpResultsToDf}. The X axis represents column of interest while Y axis
#' represents the functional category. The presence of a point indicates the enrichment of a
#' functional category in a group of interest. The data frame can be filtered beforehand or
#' setting the pCutoff parameter.
#'
#' @param x Data frame with results to plot.
#' @param idColumn Column used as id of the analysed group, will be plotted in the X axis.
#' @param splitStatus Split id Column into up and down features?
#' @param facetByStatus Facet final plot by status of the analysed features?
#' @param upLabel Label used to indicate status of up-regulated features in the idColumn.
#' @param downLabel Label used to indicate status of up-regulated features in the idColumn..
#' @param idSep Separator used to paste group and status label.
#' @param upColor Color used to plot enrichments in up regulated features.
#' @param downColor Color used to plot enrichments in down regulated features.
#' @param pCutoff Adjusted P Value cutoff used to filter results before plotting.
#'
#' @return A ggplot containing the default ora plot.
#'
#' @export
#'
#' @import ggplot2
#'
defaultOraPlot <- function(x, idColumn = "comparison", splitStatus = FALSE, facetByStatus = TRUE,
                           upLabel = "Up", downLabel = "Down", idSep = "_",
                           upColor = "red", downColor = "blue", pCutoff = NULL) {

    # filter by cutoff before plotting
    if(!is.null(pCutoff)) {
        x <- x[x$p.adjust <= pCutoff,]
        # check non empty df
        if( nrow(x)==0 ) stop("No elements in the data frame after filtering by cutoff.")
    }
    # remove the up and down label from idColumn and add new column with status
    if(splitStatus) {
        x[,"plotStatus"] <- upLabel
        x[,"plotStatus"][grepl(downLabel, x[,idColumn])] <- downLabel
        pattern <- paste0(idSep, upLabel, "|", idSep, downLabel)
        x[, idColumn] <- gsub(pattern = pattern, replacement = "", x = x[, idColumn])
        # prepare color scale for plot
        colVal <- c(upColor, downColor)
        names(colVal) <- c(upLabel, downLabel)
    }
    # create plot
    p <- ggplot(data = x, mapping = aes(x = !!sym(idColumn), y = !!sym("ID"), size = !!sym("Count")))
    # add colored points if splitUpDown = TRUE
    if(splitStatus)
        p <- p +
        geom_point(mapping = aes(color = !!sym("plotStatus"))) +
        scale_color_manual(values = colVal)
    # facet by status if required
    if(splitStatus & facetByStatus)
        p <- p + facet_grid(cols = vars(!!sym("plotStatus")))
    # otherwise, add points colored by adjusted P
    if(!splitStatus)
        p <- p + geom_point(aes(color = !!sym("p.adjust")))
    # format final plot
    p <- p + theme_bw() + ggtitle("Over Representation Analysis") + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), axis.title.x = element_blank())
    return(p)

}

#' Create default GSEA dot plot
#'
#' Creates a simple representation of the results of \link[clusterProfiler]{GSEA} formatted
#' with \link[biokit]{cpResultsToDf}. The X axis represents comparison or status of the pathway.
#' The data frame can be filtered beforehand or setting the pCutoff parameter.
#'
#' @param x Data frame with results to plot.
#' @param splitById Use a column id to separate results?
#' @param idCol Name of the column used to separate results. Placed at the X axis instead of Up and down features.
#' @param upLabel Label used to indicate status of up regulated functional categories.
#' @param downLabel Label used to indicate status of down regulated functional categories.
#' @param upColor Color used to indicate status of up regulated functional categories.
#' @param downColor Color used to indicate status of down regulated functional categories.
#' @param pCutoff Adjusted P Value cutoff used to filter results before plotting.
#'
#' @return A ggplot containing the default gsea dot plot.
#'
#' @export
#'
#' @import ggplot2
#'
defaultGseaDotPlot <- function(x, splitById = FALSE, idCol = "comparison", upLabel = "Up",
                               downLabel = "Down", upColor = "red", downColor = "blue", pCutoff = NULL) {

    # filter by cutoff before plotting
    if(!is.null(pCutoff)) {
        x <- x[x$p.adjust <= pCutoff,]
        # check non empty df
        if( nrow(x)==0 ) stop("No elements in the data frame after filtering by cutoff.")
    }
    # add status column based on NES
    x[, "plotStatus"] <- upLabel
    x[, "plotStatus"][x[, "NES"] < 0 ] <- downLabel
    # prepare color scale for plot
    colVal <- c(upColor, downColor)
    names(colVal) <- c(upLabel, downLabel)
    shapeVal <- c(24,25)
    names(shapeVal) <- c(upLabel, downLabel)
    # create base plot depending on ID
    if(splitById) {
        p <- ggplot(data = x, mapping = aes( x = !!sym(idCol)))
    } else {
        p <- ggplot(data = x, mapping = aes( x = !!sym("plotStatus")))
    }
    # add points
    p <- p + geom_point(aes(y = !!sym("ID"), fill = !!sym("plotStatus"),
                            shape = !!sym("plotStatus"), size = -log10(!!sym("p.adjust")))) +
        scale_color_manual(values = colVal) + scale_shape_manual(values = shapeVal)
    # add theme and return
    p <- p +  ggtitle("Gene Set Enrichment Analysis")+ theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), axis.title.x = element_blank())
    return(p)

}




