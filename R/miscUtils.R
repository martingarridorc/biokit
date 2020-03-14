#' Add status column to data frame using cutoff
#'
#' Uses a significance and metric cutoff to annotate a
#' data frame containing the results of a stat analysis.
#'
#' @param x Data frame to annotate.
#' @param metricCutoff Metric cutoff. If set to \code{NULL}, only the significance cutoff is used.
#' @param sigCutoff Significance cutoff. If set to \code{NULL}, only the metric cutoff is used.
#' @param splitUpDown Split generated label into up and down features? If TRUE, uses \code{-metricCutoff} to select down features.
#' @param metricCol Column used to filter by metric.
#' @param sigCol Column used to filter by significance.
#' @param noChangeLabel Label used for not-changing features.
#' @param annotationCol Name of the new column containing labels.
#' @param upLabel Label for up-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param downLabel Label for down-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param sigLabel Label for significantly altered features. Only used if \code{splitUpDown = FALSE}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
annotateByCutoff <- function(x, metricCutoff = 1, sigCutoff = 0.05, splitUpDown = TRUE, metricCol = "logFc", sigCol = "pAdj",
    noChangeLabel = "No change", annotationCol = "status", upLabel = "Up", downLabel = "Down", sigLabel = "Significant") {

    x[, annotationCol] <- noChangeLabel
    # if split, then return up and down labels, else return significant or not significant
    if (splitUpDown) {
        if(is.null(sigCutoff)) {
            up <- x[, metricCol] >= metricCutoff
            down <- x[, metricCol] <= -metricCutoff
            x[up, annotationCol] <- upLabel
            x[down, annotationCol] <- downLabel
        } else {
            up <- x[, metricCol] >= metricCutoff & x[, sigCol] <= sigCutoff
            down <- x[, metricCol] <= -metricCutoff & x[, sigCol] <= sigCutoff
            x[up, annotationCol] <- upLabel
            x[down, annotationCol] <- downLabel
        }
    } else {
        if (is.null(metricCutoff))
            sig <- x[, sigCol] <= sigCutoff
        if (is.null(sigCutoff))
            sig <- abs(x[, metricCol]) >= metricCutoff
        if (!is.null(sigCutoff) & !is.null(metricCutoff))
            sig <- abs(x[, metricCol]) >= metricCutoff & x[, sigCol] <= sigCutoff
        x[sig, annotationCol] <- sigLabel
    }
    return(x)

}

#' Add status column to comparison data frame using top N features
#'
#' Orders a data frame and adds a column to indicate the
#' status of the top N features.
#'
#' @param x Data frame to annotate.
#' @param n Top N features to be annotated
#' @param sortCol Column used to sort the data frame.
#' @param decreasing Order sortCol decreasing?
#' @param splitUpDown Split top N features into up and down regulated? Uses the metricCol argument. This will annotate 2*n features.
#' @param twoSides Pick up both extremes of the sorted data frame as the top N features? SplitUpDown should be TRUE.
#' @param metricCol Metric column used to sepparate up and down features.
#' @param metricCutoff Cutoff used to separate features into up and down. Only used if \code{splitUpDown = TRUE}.
#' @param noChangeLabel Label used for not-changing features.
#' @param annotationCol Name of the new column containing labels.
#' @param upLabel Label for up-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param downLabel Label for down-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param sigLabel Label for significantly altered features. Only used if \code{splitUpDown = FALSE}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
annotateTopN <- function(x, n, sortCol, decreasing = FALSE, splitUpDown = TRUE, twoSides = FALSE, metricCol = "logFc", metricCutoff = 0,
    noChangeLabel = "No change", annotationCol = "status", upLabel = "Up", downLabel = "Down", sigLabel = "Significant") {

    x[, annotationCol] <- noChangeLabel
    x <- x[order(x[, sortCol], decreasing = decreasing), ]
    # if split and twoSides, annotate both extremes of the df
    if (splitUpDown & twoSides) {
        x[1:n, annotationCol] <- upLabel
        x[(nrow(x) - n + 1):nrow(x), annotationCol] <- downLabel
    }
    # if split and !twoSides, pick top N
    if (splitUpDown & !twoSides) {
        x[x[, metricCol] > metricCutoff, annotationCol][1:n] <- upLabel
        x[x[, metricCol] < metricCutoff, annotationCol][1:n] <- downLabel
    }
    # else return top N
    if (!splitUpDown & !twoSides) {
        x[, annotationCol][1:n] <- sigLabel
    }
    return(x)

}

#' Annotate multiple comparisons data frame
#'
#' Uses \link[biokit]{annotateByCutoff} or \link[biokit]{annotateTopN} to annotate a
#' data frame by sepparating data for each comparison.
#'
#' @param x Multi-comparison data frame to annotate.
#' @param useCutoff Wether to use \link[biokit]{annotateByCutoff} (if TRUE) or \link[biokit]{annotateTopN} (if FALSE).
#' @param compCol Column with the comparison id.
#' @param ... Rest of arguments passed to \link[biokit]{annotateByCutoff} or \link[biokit]{annotateTopN}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#'
annotateMultiComparison <- function(x, useCutoff, compCol = "comparison", ...) {

    # split list into comparisons
    dfList <- split(x, f = x[, compCol])
    # apply function of interest
    if (useCutoff) {
        aDfList <- lapply(dfList, function(y) annotateByCutoff(x = y, ...))
    } else {
        aDfList <- lapply(dfList, function(y) annotateTopN(x = y, ...))
    }
    # coerce resultant data frames
    result <- dplyr::bind_rows(aDfList)
    return(result)

}


#' Split data frame into list of features
#'
#' Uses annotation labels to split a data frame of results into
#' a list containing a character vector of features. If an additional split column is provided,
#' then it split by the combination of the annotation and the split column.
#'
#' @param x Data frame to split.
#' @param annotCol Column containing the annotation label.
#' @param featCol Column containing the feature id. Coerced to character.
#' @param labelToRemove Label assigned to the elements which will not be splitted.
#' @param splitCol Additional column to use for splitting.
#' @param sep Separator of the new splitting variable.
#'
#' @return A list of character vectors.
#'
#' @export
#'
splitFeatures <- function(x, annotCol = "status", featCol = "feature", labelToRemove = "No change",
                          splitCol = NULL, sep = "_") {

    # remove rows based on annotation column and label to remove
    keep <- x[, annotCol] != labelToRemove
    x <- x[keep,]
    # if an additional split column is defined, then paste it to annotation column
    if(is.null(splitCol)) {
        splitted <- split(x, f = x[,annotCol])
        splitted <- lapply(splitted, function(y) as.character(y[, featCol]))
    } else {
        x$toSplitBy <- paste(x[,splitCol], x[,annotCol], sep = "_")
        splitted <- split(x, f = x$toSplitBy)
        splitted <- lapply(splitted, function(y) as.character(y[, featCol]))
    }
    return(splitted)

}

#' Get ranked vector from data frame
#'
#' Prepares a vector of the given data frame column, after naming and sorting it.
#'
#' @param x Data frame to process.
#' @param metricCol Column to extract as ranked vector.
#' @param featCol Column containing the names for the ranked vector.
#' @param decreasing Sort vector decreasing?
#'
#' @return A vector in the desired format.
#'
#' @export
#'
getRankedVector <- function(x, metricCol = "logFc", featCol = "feature", decreasing = TRUE) {

    ranked <- x[, metricCol]
    names(ranked) <- x[, featCol]
    ranked <- sort(ranked, decreasing = decreasing)
    return(ranked)

}

#' Get ranked list of vectors from data frame
#'
#' Uses a column to split a data frame and then creates a list of named vectors using
#' \link[biokit]{getRankedVector}.
#'
#' @param x Data frame to process.
#' @param splitCol Column used to split the data frame.
#' @param ... Rest of arguments passed to \link[biokit]{getRankedVector}.
#'
#' @return A list of sorted and named vectors.
#'
#' @export
#'
getRankedVectorList <- function(x, splitCol = "comparison", ...) {

    splitted <- split(x, f = x[,splitCol])
    rankedList <- lapply(splitted, function(y) getRankedVector(y, ...))
    return(rankedList)

}
