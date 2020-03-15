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
