#' Add status by cutoff
#'
#' Adds a "status" column to a data frame using a specific cutoff.
#'
#' @param resDf The data frame to annotate.
#' @param splitUpDown Split annotation into up and down features? If TRUE, splitCol and splitMidPoint should be provided.
#' @param splitCol If splitUpDown = TRUE, the column used to split features into up and down.
#' @param splitMidPoint If splitUpDown = TRUE, the splitCol value considered as midpoint.
#' @param metricCol Column containing a metric to stablish the cutoff.
#' @param metricCutoff The metric cutoff value. Rows with an absolute metricCol value above it will be annotated.
#' @param sigCol Column containing a significance value to stablish the cutoff.
#' @param sigCutoff The significance cutoff value. Rows with an sigCol value below it will be annotated.
#' @param noChangeLabel Label to be used in not annotated rows.
#' @param statusCol Name of the status column.
#' @param upLabel Label for up-regulated features.
#' @param downLabel Label for down-regulated features.
#' @param sigLabel Label for significant features. Only used if splitUpDown = FALSE.
#'
#' @return The annotated data frame.
#'
#' @export
#'
addStatusByCutoff <- function(resDf, splitUpDown,
                              splitCol = "logFc", splitMidPoint = 0,
                              metricCol = "logFc", metricCutoff = 1,
                              sigCol = "pAdj", sigCutoff = 0.05,
                              noChangeLabel = "No change", statusCol = "status",
                              upLabel = "Up", downLabel = "Down",
                              sigLabel = "Significant") {

  # add no change label column
  resDf[ , statusCol] <- noChangeLabel
  if(splitUpDown) {
    upIndex <- resDf[ , splitCol] >= splitMidPoint & abs(resDf[, metricCol]) >= metricCutoff & resDf[, sigCol] <= sigCutoff
    downIndex <- resDf[ , splitCol] < splitMidPoint & abs(resDf[, metricCol]) >= metricCutoff & resDf[, sigCol] <= sigCutoff
    resDf[upIndex, statusCol] <- upLabel
    resDf[downIndex, statusCol] <- downLabel
  } else {
    sigIndex <- abs(resDf[, metricCol]) >= metricCutoff & resDf[, sigCol] <= sigCutoff
    resDf[sigIndex, statusCol] <- sigLabel
  }
  return(resDf)

}

#' Add status by rank
#'
#' Order the input data frame using a specific column and annotates
#' the top N features.
#'
#' @param resDf The data frame to annotate.
#' @param rankCol The column used to order the data frame.
#' @param splitUpDown Split annotation into up and down features? If TRUE, splitCol and splitMidPoint should be provided.
#' @param splitCol If splitUpDown = TRUE, the column used to split features into up and down.
#' @param splitMidPoint If splitUpDown = TRUE, the splitCol value considered as midpoint.
#' @param topN The top N features to be annotated.
#' @param decreasing Order rankCol decreasing?
#' @param absolute Use absolute rankCol values?
#' @param noChangeLabel Label to be used in not annotated rows.
#' @param upLabel Label for up-regulated features.
#' @param downLabel Label for down-regulated features.
#' @param sigLabel Label for significant features. Only used if splitUpDown = FALSE.
#' @param statusCol Name of the status column.
#'
#' @return The annotated data frame.
#'
#' @export
#'
addStatusByRank <- function(resDf, rankCol, topN, splitUpDown,
                            decreasing = TRUE, absolute = TRUE,
                            splitCol = "logFc", splitMidPoint = 0,
                            noChangeLabel = "No change", statusCol = "status",
                            upLabel = "Up", downLabel = "Down",
                            sigLabel = "Significant") {

  # add no change label column
  resDf[ , statusCol] <- noChangeLabel
  # re order input data frame
  if(absolute) {
    newOrder <- order(abs(resDf[ , rankCol]), decreasing = decreasing)
  } else {
    newOrder <- order(resDf[ , rankCol] , decreasing = decreasing)
  }
  resDf <- resDf[newOrder, ]
  # add status column
  if(splitUpDown) {
    upIndex <- resDf[, splitCol] >= splitMidPoint
    upIndex <- which(upIndex)[1:topN]
    downIndex <- resDf[, splitCol] < splitMidPoint
    downIndex <- which(downIndex)[1:topN]
    resDf[upIndex, statusCol] <- upLabel
    resDf[downIndex, statusCol] <- downLabel
  } else {
    resDf[1:topN, statusCol] <- sigLabel
  }
  return(resDf)

}


