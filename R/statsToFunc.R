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
