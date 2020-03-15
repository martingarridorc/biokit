#' Over-representation analysis over a list of features
#'
#' Uses a over-representation strategy (Fisher Exact T-Test) to
#' evaluate the over-represented categories in a list of features. Uses
#' the \link[clusterProfiler]{enricher} function to perform the analyses.
#'
#' @param x The list of features to be analysed.
#' @param funCategories Two-columns data frame containing the functional categories info.
#' First column should correspond to categories and the second to features.
#' Passed to \link[clusterProfiler]{enricher} TERM2GENE argument.
#' @param minFeatures Minimum features to evaluate.
#' The function will filter all list elements with a lower number of features than this value.
#' @param ... Rest of arguments passed to \link[clusterProfiler]{enricher}
#'
#' @return A list of \link[clusterProfiler]{enricher} results, one per list element.
#'
#' @export
#'
#' @importFrom clusterProfiler enricher
#'
oraFromList <- function(x, funCategories, minFeatures = 10, ...) {

  # remove list elements lower than minFeatures
  l <- unlist(lapply(x, length))
  keep <- unlist(l) >= minFeatures
  x <- x[keep]
  # perform enrichment analysis in all elements of the list
  result <- lapply(x, function(y) enricher(gene = y, TERM2GENE = funCategories, ...))
  return(result)

}

#' Perform Gene Set Enrichment Analysis from list of ranked vectors
#'
#' Uses \link[clusterProfiler]{GSEA} to perform the gene set enrichment analysis
#' over a list of ranked vectors.
#'
#' @param x The list of ranked vectors.
#' @param funCategories Two-columns data frame containing the functional categories info.
#' First column should correspond to categories and the second to features.
#' Passed to \link[clusterProfiler]{GSEA} TERM2GENE argument.
#' @param ... Rest of arguments passed to \link[clusterProfiler]{GSEA}.
#'
#' @return A list of \link[clusterProfiler]{GSEA} results.
#'
#' @export
#'
#' @importFrom clusterProfiler GSEA
#'
gseaFromList <- function(x, funCategories, ...) {

  result <- lapply(x, function(y) GSEA(geneList = y, TERM2GENE = funCategories, ...))
  return(result)

}

#' Transform list of clusterProfiler results into tidy data frame
#'
#' Binds the rows of a list of clusterProfiler results annotating the
#' resulting data frame with the list names.
#'
#' @param x The list of clusterProfiler results
#' @param splitName The name for the new id column.
#'
#' @return A data frame with the clusterProfiler results in a tidy format.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#'
cpResultsToDf <- function(x, splitName = "comparison") {

  dfList <- lapply(x, function(y)y@result)
  df <- dplyr::bind_rows(dfList, .id = splitName)
  return(df)

}
