#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep The separator for the resulting contrasts.
#'
#' @return A character vector containing all the possible pairwise contrasts.
#' @export
#'
pairwiseContrasts <- function(x, sep = "-") {

  # create all possible pairwise comparison matrix
  pairMat <- combn(unique(x), 2)

  # return collapsed columns
  contrasts <- apply(pairMat, 2, function(y) paste(y[1], y[2], sep = sep))

  return(contrasts)

}

#' Perform PCA and create list with results
#'
#' Performs principal component analysis for a given matrix
#' and returns a list with the formatted results.
#'
#'
#' @param m The matrix to analyze.
#' @param transpose Wether to transpose or not the input matrix.
#' @param roundDigits Number of digits for the rounded percent of
#' variance explained by each PC.
#' @param ... Other parameters passed to the prcomp function
#'
#' @return A list with three elements:
#' \item{result}{The prcomp object containing the results of the PCA.}
#' \item{summary}{The summary objects containing the matrix with PC values.}
#' \item{pcts}{A character vector with the variance explained by each PC.}
#' @export
#'
pcaToList <- function(m, transpose = FALSE, roundDigits = 2, ...) {

  # transpose if it is not transposed
  if(transpose) m <- t(m)

  # perform prcomp analysis
  pcaRes <- prcomp(m , ...)
  pcaResSum <- summary(pcaRes)

  # add proportion of variance explained
  variancePcts <- round(pcaResSum$importance[2, ]*100, roundDigits)
  pcPcts <- paste0(colnames(pcaRes$x), " ( ", variancePcts, "% )")

  return(list(result = pcaRes, summary = pcaResSum, pcts = pcPcts))

}
