#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep Separator for the resulting contrasts.
#'
#' @return A character vector containing all the possible pairwise contrasts.
#'
#' @export
#'
#' @importFrom utils combn
#'
pairwiseContrasts <- function(x, sep = "-") {

  # create all possible pairwise comparison matrix
  pairMat <- combn(unique(x), 2)
  # return collapsed columns
  contrasts <- apply(pairMat, 2, function(y) paste(y[1], y[2], sep = sep))
  return(contrasts)

}

#' Create contrast matrix given a design matrix.
#'
#' Uses levels defined in the column names of the design matrix
#' to create a contrast matrix with \link[limma]{makeContrasts}.
#'
#' @param x The design matrix.
#'
#' @return A contrast matrix with all the possible pairwise comparisons.
#' @export
#'
#' @importFrom limma makeContrasts
#'
contrastsFromDesign <- function(x) {

  # get contrasts from interesting column and build contrast matrix with limma
  contrasts <- pairwiseContrasts(colnames(x))
  contrastMatrix <- limma::makeContrasts(contrasts = contrasts, levels = x)
  return(contrastMatrix)

}
