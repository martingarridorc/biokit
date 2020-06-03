#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep Separator of resulting contrasts.
#' @param reverseLevs Reverse levels before contrasts? Defaults to TRUE
#'
#' @return A character vector containing all the possible pairwise contrasts.
#'
#' @export
#'
#' @importFrom utils combn
#'
createPairwiseContrasts <- function(x, sep = "-", reverseLevs = TRUE) {

  # subset to unique elements of the vector
  difNames <- unique(x)
  if(length(difNames) != length(x)) {
    message("Found not unique levels on input vectors when making contrasts. They will be removed.")
  }
  # check NAs
  if(any(is.na(difNames))) {
    difNames <- difNames[!is.na(difNames)]
    message("Found NA on input vector when making contrasts. They will be removed.")
  }
  # reverse levels
  if(reverseLevs) {
    difNames <- rev(difNames)
  }
  # create all possible pairwise comparison matrix
  pairMat <- combn(difNames, 2)
  # return collapsed columns
  contrasts <- apply(pairMat, 2, function(y) paste(y[1], y[2], sep = sep))
  return(contrasts)

}

#' Create contrast matrix from design matrix.
#'
#' Uses the levels defined in the column names of the design matrix
#' to create a contrast matrix with \link[limma]{makeContrasts}.
#'
#' @param desMat The design matrix.
#'
#' @return A contrast matrix with all the possible pairwise comparisons.
#'
#' @export
#'
#' @importFrom limma makeContrasts
#'
contrastsFromDesign <- function(desMat) {

  # get contrasts from interesting column and build contrast matrix with limma
  contrasts <- createPairwiseContrasts(colnames(desMat))
  contrastMatrix <- limma::makeContrasts(contrasts = contrasts, levels = desMat)
  return(contrastMatrix)

}
