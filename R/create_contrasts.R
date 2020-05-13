#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep Separator of resulting contrasts.
#'
#' @return A character vector containing all the possible pairwise contrasts.
#'
#' @export
#'
#' @importFrom utils combn
#'
pairwiseContrasts <- function(x, sep = "-") {

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
  # create all possible pairwise comparison matrix
  pairMat <- combn(difNames, 2)
  # return collapsed columns
  contrasts <- apply(pairMat, 2, function(y) paste(y[1], y[2], sep = sep))
  return(contrasts)

}
