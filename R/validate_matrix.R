#' Validate matrix
#'
#' Evaluates integrity of the input matrix.
#'
#' @param mat The matrix to validate.
#' @param checkRowNames Check rownames are present? Stop execution if TRUE and mat does not contain rownames.
#' @param allowNa Allow NA presence in mat?  Stop execution if NAs are found in mat.
#'
#' @return TRUE if no errors are found.
#'
validateMatrix <- function(mat, checkRowNames = TRUE, allowNa = TRUE) {

  # check class
  matClass <- class(mat)
  if(! all(matClass == c("matrix", "array"))) stop("Input is not a matrix.")
  # check rownames
  if(checkRowNames & is.null(rownames(mat))) stop("Matrix does not have rownames.")
  # check NAs
  naCount <- sum(is.na(mat))
  if(naCount != 0 & allowNa) message("Matrix contain NAs. This may produce unexpected results.")
  if(naCount != 0 & (! allowNa)) stop("Matrix contain NAs.")
  return(TRUE)

}
