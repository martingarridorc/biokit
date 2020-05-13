#' Create design matrix from sample metadata
#'
#' Uses a data frame containing sample metadata to create a design matrix
#' for subsequent analyses.
#'
#' @param sampInfo The data frame that contains sample metadata.
#' @param column The sampInfo column that is used to create the design matrix.
#' @param checkNames Check if vector contains syntactically valid names? If TRUE, it will fix not valid names.
#' @param removeColName Remove column name from design matrix?
#'
#' @return The design matrix.
#'
#' @export
#'
#' @importFrom stats formula model.matrix
#'
designFromSampInfo <- function(sampInfo, column, checkNames = TRUE, removeColName = TRUE) {

  intCol <- sampInfo[ , column]
  # check not syntactically valid names
  if(checkNames) {
    goodNames <- make.names(intCol)
    if(any(goodNames != intCol)) {
      message("Some elements are not syntactically valid names. They will be fixed.")
      sampInfo[ , column] <- goodNames
    }
  }
  # prepare formula
  groupFormula <- stats::formula(paste0("~ 0 + ", column))
  # create design matrix
  designMatrix <- stats::model.matrix(groupFormula, sampInfo)
  # remove column name in design colnames
  if(removeColName) {
    colnames(designMatrix) <- gsub(pattern = column, replacement = "", x = colnames(designMatrix))
  }
  return(designMatrix)

}
