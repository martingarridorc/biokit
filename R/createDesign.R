#' Create design matrix from sample metadata
#'
#' Uses a data frame containing sample metadata to create a design matrix
#' for subsequent analyses.
#'
#' @param x The data frame containing the column metadata (e.g sample grouping).
#' @param column The data frame column that is used to compose the design matrix.
#'
#' @return The composed design matrix.
#'
#' @export
#'
#' @importFrom stats formula model.matrix
#'
designFromSampInfo <- function(x, column) {

  # prepare formula
  groupFormula <- formula(paste0("~ 0 + ", column))
  # create design matrix
  designMatrix <- model.matrix(groupFormula, x)
  # remove column name in design colnames
  colnames(designMatrix) <- gsub(pattern = column, replacement = "", x = colnames(designMatrix))
  return(designMatrix)

}
