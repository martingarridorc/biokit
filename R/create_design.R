#' Create design matrix from sample metadata
#'
#' Uses a data frame containing sample metadata to create a design matrix
#' for subsequent analyses.
#'
#' @param sampInfo The data frame that contains sample metadata.
#' @param groupCol The sampInfo column that is used to create the design matrix.
#' @param removeColName Remove column name from design matrix?
#'
#' @return The design matrix.
#'
#' @export
#'
#' @importFrom stats formula model.matrix
#'
designFromSampInfo <- function(sampInfo, groupCol, removeColName = TRUE) {

  # prepare formula
  groupFormula <- stats::formula(paste0("~ 0 + ", groupCol))
  # create design matrix
  designMatrix <- stats::model.matrix(groupFormula, sampInfo)
  # remove column name in design colnames
  if(removeColName) {
    colnames(designMatrix) <- gsub(pattern = groupCol, replacement = "", x = colnames(designMatrix))
  }
  return(designMatrix)

}
