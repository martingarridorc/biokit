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
