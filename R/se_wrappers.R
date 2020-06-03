#' Translate SummarizedExperiment
#'
#' Translates the first assay of a SummarizedExperiment.
#'
#' @param se The SummarizedExperiment object that contains the matrix to be translated.
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return The translated version of the summarized experiment.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay
#'
translateSE <- function(se, df, sourceKey, targetKey, summariseFun) {

  if(class(se) != "SummarizedExperiment") {
    stop("se is not an object of the SummarizedExperiment class.")
  }
  # get original matrix and translate it
  originalMat <- SummarizedExperiment::assay(se)
  newMat <- translateMatrix(mat = originalMat, df = df, sourceKey = sourceKey,
                            targetKey = targetKey, summariseFun = summariseFun)
  SummarizedExperiment::assay(se) <- newMat
  return(se)

}

#' Translate SummarizedExperiment with annotation package
#'
#' @param se The SummarizedExperiment object that contains the matrix to be translated.
#' @param db Annotation package object to use for the translation.
#' @param sourceKey Source key in the DB.
#' @param targetKey Target key in the DB.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return The translated version of the summarized experiment.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay
#'
translateSEWithDb <- function(se, db, sourceKey, targetKey, summariseFun) {

  if(class(se) != "SummarizedExperiment") {
    stop("se is not an object of the SummarizedExperiment class.")
  }
  # get original matrix and translate it
  originalMat <- SummarizedExperiment::assay(se)
  newMat <- translateMatrixWithDb(mat = originalMat, db = db, sourceKey = sourceKey,
                                  targetKey = targetKey, summariseFun = summariseFun)
  SummarizedExperiment::assay(se) <- newMat
  return(se)

}
