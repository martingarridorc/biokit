#' Normalize counts with TMM
#'
#' Normalizes a raw counts matrix into log transformed counts per millions using
#' the edgeR 'Trimmed Mean of M-Values' (TMM) normalization method. For more information
#' please refer to the edgeR package (https://bioconductor.org/packages/release/bioc/html/edgeR.html).
#'
#' @param counts The matrix of counts.
#' @param log Passed to \link[edgeR]{cpm}.
#' @param prior.count Passed to \link[edgeR]{cpm}.
#'
#' @return The normalized expression matrix
#'
#' @export
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'
countsToTmm <- function(counts, log = TRUE, prior.count = 3) {

  # check matrix is ok with no NAs
  inputCheck <- validateMatrix(mat = counts, allowNa = FALSE)
  # transform into edgeR DGEList
  dgeList <- edgeR::DGEList(counts = counts)
  # calculate normalization factors
  dgeList <- edgeR::calcNormFactors(object = dgeList, method = "TMM")
  # obtain normalized logCpm matrix
  logCpm <- edgeR::cpm(dgeList, log = log, prior.count = prior.count)
  return(logCpm)

}
