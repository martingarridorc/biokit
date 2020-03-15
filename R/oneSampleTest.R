
#' One sample T-Test over log ratio matrix
#'
#' Applies One Sample T-Test at row level in a matrix containing log transformed fold changes.
#'
#' @param x logFC matrix.
#' @param adjustMethod P value adjustment method.
#' @param idName Name for the resulting feature id column.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#' @param ... Rest of arguments passed to \link[stats]{t.test} function.
#'
#' @return A data frame with:
#' \item{id}{From row names.}
#' \item{logFC}{From row means.}
#' \item{pValue}{From the One Sample T-Test.}
#' \item{pAdj}{From multiple testing correction.}
#'
#' @export
#'
#' @importFrom stats p.adjust
#'
osTestMatrix <- function(x, adjustMethod = "BH", idName = "id", fcName = "logFc", pName = "pValue", pAdjName = "pAdj", ...) {

  # create statistics
  logFc <- rowMeans(x, na.rm = TRUE)
  pValue <- apply(x, 1, function(y) nsTest(y, ...))
  pAdj <- p.adjust(p = pValue, method = adjustMethod)
  # prepare id from rownames
  if (is.null(rownames(x)))
    outId <- 1:nrow(x)
  if (!is.null(rownames(x)))
    outId <- rownames(x)
  # prepare out df
  outDf <- data.frame(outId, logFc, pValue, pAdj, stringsAsFactors = FALSE)
  colnames(outDf) <- c(idName, fcName, pName, pAdjName)
  return(outDf)

}
