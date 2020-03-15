#' Row level statistics of matrix
#'
#' Applies several functions over the rows of a matrix to obtain a tidy
#' data frame containing the results. Defaults to a One-Sample T-Test over
#' the input matrix
#'
#' @param x Matrix to analyze.
#' @param metricFun Function to obtain the metric.
#' @param statisticFun Function to obtain the statistic.
#' @param pValueFun Function to obtain significance as p-value.
#' @param adjustMethod Method to perform p value adjustment.
#' @param idName Name for the id column, which is created using \link[base]{rownames}.
#' @param metricName Name for the metric column.
#' @param statName Name for the statistic column.
#' @param pName Name for the p value column.
#' @param pAdjName Name for the adjusted p value column.
#'
#' @return A tidy data frame with results.
#'
#' @export
#'
singleMatrixTest <- function(x,
                         metricFun = function(x) mean(x, na.rm = TRUE),
                         statisticFun = function(x) nsTestT(x),
                         pValueFun = function(x) nsTestPValue(x),
                         adjustMethod = "BH",
                         idName = "id", metricName = "logFc",statName = "t", pName = "pValue", pAdjName = "pAdj") {

  # create statistics
  metric <- apply(x, 1, metricFun)
  statistic <- apply(x, 1, statisticFun)
  pValue <- apply(x, 1, pValueFun)
  pAdj <- p.adjust(p = pValue, method = adjustMethod)
  # prepare id from rownames
  if (is.null(rownames(x)))
    outId <- 1:nrow(x)
  if (!is.null(rownames(x)))
    outId <- rownames(x)
  # prepare out df
  outDf <- data.frame(outId, metric, statistic, pValue, pAdj, stringsAsFactors = FALSE)
  colnames(outDf) <- c(idName, metricName, statName, pName, pAdjName)
  return(outDf)

}
