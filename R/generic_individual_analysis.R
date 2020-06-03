#' Generic single matrix analysis
#'
#' Evaluates matrix data row-wise. This function is useful when the matrix
#' consists on relative measures (i.e a matrix of log2 fold changes) or to analyze
#' other features such as normal distribution.
#'
#' @param mat Matrix containing data to be analyzed.
#' @param metricFun Function to calculate metric of change. Defaults to log2(ratio of means).
#' @param statFun Function to calculate statistic of change. Defaults to a non-sensitive T-Test T value.
#' @param sigFun Function to calculate significance of change. Defaults to a non-sensitive T-Test P value.
#' @param pAdjustMethod Method used to adjust P values.
#' @param allowNa Allow NAs on input matrix?
#' @param featName Name of the column containing the compared features (rows).
#' @param metricName Name of the column containing the metric.
#' @param statName Name of the column containing the statistic.
#' @param pValName Name of the column containing the p values.
#' @param pAdjName Name of the column containing the adjusted p values.
#'
#' @return A data frame with the row-wise statistical results.
#'
#' @export
#'
singleMatrixTest <- function(mat, metricFun = function(x) mean(x, na.rm = TRUE),
                             statFun = nsTestT, sigFun = nsTestPValue,
                             allowNa = TRUE, pAdjustMethod = "BH", featName = "feature", metricName = "logFc",
                             statName = "t", pValName = "pValue", pAdjName = "pAdj") {

  # validate input matrix
  checkInput <- validateMatrix(mat = mat, checkRowNames = TRUE, allowNa = allowNa)
  # create statistics
  metric <- apply(mat, 1, metricFun)
  statistic <- apply(mat, 1, statFun)
  pValue <- apply(mat, 1, sigFun)
  pAdj <-stats:: p.adjust(p = pValue, method = pAdjustMethod)
  outId <- rownames(mat)
  # prepare out df
  outDf <- data.frame(outId, metric, statistic, pValue, pAdj)
  colnames(outDf) <- c(featName, metricName, statName, pValName, pAdjName)
  return(outDf)

}
