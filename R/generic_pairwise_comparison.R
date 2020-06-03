#' Generic pairwise matrix test
#'
#' Compare two pieces of data from a matrix given a
#' sample information data frame and a grouping column.
#'
#' @param mat Matrix containing data to be analyzed.
#' @param sampInfo Sample information data frame.
#' @param groupCol Name of the sampInfo column containing the grouping variable.
#' @param numerator Label of the grouping variable to be used as numerator.
#' @param denominator Label of the grouping variable to be used as denominator.
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
pairwiseMatrixTest <- function(mat, sampInfo, groupCol, numerator, denominator,
                               metricFun = defaultLogFc, statFun = nsTestT, sigFun = nsTestPValue,
                               pAdjustMethod = "BH", allowNa = TRUE,
                               featName = "feature", metricName = "logFc", statName = "t",
                               pValName = "pValue", pAdjName = "pAdj") {

  # validate input matrix
  checkInput <- validateMatrix(mat = mat, checkRowNames = TRUE, allowNa = allowNa)
  # get numerator and denominator indexs
  numIndex <- sampInfo[, groupCol] == numerator
  denIndex <- sampInfo[, groupCol] == denominator
  # create pairwise statistics
  metric <- apply(mat, 1, function(rData) metricFun(x = rData[numIndex], y = rData[denIndex]))
  statistic <- apply(mat, 1, function(rData) statFun(x = rData[numIndex], y = rData[denIndex]))
  pValue <- apply(mat, 1, function(rData) sigFun(x = rData[numIndex], y = rData[denIndex]))
  pAdj <- stats::p.adjust(p = pValue, method = pAdjustMethod)
  # prepare id from rownames
  outId <- rownames(mat)
  # prepare out df
  outDf <- data.frame(outId, metric, statistic, pValue, pAdj)
  colnames(outDf) <- c(featName, metricName, statName, pValName, pAdjName)
  return(outDf)

}

#' Automatic pairwise contrast
#'
#' Perform the automatic pairwise comparison
#' of matrix data using a sample information
#' data frame column.
#'
#' @param mat Matrix containing data to be analyzed.
#' @param sampInfo Sample information data frame.
#' @param groupCol Name of the sampInfo column containing the grouping variable.
#' @param sep Separator while creating contrasts.
#' @param compName Name of the output column that identifies the comparison.
#' @param ... Other arguments for \link[biokit]{pairwiseMatrixTest}.
#'
#' @return A data frame containing the results for all the comparisons
#'
#' @export
#'
#' @importFrom dplyr bind_rows %>%
#'
autoPairwiseMatrixTest <- function(mat, sampInfo, groupCol, sep = "-",
                                   compName = "comparison", ...) {

  # validate sample information fixing broken names
  sampInfo <- validateSampInfo(sampInfo = sampInfo, groupCol = groupCol, mat = mat, checkNames = TRUE)
  # create contrasts
  levs <- as.character(sampInfo[ , groupCol]) %>%
    unique()
  contrasts <- createPairwiseContrasts(levs, sep = sep)
  # obtain results by contrast
  dfList <- lapply(contrasts, function(c) {
    num <- strsplit(x = c, split = sep)[[1]][1]
    den <- strsplit(x = c, split = sep)[[1]][2]
    resDf <- pairwiseMatrixTest(mat = mat, sampInfo = sampInfo, groupCol = groupCol,
                                numerator = num, denominator = den, ...)
    return(resDf)
  })
  names(dfList) <- contrasts
  # bind rows of df list
  outDf <- dplyr::bind_rows(dfList, .id = compName)
  return(outDf)

}
