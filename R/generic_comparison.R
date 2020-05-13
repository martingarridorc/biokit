compareMatrixPairwise <- function(mat, sampInfo, groupCol, numerator, denominator,
                                  metricFun, statFun, pValueFun, pAdjustMethod = "BH") {

  # get numerator and denominator indexs
  numIndex <- rownames(sampInfo)[sampInfo[, groupCol] == numerator]
  denIndex <- rownames(sampInfo)[sampInfo[, groupCol] == denominator]
  # create pairwise statistics
  metric <- apply(x, 1, function(r) metricFun(x = r[numIndex], y = r[denIndex]))
  statistic <- apply(x, 1, function(r) statisticFun(x = r[numIndex], y = r[denIndex]))
  pValue <- apply(x, 1, function(r) pValueFun(x = r[numIndex], y = r[denIndex]))
  pAdj <- p.adjust(p = pValue, method = adjustMethod)
  # prepare id from rownames
  if (is.null(rownames(x)))
    outId <- 1:nrow(x)
  if (!is.null(rownames(x)))
    outId <- rownames(x)
  # prepare out df
  outDf <- data.frame(outId, metric, statistic, pValue, pAdj, stringsAsFactors = FALSE)
  colnames(outDf) <- c(featName, metricName, statName, pName, pAdjName)
  return(outDf)

}


