#' Row level statistics of matrix
#'
#' Applies several functions over the rows of a matrix to obtain a tidy
#' data frame containing the results. It uses three functions to obtain: 1: a metric,
#' which defaults to \link[stat]{mean}, 2: a statistic, which defaults to
#' \link[biokit]{nsTestT} and 3: a p value, which defaults to \link[stats]{nsTestPValue}.
#'
#' @param x Input matrix.
#' @param metricFun Function to obtain metric.
#' @param statisticFun Function to obtain statistic.
#' @param pValueFun Function to obtain p-value.
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
matrixSingleTest <- function(x,
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

#' Row level statistics from the pairwise comparison of matrix columns
#'
#' Divides the matrix into numerator and denominator using a column information data frame (sampInfoDf)
#' and applies several functions at row level to obtain a tidy data frame containing the results. It uses three functions to obtain:
#' 1: a metric, which defaults to \link[stat]{mean} division, 2: a statistic, which defaults to
#' \link[biokit]{nsTestT} and 3: a p value, which defaults to \link[stats]{nsTestPValue}.
#'
#' @param x Input matrix.
#' @param sampInfoDf Column information data frame. Rownames of this data frame should match matrix columns.
#' @param groupCol Column from the sampInfoDf that indicates groups of columns (samples?).
#' @param numerator Group of columns used as the numerator (first element) of functions.
#' @param denominato rGroup of columns used as the denominator (second element) of functions.
#' @param metricFun Function to obtain metric.
#' @param statisticFun Function to obtain statistic.
#' @param pValueFun Function to obtain p-value.
#' @param adjustMethod Method to perform p value adjustment.
#' @param idName Name for the id column, which is created using \link[base]{rownames}.
#' @param metricName Name for the metric column.
#' @param statName Name for the statistic column.
#' @param pName Name for the p value column.
#' @param pAdjName Name for the adjusted p value column.
#'
#' @return  A tidy data frame with results.
#'
#' @export
#'
matrixPairwiseTest <- function(x, sampInfoDf, groupCol, numerator, denominator,
                               metricFun = function(x,y) log2(mean(x, na.rm = TRUE) / mean(y, na.rm = TRUE)),
                               statisticFun = function(x,y) nsTestT(x = x, y = y),
                               pValueFun = function(x,y) nsTestPValue(x = x, y = y),
                               adjustMethod = "BH",
                               idName = "id", metricName = "logFc",statName = "t", pName = "pValue", pAdjName = "pAdj") {

  # check that sampInfoDf rownames match colnames of the matrix to analyse
  if (!all(rownames(sampInfoDf) == colnames(x)))
    stop("Rownames of the sample information data frame could not be mapped to matrix columns.")
  # get numerator and denominator indexs
  numIndex <- rownames(sampInfoDf)[sampInfoDf[,groupCol] == numerator]
  denIndex <- rownames(sampInfoDf)[sampInfoDf[,groupCol] == denominator]
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
  colnames(outDf) <- c(idName, metricName, statName, pName, pAdjName)
  return(outDf)

}

#' Pairwise matrix test from contrasts
#'
#' Performs the pairwise statistical analysis of the groups of columns of interest
#' using an input vector of contrasts.
#'
#' @param x Input matrix.
#' @param sampInfoDf Column information data frame. Rownames of this data frame should match matrix columns.
#' @param groupCol Column from the sampInfoDf that indicates groups of columns.
#' @param contrasts Vector containing the contrasts to analyze.
#' @param sep Character used to split contrasts into numerator and denominator.
#' @param compName Name for the column that indicates the pairwise comparison in the tidy data frame.
#' @param ... Rest of arguments passed to \link[biokit]{matrixPairwiseTest}.
#'
#' @return  A tidy data frame with results.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#'
matrixTestFromContrasts <- function(x, sampInfoDf, groupCol, contrasts, sep = "-", compName = "comparison", ...) {

  # iterate through contrasts, using first element of the splitted contrast as numerator and the second as denominator
  statList <- lapply(contrasts, function(y) {
    splitted <- unlist(strsplit(y, split = sep))
    # check length of splitted not higher of two
    if(length(splitted) != 2)
      stop("Contrasts should be splitted in only two elements when using the separator character")
    num <- splitted[1]
    den <- splitted[2]
    df <- matrixPairwiseTest(x = x, sampInfoDf = sampInfoDf, groupCol = groupCol, numerator = num, denominator = den, ...)
    return(df)
  })
  # set list names
  names(statList) <- contrasts
  # bind rows of resulting list
  result <- dplyr::bind_rows(statList, .id = compName)
  return(result)

}

#' Automatic pairwise comparison of SummarizedExperiment matrix
#'
#' Uses a matrix stored in a SummarizedExperiment object to
#' automatically perform all the possible pairwise comparisons
#' between the groups defined in a \link[SummarizedExperiment]{colData} column.
#'
#' @param se The SummarizedExperiment object to analyze.
#' @param groupCol The id of the \link[SummarizedExperiment]{colData} column that is used to create the pairwise contrasts.
#' @param ... Rest of arguments passed to \link[biokit]{matrixTestFromContrasts}.
#'
#' @return  A tidy data frame with results.
#'
#' @export
#'
autoMatrixPairwise <- function(se, groupCol, ...) {

  contrasts <- pairwiseContrasts(colData(se)[, groupCol])
  mat <- assay(se)
  res <- matrixTestFromContrasts(x = mat, sampInfoDf = colData(se), contrasts = contrasts, groupCol = groupCol, ...)
  return(res)

}












