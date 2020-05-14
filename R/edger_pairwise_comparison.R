#' Obtain edgeR results from counts, design and contrasts matrix
#'
#' Applies recommended edgeR steps to obtain differential expression tables from
#' an initial count matrix object. As in \link[biokit]{limmaDfFromContrasts}, it uses
#' a design matrix and a contrast matrix as input.
#'
#' @param counts The count matrix to analyze.
#' @param desMat Design matrix.
#' @param conMat Contrasts matrix.
#' @param compName Name for the column that indicates the pairwise comparison in the tidy data frame.
#' @param featName Name for the column that indicates the analyzed row feature in the tidy data frame.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#' @param exprName Name for the resulting average expresion column.
#' @param filterByExpr Apply \link[edgeR]{filterByExpr}?
#'
#' @return A data frame with the results in a tidy format.
#'
#' @export
#'
#' @importFrom edgeR calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr %>% bind_rows
#' @importFrom tibble rownames_to_column
#'
edgerDfFromContrasts <- function(counts, desMat, conMat, filterByExpr = FALSE,
                                 compName = "comparison", featName = "feature",
                                 exprName = "AveExpr", fcName = "logFc",
                                 pName = "pValue", pAdjName = "pAdj") {

  # validate input matrix
  inputCheck <- validateMatrix(mat = counts, checkRowNames = TRUE, allowNa = FALSE)
  # create DGElist, calculate norm factor and estimate dispersion
  y <- DGEList(counts)
  # filter by expr if required
  if(filterByExpr) {
    keep <- edgeR::filterByExpr(y, design = desMat)
    y <- y[keep, ]
  }
  # calculate normalization factor and estimate dispersion
  y <- edgeR::calcNormFactors(object = y) %>%
    edgeR::estimateDisp(y = ., design = desMat)
  # fit ql netgative binomial gl model
  fit <- edgeR::glmQLFit(y, desMat)
  # obtain DE table for each contrast
  dfList <- apply(conMat, 2, function(c) {
    t <- edgeR::glmQLFTest(fit, contrast = c)
    df <- edgeR::topTags(t, n = Inf) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = featName)
    return(df)
  })
  # bind rows to obtain output df
  outDf <- dplyr::bind_rows(dfList, .id = compName)
  # rename columns and return
  colnames(outDf) <- c(compName, featName, fcName, exprName, "F", pName, pAdjName)
  return(outDf)

}

#' Automatic edgeR analysis from counts
#'
#' Perform the default edgeR analysis given a count matrix, a sample
#' information data frame and a grouping column.
#'
#' @param counts The count matrix to analyze.
#' @param sampInfo sampInfo The sample information data frame.
#' @param groupCol groupCol The column containing the grouping variable
#' @param removeColName Remove column name from design matrix?
#' @param ... Other arguments passed to \link[biokit]{edgerDfFromContrasts}.
#' @param filterByExpr Apply \link[edgeR]{filterByExpr}?
#'
#' @return A data frame with the results in a tidy format.
#'
#' @export
#'
autoEdgerComparison <- function(counts, sampInfo, groupCol, removeColName = TRUE,
                                filterByExpr = FALSE, ...) {

  # check sampInfo
  sampInfo <- validateSampInfo(sampInfo, groupCol, counts, checkNames = TRUE)
  # creat design and contrast matrix and perform analysis
  desMat <- designFromSampInfo(sampInfo = sampInfo, groupCol = groupCol, removeColName = removeColName)
  conMat <- contrastsFromDesign(desMat = desMat)
  results <- edgerDfFromContrasts(counts = counts, desMat = desMat,
                                  conMat = conMat, filterByExpr = filterByExpr, ...)
  return(results)

}
