#' Obtain limma results from data, design and contrasts matrix.
#'
#' Creates a linear model with limma and fits it to the provided contrasts. Then uses \link[limma]{eBayes}
#' to smooth standard errors and obtain the list of pairwise comparison results with \link[limma]{topTable}.
#'
#' @param mat Matrix with data to analyse.
#' @param desMat Design matrix.
#' @param conMat Contrasts matrix.
#' @param compName Name for the column that indicates the pairwise comparison in the tidy data frame.
#' @param featName Name for the column that indicates the analyzed row feature in the tidy data frame.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#' @param exprName Name for the resulting average expresion column.
#'
#' @return A data frame with the results in a tidy format.
#' @export
#'
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @importFrom dplyr bind_rows %>%
#' @importFrom tibble rownames_to_column
#'
limmaDfFromContrasts <- function(mat, desMat, conMat, compName = "comparison",
                                 featName = "feature", exprName = "AveExpr",
                                 fcName = "logFc", pName = "pValue", pAdjName = "pAdj") {

  # check input matrix
  inputCheck <- validateMatrix(mat = mat, checkRowNames = TRUE, allowNa = FALSE)
  # fit model
  fit <- limma::lmFit(object = mat, design = desMat)
  fit2 <- limma::contrasts.fit(fit, contrasts = conMat)
  # empirical bayes smoothing to the standard errors
  fit2 <- limma::eBayes(fit2)
  # obtain list of dataframes from contrasts column names
  dfList <- lapply(colnames(conMat), function(y) {
    resDf <- limma::topTable(fit2, coef = y, number = Inf)
    # if the ID column is present, then return directly
    if("ID" %in% colnames(resDf)) {
      return(resDf)
    } else {
      resDf <- tibble::rownames_to_column(resDf, var = featName)
      return(resDf)
    }
  })
  # set names to list
  names(dfList) <- colnames(conMat)
  # bind rows
  outDf <- dplyr::bind_rows(dfList, .id = compName)
  # set column names as indicated and return formatted df
  colnames(outDf) <- c(compName, featName, fcName, exprName, "t", pName, pAdjName, "B")
  return(outDf)

}

#' Automatic limma analysis
#'
#' Perform the default limma analysis given a data matrix, a sample
#' information data frame and a grouping column.
#'
#' @param mat Matrix with data to analyse.
#' @param sampInfo sampInfo The sample information data frame.
#' @param groupCol groupCol The column containing the grouping variable
#' @param removeColName Remove column name from design matrix?
#' @param ... Other arguments passed to \link[biokit]{limmaDfFromContrasts}.
#'
#' @return A data frame with the results in a tidy format.
#'
#' @export
#'
autoLimmaComparison <- function(mat, sampInfo, groupCol, removeColName = TRUE, ...) {

  # check sampInfo
  sampInfo <- validateSampInfo(sampInfo, groupCol, mat, checkNames = TRUE)
  # creat design and contrast matrix and perform analysis
  desMat <- designFromSampInfo(sampInfo = sampInfo, groupCol = groupCol, removeColName = removeColName)
  conMat <- contrastsFromDesign(desMat = desMat)
  results <- limmaDfFromContrasts(mat = mat, desMat = desMat, conMat = conMat, ...)
  return(results)

}
