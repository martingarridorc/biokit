#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep The separator for the resulting contrasts.
#'
#' @return A character vector containing all the possible pairwise contrasts.
#'
#' @export
#'
#' @importFrom utils combn
#'
pairwiseContrasts <- function(x, sep = "-") {

  # create all possible pairwise comparison matrix
  pairMat <- combn(unique(x), 2)
  # return collapsed columns
  contrasts <- apply(pairMat, 2, function(y) paste(y[1], y[2], sep = sep))
  return(contrasts)

}

#' Perform PCA and create list with results
#'
#' Performs principal component analysis for a given matrix
#' and returns a list with the formatted results.
#'
#'
#' @param x The matrix to analyze.
#' @param transpose Wether to transpose or not the input matrix.
#' @param roundDigits Number of digits for the rounded percent of
#' variance explained by each PC.
#' @param ... Other parameters passed to the prcomp function
#'
#' @return A list with three elements:
#' \item{result}{The prcomp object containing the results of the PCA.}
#' \item{summary}{The summary objects containing the matrix with PC values.}
#' \item{pcts}{A character vector with the variance explained by each PC.}
#'
#' @export
#'
#' @importFrom stats prcomp
#'
pcaToList <- function(x, transpose = TRUE, roundDigits = 2, ...) {

  # transpose if it is not transposed
  if(transpose) x <- t(x)
  # perform prcomp analysis
  pcaRes <- prcomp(x , ...)
  pcaResSum <- summary(pcaRes)
  # add proportion of variance explained
  variancePcts <- round(pcaResSum$importance[2, ]*100, roundDigits)
  pcPcts <- paste0(colnames(pcaRes$x), " ( ", variancePcts, "% )")
  return(list(result = pcaRes, summary = pcaResSum, pcts = pcPcts))

}

#' Not sensitive T Test
#'
#' Applies the T Test function returning NA
#' instead of error when problems appear (as continuous). Additionally,
#' it extracts the P value from the resulting object when available.
#'
#' @param ... Arguments for \code{t.test} function.
#'
#' @return The P value resulting from the T Test or \code{NA} when any problem appears.
#' @export
#'
#' @importFrom stats t.test
#'
nsTest <- function(...) {

  pValue <- tryCatch(t.test(...)$p.value, error= function(y) NA)
  return(pValue)

}

#' One sample T-Test over log ratio matrix
#'
#' Applies One Sample T-Test at row level in a matrix containing log transformed fold changes.
#'
#' @param x The logFC matrix.
#' @param adjustMethod The P value adjustment method.
#' @param idName Name for the resulting feature id column.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#' @param ... Rest of arguments for \code{t.test} function.
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
  if(is.null(rownames(x))) outId <- 1:nrow(x)
  if(!is.null(rownames(x))) outId <- rownames(x)
  # prepare out df
  outDf <- data.frame(outId, logFc, pValue, pAdj, stringsAsFactors = FALSE)
  colnames(outDf) <- c(idName, fcName, pName, pAdjName)
  return(outDf)

}

#' Create design matrix from sample metadata
#'
#' Uses a data frame containing sample metadata to create a design matrix
#' for subsequent analyses.
#'
#' @param x The data frame containing the column metadata (e.g sample grouping).
#' @param column The data frame column that is used to compose the design matrix.
#'
#' @return The composed design matrix.
#'
#' @export
#'
#' @importFrom stats formula model.matrix
#'
designFromSampInfo <- function(x, column) {

  # prepare formula
  groupFormula <- formula(paste0("~ 0 + ", column))
  # create design matrix
  designMatrix <- model.matrix(groupFormula, x)
  # remove column name in design colnames
  colnames(designMatrix) <- gsub(pattern = column, replacement = "", x = colnames(designMatrix))
  return(designMatrix)

}

#' Create contrast matrix given a design matrix.
#'
#' Uses levels defined in the column names of the design matrix
#' to create a contrast matrix with all the possible pairwise comparisons.
#'
#' @param x The design matrix.
#'
#' @return A contrast matrix with all the possible pairwise comparisons.
#' @export
#'
#' @importFrom limma makeContrasts
#'
contrastsFromDesign <- function(x) {

  # get contrasts from interesting column and build contrast matrix with limma
  contrasts <- pairwiseContrasts(colnames(x))
  contrastMatrix <- limma::makeContrasts(contrasts = contrasts, levels = x)
  return(contrastMatrix)

}

#' Obtain limma results from data, design and contrasts matrix.
#'
#' Creates a linear model with limma and fits it to the provided contrasts. Then uses \code{eBayes()}
#' to smooth standard errors and obtain the list of pairwise comparison results with \code{topTable()}.
#'
#' @param x Matrix with data to analyse.
#' @param desMat Design matrix.
#' @param conMat Contrasts matrix.
#' @param compName Name for the column that indicates the pairwise comparison in the tidy data frame.
#' @param featName Name for the column that indicates the analyzed row feature the tidy data frame.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#'
#'
#' @return A data frame with the results in a tidy format.
#' @export
#'
#' @importFrom limma lmFit contrasts.fit eBayes topTable
#' @importFrom dplyr bind_rows
#' @importFrom tibble rownames_to_column
#'
limmaDfFromContrasts <- function(x, desMat, conMat,
                                 compName = "comparison", featName = "feature",
                                 fcName = "logFc", pName = "pValue", pAdjName = "pAdj") {

  # fit model
  fit <- limma::lmFit(object = x, design = desMat)
  fit2 <- limma::contrasts.fit(fit, contrasts = conMat)
  # empirical bayes smoothing to the standard errors
  fit2 <- limma::eBayes(fit2)
  # obtain list of dataframes from contrasts column names
  dfList <- lapply(colnames(conMat), function(y) limma::topTable(fit2, coef = y, number = Inf))
  # rownames to feature column
  dfList <- lapply(dfList, function(y) tibble::rownames_to_column(y, var = featName))
  # set names to list
  names(dfList) <- colnames(conMat)
  # bind rows
  outDf <- dplyr::bind_rows(dfList, .id = compName)
  # set column names as indicated and return formatted df
  colnames(outDf) <- c(compName, featName, fcName, "AveExpr", "t", pName, pAdjName, "B")
  return(outDf)

}

#' Automatic limma analysis from SummarizedExperiment
#'
#' Uses data stored in a SummarizedExperiment object to
#' automatically perform all the possible pairwise comparisons
#' between the groups defined in a \code{colData()} column.
#'
#' @param se The SummarizedExperiment object to analyze.
#' @param groupColumn The id of the \code{colData()} column that is used to
#' @param compName Name for the column that indicates the pairwise comparison in the tidy data frame.
#' @param featName Name for the column that indicates the analyzed row feature the tidy data frame.
#' @param fcName Name for the resulting fold change column.
#' @param pName Name for the resulting p value column.
#' @param pAdjName Name for the resulting adjusted p value column.
#'
#' @return A data frame with the results in a tidy format.
#' @export
#'
#' @import SummarizedExperiment
#'
autoLimma <- function(se, groupColumn, compName = "comparison", featName = "feature",
                      fcName = "logFc", pName = "pValue", pAdjName = "pAdj") {

  designMatrix <- designFromSampInfo(x = colData(se), column = groupColumn)
  contrastMatrix <- contrastsFromDesign(x = designMatrix)
  results <- limmaDfFromContrasts(x = assay(se), desMat = designMatrix, conMat = contrastMatrix,
                                  compName = compName, featName =  featName,
                                  fcName = fcName, pName = pName, pAdjName = pAdjName)
  return(results)

}

#' Normalize counts with TMM
#'
#' Normalizes a raw counts matrix into log transformed counts per millions using
#' the edgeR "Trimmed Mean of M-Values" (TMM) normalization method. For more information
#' please refer to the edgeR package (https://bioconductor.org/packages/release/bioc/html/edgeR.html).
#'
#' @param x The matrix of counts.
#' @param log Passed to \code{edgeR::cpm()}.
#' @param prior.count Passed to \code{edgeR::cpm()}.
#'
#' @return The normalized expression matrix
#'
#' @export
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'
countsToTmm <- function(x, log = TRUE, prior.count = 3) {

  # transform into edgeR DGEList
  dgeList <- edgeR::DGEList(counts = x)
  # calculate normalization factors
  dgeList <- edgeR::calcNormFactors(object = dgeList, method = "TMM")
  # obtain normalized logCpm matrix
  logCpm <- edgeR::cpm(dgeList, log = log, prior.count = prior.count)
  return(logCpm)

}


autoEdgeR <- function(se, groupColumn, useFilterByExpr = TRUE,
                      compName = "comparison", featName = "feature",
                      fcName = "logFc", pName = "pValue", pAdjName = "pAdj") {






}




