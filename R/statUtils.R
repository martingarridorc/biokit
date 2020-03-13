#' Create pariwise contrasts
#'
#' Creates all the possible pairwise contrasts between
#' the unique elements of a vector.
#'
#' @param x Vector to be combined in pairs.
#' @param sep Separator for the resulting contrasts.
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
#' @param x Matrix to analyze.
#' @param transpose Transpose input matrix?
#' @param roundDigits Number of digits for the rounded percent of
#' variance explained by each PC.
#' @param ... Other parameters passed to the \link[stats]{prcomp} function.
#'
#' @return A list with three elements:
#' \item{result}{Prcomp object containing the results of the PCA.}
#' \item{summary}{Summary object containing the matrix with PC values.}
#' \item{pcts}{Character vector with the variance explained by each PC.}
#'
#' @export
#'
#' @importFrom stats prcomp
#'
pcaToList <- function(x, transpose = TRUE, roundDigits = 2, ...) {
    
    # transpose if it is not transposed
    if (transpose) 
        x <- t(x)
    # perform prcomp analysis
    pcaRes <- prcomp(x, ...)
    pcaResSum <- summary(pcaRes)
    # add proportion of variance explained
    variancePcts <- round(pcaResSum$importance[2, ] * 100, roundDigits)
    pcPcts <- paste0(colnames(pcaRes$x), " (", variancePcts, "%)")
    return(list(result = pcaRes, summary = pcaResSum, pcts = pcPcts))
    
}

#' Not sensitive T Test
#'
#' Applies the T Test function returning NA
#' instead of error when problems appear. Additionally,
#' it extracts the P value from the resulting object when available.
#'
#' @param ... Arguments for the \link[stats]{t.test} function.
#'
#' @return P value resulting from the Test or \code{NA} when any problem appears.
#'
#' @export
#'
#' @importFrom stats t.test
#'
nsTest <- function(...) {
    
    pValue <- tryCatch(t.test(...)$p.value, error = function(y) NA)
    return(pValue)
    
}

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
