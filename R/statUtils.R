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
#' @param transpose Wether to transpose or not the input matrix.
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
    pcPcts <- paste0(colnames(pcaRes$x), " ( ", variancePcts, "% )")
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

#' Add status column to comparison data frame using cutoff
#'
#' Uses a fold change and P value cutoff to annotate a
#' data frame containing the results of a stat comparison.
#'
#' @param x Data frame to annotate.
#' @param fcCutoff Fold change cutoff. Should be positive, as it is transformed into negative when annotating 'downregulated' features.
#' @param pCutoff P Value cutoff.
#' @param splitUpDown Wethter to split generated label into up and down features.
#' @param fcCol Column used to filter by fc.
#' @param pCol Column used to filter by p value.
#' @param noChangeLabel Label used for not-changing features.
#' @param statusCol Name of the new column containing labels.
#' @param upLabel Label for up-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param downLabel Label for down-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param sigLabel Label for significantly altered features. Only used if \code{splitUpDown = FALSE}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
annotateByCutoff <- function(x, fcCutoff = 1, pCutoff = 0.05, splitUpDown = TRUE, fcCol = "logFc", pCol = "pAdj", noChangeLabel = "No change", statusCol = "status",
    upLabel = "Up", downLabel = "Down", sigLabel = "Significant") {

    x[, statusCol] <- noChangeLabel
    # if split, then return up and down labels, else return significant or not significant
    if (splitUpDown) {
        up <- x[, fcCol] >= fcCutoff & x[, pCol] <= pCutoff
        down <- x[, fcCol] <= -fcCutoff & x[, pCol] <= pCutoff
        x[up, statusCol] <- upLabel
        x[down, statusCol] <- downLabel
    } else {
        sig <- abs(x[, fcCol]) >= fcCutoff & x[, pCol] <= pCutoff
        x[sig, statusCol] <- sigLabel
    }
    return(x)

}

#' Add status column to comparison data frame using top N features
#'
#' Orders a data frame and adds a column to indicate the
#' status of the top N features.
#'
#' @param x Data frame to annotate.
#' @param n Top N features to be annotated
#' @param sortCol Column used to sort the data frame.
#' @param splitUpDown Wether to split top N features into up and down regulated using a fold change column. This will annotate 2*n features.
#' @param decreasing Order sortCol decreasing.
#' @param twoSides Pick up both extremes of the sorted data frame as the top N features. SplitUpDown should be TRUE.
#' @param fcCol Fold change column used to sepparate up and down features.
#' @param noChangeLabel Label used for not-changing features.
#' @param statusCol Name of the new column containing labels.
#' @param upLabel Label for up-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param downLabel Label for down-regulated features. Only used if \code{splitUpDown = TRUE}.
#' @param sigLabel Label for significantly altered features. Only used if \code{splitUpDown = FALSE}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
annotateTopN <- function(x, n, sortCol, decreasing = FALSE, splitUpDown = TRUE, twoSides = FALSE, fcCol = "logFc", noChangeLabel = "No change", statusCol = "status",
    upLabel = "Up", downLabel = "Down", sigLabel = "Significant") {

    x[, statusCol] <- noChangeLabel
    x <- x[order(x[, sortCol], decreasing = decreasing), ]
    # if split and twoSides, annotate both extremes of the df
    if (splitUpDown & twoSides) {
        x[1:n, statusCol] <- upLabel
        x[(nrow(x) - n + 1):nrow(x), statusCol] <- downLabel
    }
    # if split and !twoSides, pick top N
    if (splitUpDown & !twoSides) {
        x[x[, fcCol] > 0, statusCol][1:n] <- upLabel
        x[x[, fcCol] < 0, statusCol][1:n] <- downLabel
    }
    # else return top N
    if (!splitUpDown & !twoSides) {
        x[, statusCol][1:n] <- sigLabel
    }
    return(x)

}

#' Annotate multiple comparisons data frame
#'
#' Uses \link[biokit]{annotateByCutoff} or \link[biokit]{annotateTopN} to annotate a
#' data frame by sepparating data for each comparison.
#'
#' @param x Multi-comparison data frame to annotate.
#' @param useCutoff Wether to use \link[biokit]{annotateByCutoff} or \link[biokit]{annotateTopN}.
#' @param compCol Column with the comparison id.
#' @param ... Rest of arguments passed to \link[biokit]{annotateByCutoff} or \link[biokit]{annotateTopN}.
#'
#' @return Annotated data frame.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#'
annotateMultiComparison <- function(x, useCutoff, compCol = "comparison", ...) {

    # split list into comparisons
    dfList <- split(x, f = x[, compCol])
    # apply function of interest
    if (useCutoff) {
        aDfList <- lapply(dfList, function(y) annotateByCutoff(x = y, ...))
    } else {
        aDfList <- lapply(dfList, function(y) annotateTopN(x = y, ...))
    }
    # coerce resultant data frames
    result <- dplyr::bind_rows(aDfList, .id = compCol)
    return(result)

}
