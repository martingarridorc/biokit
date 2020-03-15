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



