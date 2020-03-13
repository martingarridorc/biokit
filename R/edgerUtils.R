#' Normalize counts with TMM
#'
#' Normalizes a raw counts matrix into log transformed counts per millions using
#' the edgeR 'Trimmed Mean of M-Values' (TMM) normalization method. For more information
#' please refer to the edgeR package (https://bioconductor.org/packages/release/bioc/html/edgeR.html).
#'
#' @param x The matrix of counts.
#' @param log Passed to \link[edgeR]{cpm}.
#' @param prior.count Passed to \link[edgeR]{cpm}.
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

#' Obtain edgeR results from object, design and contrasts matrix
#'
#' Applies recommended edgeR steps to obtain differential expression tables from
#' an initial \link[edgeR]{DGEList} object. As in \link[biokit]{limmaDfFromContrasts}, it uses
#' a design matrix and a contrast matrix as input.
#'
#' @param object The edgeR object to analyze. Generated with: \link[edgeR]{DGEList}.
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
#'
#' @export
#'
#' @import edgeR
#' @importFrom dplyr %>% bind_rows
#' @importFrom tibble rownames_to_column
#'
edgeRFromContrasts <- function(object, desMat, conMat, compName = "comparison", featName = "feature", exprName = "AveExpr", 
    fcName = "logFc", pName = "pValue", pAdjName = "pAdj") {
    
    # calculate normalization factors and estimate dispersion
    y <- edgeR::calcNormFactors(object)
    y <- edgeR::estimateDisp(y, desMat)
    # fit ql netgative binomial gl model
    fit <- edgeR::glmQLFit(y, desMat)
    # obtain DE table for each contrast
    dfList <- apply(conMat, 2, function(c) {
        t <- edgeR::glmQLFTest(fit, contrast = c)
        df <- edgeR::topTags(t, n = Inf) %>% as.data.frame() %>% tibble::rownames_to_column(var = featName)
        return(df)
    })
    # bind rows to obtain output df
    outDf <- dplyr::bind_rows(dfList, .id = compName)
    # rename columns and return
    colnames(outDf) <- c(compName, featName, fcName, exprName, "F", pName, pAdjName)
    return(outDf)
    
}

#' Automatic edgeR analysis from SummarizedExperiment with counts
#'
#' Uses a count matrix stored in a SummarizedExperiment object to
#' automatically perform all the possible pairwise comparisons
#' between the groups defined in a \link[SummarizedExperiment]{colData} column.
#'
#' @param se The SummarizedExperiment object to analyze.
#' @param groupColumn The id of the \link[SummarizedExperiment]{colData} column that is used to create the pairwise contrasts.
#' @param minCounts The minimum number of counts per feature that will be used to filter the count matrix.
#' @param useFilterByExpr Apply \link[edgeR]{filterByExpr}?
#' @param ... Arguments passed to the \link[biokit]{edgeRFromContrasts} function.
#'
#' @return A data frame with the results in a tidy format.
#' @export
#'
#' @importFrom edgeR DGEList filterByExpr
#'
autoEdgeR <- function(se, groupColumn, minCounts = 0, useFilterByExpr = TRUE, ...) {
    
    # obtain design matrix
    desMat <- designFromSampInfo(x = colData(se), column = groupColumn)
    # obtain contrasts
    conMat <- contrastsFromDesign(x = desMat)
    # get counts matrix and filter by minCounts
    counts <- assay(se)
    if (minCounts != 0) 
        counts <- counts[rowSums(counts) >= minCounts, ]
    # prepare DGEList object
    y <- edgeR::DGEList(counts = counts, group = colData(se)[, groupColumn])
    # filter by expr
    if (useFilterByExpr) {
        keep <- edgeR::filterByExpr(y)
        y <- y[keep, , keep.lib.sizes = FALSE]
    }
    # get results df
    results <- edgeRFromContrasts(object = y, desMat = desMat, conMat = conMat, ...)
    return(results)
    
}
