#' Create contrast matrix given a design matrix.
#'
#' Uses levels defined in the column names of the design matrix
#' to create a contrast matrix with \link[limma]{makeContrasts}.
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
#' Creates a linear model with limma and fits it to the provided contrasts. Then uses \link[limma]{eBayes}
#' to smooth standard errors and obtain the list of pairwise comparison results with \link[limma]{topTable}.
#'
#' @param x Matrix with data to analyse.
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
#' @importFrom dplyr bind_rows
#' @importFrom tibble rownames_to_column
#'
limmaDfFromContrasts <- function(x, desMat, conMat, compName = "comparison", featName = "feature", exprName = "AveExpr", 
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
    colnames(outDf) <- c(compName, featName, fcName, exprName, "t", pName, pAdjName, "B")
    return(outDf)
    
}

#' Automatic limma analysis from SummarizedExperiment
#'
#' Uses data stored in a SummarizedExperiment object to
#' automatically perform all the possible pairwise comparisons
#' between the groups defined in a \link[SummarizedExperiment]{colData} column.
#'
#' @param se The SummarizedExperiment object to analyze.
#' @param groupCol The id of the \link[SummarizedExperiment]{colData} column that is used to create the pairwise contrasts.
#' @param ... Arguments passed to the \link[biokit]{limmaDfFromContrasts} function.
#'
#' @return A data frame with the results in a tidy format.
#' @export
#'
#' @import SummarizedExperiment
#'
autoLimma <- function(se, groupCol, ...) {
    
    desMat <- designFromSampInfo(x = colData(se), column = groupCol)
    conMat <- contrastsFromDesign(x = desMat)
    results <- limmaDfFromContrasts(x = assay(se), desMat = desMat, conMat = conMat, ...)
    return(results)
    
}
