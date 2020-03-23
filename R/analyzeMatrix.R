#' Default analysis of counts in SummarizedExperiment
#'
#' @param se Summarized Experiment.
#' @param groupCol Column from \link[SummarizedExperiment]{colData} used to group samples.
#' @param funCategories Functional categories to use in \link[clusterProfiler]{enricher} and \link[clusterProfiler]{GSEA}.
#' @param ... Rest of arguments passed to \link[biokit]{defaultFromTable}.
#'
#' @return A list with results
#'
#' @export
#'
analyzeCounts <- function(se, groupCol, funCategories, ...) {

  # obtain TMMs
  normExpr <- countsToTmm(assay(se))
  # obtain deTable
  deTable <- autoEdgeR(se, groupCol)
  # obtain deTable
  result <- analyzeDeTable(deTable, compCol = "comparison", funCategories = mouseHallmarks, ...)
  # get pcaPlot
  pcaRes <- pcaToList(normExpr)
  # add elements to result list
  result[["normExpr"]] <- normExpr
  return(result)

}

#' Default analysis of normalized matrix in SummarizedExperiment
#'
#' @param se Summarized Experiment.
#' @param groupCol Column from \link[SummarizedExperiment]{colData} used to group samples.
#' @param funCategories Functional categories to use in \link[clusterProfiler]{enricher} and \link[clusterProfiler]{GSEA}.
#' @param useLimma use \link[biokit]{autoLimma}? If not, uses \link[package]{autoMatrixPairwise}.
#' @param ... Rest of arguments passed to \link[biokit]{defaultFromTable}.
#'
#' @return A list with results
#'
#' @export
#'
analyzeMatrix <- function(se, groupCol, funCategories, useLimma = TRUE, ...) {

  # use autoLimma or autoMatrixPairwise to obtain deTable
  if(useLimma) deTable <- autoLimma(se, groupCol)
  if(!useLimma) deTable <- autoMatrixPairwise(se, groupCol)
  # obtain deTable
  result <- analyzeDeTable(deTable, compCol = "comparison", funCategories = funCategories, ...)
  return(result)

}




