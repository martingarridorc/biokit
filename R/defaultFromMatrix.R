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
defaultFromCounts <- function(se, groupCol, funCategories, ...) {

  # obtain TMMs
  normExpr <- countsToTmm(assay(se))
  # obtain deTable
  deTable <- autoEdgeR(se, groupCol)
  # obtain deTable
  result <- defaultFromTable(deTable, compCol = "comparison", funCategories = mouseHallmarks, gseaPCutoff = 1, oraPCutoff = 1, ...)
  # get pcaPlot
  pcaRes <- pcaToList(normExpr)
  # plot results
  pcaPlot <- defaultPcaPlot(x = pcaRes, sampInfoDf = as.data.frame(colData(se)), groupCol = "group")
  # add elements to result list
  result[["normExpr"]] <- normExpr
  result[["pcaPlot"]] <- pcaPlot
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
defaultFromNormMatrix <- function(se, groupCol, funCategories, useLimma = TRUE, ...) {

  # use autoLimma or autoMatrixPairwise to obtain deTable
  if(useLimma) deTable <- autoLimma(se, groupCol)
  if(!useLimma) deTable <- autoMatrixPairwise(se, groupCol)
  # obtain deTable
  result <- defaultFromTable(deTable, compCol = "comparison", funCategories = funCategories, gseaPCutoff = 1, oraPCutoff = 1, ...)
  # get pcaPlot
  pcaRes <- pcaToList(assay(se))
  # plot results
  pcaPlot <- defaultPcaPlot(x = pcaRes, sampInfoDf = as.data.frame(colData(se)), groupCol = "group")
  # add elements to result list
  result[["pcaPlot"]] <- pcaPlot
  return(result)

}




