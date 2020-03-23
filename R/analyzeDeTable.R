#' Default analysis from comparison table
#'
#' Uses a column to identify if it is a single-comparison or a multi-comparison data frame and perform
#' subsequent analyses.
#'
#' @param deTable Data frame with stat results.
#' @param compCol Column containing the comparison id.
#' @param ... Rest of arguments passed to \link[biokit]{defaultSingleComparison} or \link[biokit]{defaultMultiComparison}
#'
#' @return
#'
#' @export
#'
analyzeDeTable <- function(deTable, compCol = "comparison", ...) {

  # get number of comparisons
  nComp <- length(unique(deTable[,compCol]))
  if(nComp == 1) {
    outList <- analyzeSingleComp(deTable = deTable, ...)
  } else {
    outList <- analyzeMultiComp(deTable = deTable, ...)
  }
  return(outList)

}

#' Default analysis from multi comparison data frame.
#'
#' Performs the "default" analysis given a tidy data frame with the statistical results of several comparisons.
#'
#' @param deTable Data frame with stat results.
#' @param funCategories Functional categories to use in \link[clusterProfiler]{enricher} and \link[clusterProfiler]{GSEA}.
#' @param useCutoff Use cutoff to annotate up and down features?
#' @param fcCutoff Fc cutoff. Passed to \link[biokit]{annotateByCutoff}.
#' @param pValueCutoff P value cutoff. Passed to \link[biokit]{annotateByCutoff}.
#' @param n N features to annotate. Passed to \link[biokit]{annotateTopN}.
#' @param sortCol Column used to sort data.  Passed to \link[biokit]{annotateTopN}.
#'
#' @return Data frame containing results of annotating deTable, ORA and GSEA and several default plots.
#'
#' @importFrom dplyr %>%
#'
#' @export
#'
analyzeMultiComp <- function(deTable, funCategories, useCutoff = TRUE, fcCutoff = 1,
                             pValueCutoff = 0.05, n = 50,  sortCol = "pAdj") {

  # annotate
  if(useCutoff) deTable <- annotateMultiComparison(x = deTable, useCutoff = useCutoff, metricCutoff = fcCutoff, sigCutoff = pValueCutoff)
  if(!useCutoff) deTable <- annotateMultiComparison(x = deTable, useCutoff = useCutoff, n = n, sortCol = sortCol)
  # ora
  oraTable <- splitFeatures(deTable, splitCol = "comparison") %>%
    oraFromList(x = ., funCategories = funCategories, pvalueCutoff = 1, qvalueCutoff = 1) %>%
    cpResultsToDf()
  # gsea
  gseaTable <- getRankedVectorList(deTable, splitCol = "comparison") %>%
    gseaFromList(x = . , funCategories = funCategories, pvalueCutoff = 1) %>%
    cpResultsToDf()
  # create output list with all tables and plots
  outList <- list(deTable = deTable, oraTable = oraTable, gseaTable = gseaTable)
  return(outList)

}

#' Default analysis from single comparison data frame.
#'
#' Performs the "default" analysis given a tidy data frame with the statistical results.
#'
#' @param deTable Data frame with stat results.
#' @param funCategories Functional categories to use in \link[clusterProfiler]{enricher} and \link[clusterProfiler]{GSEA}.
#' @param useCutoff Use cutoff to annotate up and down features?
#' @param fcCutoff Fc cutoff. Passed to \link[biokit]{annotateByCutoff}.
#' @param pValueCutoff P value cutoff. Passed to \link[biokit]{annotateByCutoff}.
#' @param n N features to annotate. Passed to \link[biokit]{annotateTopN}.
#' @param sortCol Column used to sort data.  Passed to \link[biokit]{annotateTopN}.
#'
#' @return Data frame containing results of annotating deTable, ORA and GSEA and several default plots.
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom clusterProfiler GSEA
#'
analyzeSingleComp <- function(deTable, funCategories, useCutoff = TRUE, fcCutoff = 1,
                              pValueCutoff = 0.05, n = 50,  sortCol = "pAdj") {

  # annotate
  if(useCutoff) deTable <- annotateByCutoff(deTable, metricCutoff = fcCutoff, sigCutoff = pValueCutoff)
  if(!useCutoff) deTable <- annotateTopN(deTable, n = n, sortCol = sortCol)
  # ora
  oraTable <- splitFeatures(deTable) %>%
    oraFromList(x = ., funCategories = funCategories, pvalueCutoff = 1, qvalueCutoff = 1) %>%
    cpResultsToDf()
  # gsea
  gseaTable <- getRankedVector(deTable) %>%
    # use cp functions
    clusterProfiler::GSEA(geneList = . , TERM2GENE =  funCategories, pvalueCutoff = 1) %>%
    .@result
  # create output list with all tables and plots
  outList <- list(deTable = deTable, oraTable = oraTable, gseaTable = gseaTable)
  return(outList)

}









