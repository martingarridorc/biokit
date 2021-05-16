#' Default automatic biokit analysis
#'
#' This wrapper function performs the  "default" analyses of the package.
#'
#' @param mat The input matrix containing normalized data or counts
#' @param sampInfo The sample information data frame.
#' @param groupCol The column containing the grouping variable
#' @param funCatList A list of character vectors with the functional categories to analyze.
#' @param statMode Type of statistical analysis to perform the row-wise comparison.
#' Allowed values are: "edgeR", "limma" "tTest",
#' which will apply \link[biokit]{autoEdgerComparison}, \link[biokit]{autoLimmaComparison} and
#' \link[biokit]{autoPairwiseMatrixTest}, respectively.
#' @param filterByExpr If "edgeR" is selected, use \link[edgeR]{filterByExpr}?
#' @param annotationMode Approach to annotate differential features. Allowed values are: "byCutoff or "byRank",
#' which will apply \link[biokit]{addStatusByCutoff} or \link[biokit]{addStatusByRank}.
#' @param splitUpDown Split annotation into up and down features? If TRUE, splitCol and splitMidPoint should be provided.
#' @param splitCol If splitUpDown is set to TRUE, the column used to split features into up and down.
#' @param splitMidPoint If splitUpDown is set to TRUE, the splitCol value considered as midpoint.
#' @param sigCutoff If annotation mode is "byCutoff", the adjusted p value cutoff to use.
#' @param metricCutoff If annotation mode is "byCutoff", the absolute log FC cutoff to use.
#' @param rankCol If annotation mode is "byRank", the column used to rank results.
#' @param topN If annotation mode is "byRank", the top N features to annotate.
#' @param decreasing If annotation mode is "byRank", order rankCol decreasing?
#' @param absolute If annotation mode is "byRank", use absolute rankCol values?
#' @param minFunCatSize Minimal size of a functional category to test. All functional categories below the threshold are excluded.
#' @param maxFunCatSize Minimal size of a functional category to test. All functional categories aboce the threshold are excluded.
#' @param gseaRankCol Column used to create pre-ranked gsea list.
#' @param gseaNPerm Number of permutations to perform in the pre-ranked gsea analysis.
#'
#' @return A list containing three data frames: Statistical results, over-representation analysis results and gene set enrichment analysis results.
#'
#' @export
#'
autoBiokitAnalysis <- function(mat, sampInfo, groupCol, funCatList, statMode, annotationMode,
                               filterByExpr, splitUpDown = TRUE, splitCol = "logFc", splitMidPoint = 0,
                               sigCutoff = 0.05, metricCutoff = 1,
                               rankCol = "logFc", topN = 50, decreasing = TRUE, absolute = TRUE,
                               minFunCatSize = 1, maxFunCatSize = Inf,
                               gseaRankCol = "logFc", gseaNPerm = 10000) {

  # check allowed stat mode
  if( ! statMode %in% c("limma", "edgeR", "tTest")) stop("Supplied statMode does not exist.")
  if( ! annotationMode %in% c("byCutoff", "byRank", "tTest")) stop("Supplied annotationMode does not exist.")
  # perform stat analysis
  if(statMode == "limma") {
    statRes <- autoLimmaComparison(mat = mat, sampInfo = sampInfo, groupCol = groupCol, removeColName = TRUE)
  }
  if(statMode == "edgeR") {
    statRes <- autoEdgerComparison(counts = mat, sampInfo = sampInfo, groupCol = groupCol, removeColName = TRUE, filterByExpr = filterByExpr)
  }
  if(statMode == "tTest") {
    statRes <- autoPairwiseMatrixTest(mat = mat, sampInfo = sampInfo, groupCol = groupCol)
  }
  # annotate results with up and down features
  if(annotationMode == "byCutoff") {
    statRes <- splitFunMerge(statRes, splitCol = "comparison", fun = function(x) addStatusByCutoff(resDf = x, splitUpDown = splitUpDown,
                                                                                                   splitCol = splitCol, splitMidPoint = splitMidPoint,
                                                                                                   sigCutoff = sigCutoff, metricCutoff = metricCutoff))
  }
  if(annotationMode == "byRank") {
    statRes <- splitFunMerge(statRes, splitCol = "comparison", fun = function(x) addStatusByRank(resDf = x, splitUpDown = splitUpDown,
                                                                                                 splitCol = splitCol, splitMidPoint = splitMidPoint,
                                                                                                 topN = 50, rankCol = rankCol,
                                                                                                 decreasing = decreasing, absolute = absolute))
  }
  # perform over-representation analysis
  oraRes <- splitFunMerge(statRes, splitCol = "comparison", addSplitCol = TRUE, fun = function(x) oraFromStats(df = x, funCatList = funCatList,
                                                                                                               statusCol = "status",
                                                                                                               minSize = minFunCatSize,
                                                                                                               maxSize = maxFunCatSize))
  # perform gene set enrichment analysis
  gseaRes <- splitFunMerge(statRes, splitCol = "comparison", addSplitCol = TRUE, fun = function(x) gseaFromStats(df = x, funCatList = funCatList,
                                                                                             minSize = minFunCatSize, maxSize = maxFunCatSize,
                                                                                             rankCol = gseaRankCol , nperm = gseaNPerm))
  # prepare and return output
  outList <- list(statResults = statRes, oraResults = oraRes, gseaResults = gseaRes)
  return(outList)

}






