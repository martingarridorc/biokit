#' Over representation analysis from statistical results
#'
#' Performs an over-representation analysis from a data frame containing the results of
#' the statistical analysis of omic data. Uses the annotation column to split
#' the features and perform the over-representation analysis.
#'
#' @param df Data frame with results to analyze.
#' @param funCatList A list of character vectors with the functional categories to analyze.
#' @param statusCol Column containing the annotation of features.
#' @param noChangeLabel Label used to indicate not changing features. Rows with a statusCol value equal to this will be removed.
#' @param featCol Column containing the feature to submit for the over-representation analysis (i.e gene symbol).
#' @param ... Other parameters passed to \link[biokit]{overRepresentationAnalysis}.
#'
#' @return A tidy data frame with the results.
#'
#' @export
#'
oraFromStats <- function(df, funCatList, statusCol, noChangeLabel = "No change",
                         featCol = "feature", ...) {

  # remove rows annotated as "not changing"
  df <- subset(df, df[ , statusCol] != noChangeLabel)
  # split data frame, apply ora and return binded df
  resDf <- splitFunMerge(df = df, splitCol = statusCol, fun = function(x) {
    toOra <- as.character(x[ , featCol])
    oraRes <- overRepresentationAnalysis(features = toOra, funCatList = funCatList, ...)
    return(oraRes)
  })
  return(resDf)

}

#' Gene set enrichment analysis from statistical results
#'
#' @param df Data frame with results to analyze.
#' @param funCatList A list of character vectors with the functional categories to analyze.
#' @param rankCol Column used to perform the pre-ranked gsea analysis.
#' @param featCol Column used to name the rankCol values.
#' @param seed Seed used to avoid reproducibility problems between runs.
#' @param ... Rest of arguments passed to \link[fgsea]{fgseaSimple}.
#'
#' @return A tidy data frame with the results.
#'
#' @export
#'
#' @importFrom fgsea fgseaSimple
#'
gseaFromStats <- function(df, funCatList, rankCol, featCol = "feature", seed =149, ...) {

  toGsea <- df[ , rankCol]
  names(toGsea) <- df[ , featCol]
  toGsea <- toGsea[order(toGsea, decreasing = TRUE)]
  # set seed to gain reproducibility between analyses
  set.seed(seed)
  gseaRes <- fgsea::fgseaSimple(pathways = funCatList, stats = toGsea, ...)
  return(gseaRes)

}
