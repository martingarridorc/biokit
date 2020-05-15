#' Measure over-representation using Fisher's Exact Test
#'
#' Uses \link[stats]{fisher.test} to calculate the over-representation of
#' a set of features on a functional category given an "universe" that is used
#' as statistical background. Based on the Fisher's Exact Test over a contingency
#' table.
#'
#' @param features Input features.
#' @param funCat All the features in the functional category.
#' @param universe All the universe features.
#'
#' @return The over-representation significance as a P value.
#'
#' @export
#'
#' @importFrom stats fisher.test
#'
#' @examples
#'
#' g <- letters[1:5]
#' f <- letters[2:6]
#' f2 <- letters[c(1, 10:15)]
#' u <- letters
#'
#' # this should give a low p value (high over-representation)
#' fisherExactTest(g, f, u)
#' # this should give a high p value (low over-representation)
#' fisherExactTest(g, f2, u)
#'
fisherExactTest <- function(features, funCat, universe) {

  ### adapted from http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
  contingencyTable <- matrix(nrow = 2, ncol = 2)
  rownames(contingencyTable) <- c("inCategory", "outCategory")
  colnames(contingencyTable) <- c("altered", "notAltered")
  # prepare values
  x <- sum(features %in% funCat)
  k <- length(funCat)
  m <- length(features)
  N <- length(universe)
  n <- N - m
  # fill contingency table
  contingencyTable["inCategory", "altered"] <- x
  contingencyTable["inCategory", "notAltered"] <- k - x
  contingencyTable["outCategory", "altered"] <- m - x
  contingencyTable["outCategory", "notAltered"] <- n - (k - x)
  # apply fisher exact test
  p <- stats::fisher.test(x = contingencyTable, alternative = "greater")$p.value
  return(p)

}

#' Over-representation analysis
#'
#' Uses two different approaches to assess the significance of the over-representation
#' of a feature vector over a list of functional categories. Uses \link[biokit]{fisherExactTest}.
#'
#' @param features A vecor of characters that contains the features to analyze.
#' @param funCatList A list of character vectors with the functional categories to analyze.
#' @param universe A vector of characters to be used as the statistical background for the different tests. If not supplied, defaults to all the features in funCatList.
#' @param minSize Minimal size of a functional category to test. All functional categories below the threshold are excluded.
#' @param maxSize Minimal size of a functional category to test. All functional categories aboce the threshold are excluded.
#' @param pAdjustMethod The method used to correct estimated p values. Passed to \link[stats]{p.adjust}.
#'
#' @return A data frame with results in a tidy format
#'
#' @export
#'
#' @importFrom dplyr bind_rows %>%
#'
#' @examples
#'
#' data("humanHallmarks")
#' f <- sample(humanHallmarks$gene_symbol, 500)
#' fList <- split(humanHallmarks, humanHallmarks$gs_name)
#' fList <- lapply(fList, function(x) as.character(x$gene_symbol))
#' overRepresentationAnalysis(features = f, funCatList = fList)
#'
overRepresentationAnalysis <- function(features, funCatList, universe,
                                       minSize = 1, maxSize = Inf, pAdjustMethod = "BH") {

  # if missing universe, use all features in the functional category lis
  if(missing(universe)) {
    universe <- unique(unlist(funCatList))
  }
  # remove categories out of size
  size <- lapply(funCatList, length) %>%
    unlist()
  funCatList <- funCatList[size >= minSize & size <= maxSize]
  # calculate overlap and p value for each funCat, returning a one-row data frame with results
  dfList <- lapply(funCatList, function(funCat) {

    overlappingFeatures <- features[features %in% funCat]
    if(length(overlappingFeatures) == 0) {
      overlappingFeatures <- "No overlap"
    } else {
      overlappingFeatures <- paste(overlappingFeatures, collapse = ", ")
    }
    # get enrichment p values
    fisherP <- fisherExactTest(features = features, funCat = funCat, universe = universe)
    # return out df
    outDf <- data.frame(overlap = overlappingFeatures, pValue = fisherP)
    return(outDf)

  })
  # bind rows
  resDf <- dplyr::bind_rows(dfList, .id = "functionalCategory")
  # adjust p values
  resDf[ , "pAdj"] <- p.adjust(resDf$pValue, method = pAdjustMethod)
  # return data frame
  return(resDf)

}
