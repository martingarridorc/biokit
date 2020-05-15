#' Not sensitive T test p value
#'
#' @param ... Arguments passed to \link[stats]{t.test}.
#'
#' @return P value resulting from the T Test or \code{NA} when any problem appears.
#'
#' @importFrom stats t.test
#'
nsTestPValue <- function(...) {

  pValue <- tryCatch(stats::t.test(...)$p.value, error = function(y) NA)
  return(pValue)

}

#' Non sensitive T Test statistic
#'
#' @param ... Arguments passed to \link[stats]{t.test}.
#'
#' @return Statistic resulting from the T Test or \code{NA} when any problem appears.
#'
#' @importFrom stats t.test
#'
nsTestT <- function(...) {

  t <- tryCatch(stats::t.test(...)$statistic, error = function(y) NA)
  return(t)

}

#' Log 2 of means ratio
#'
#' @param x Numerator vector.
#' @param y Denominator vector.
#' @param na.rm Remove NAs?
#'
#' @return Log 2  of means ratio.
#'
defaultLogFc <- function(x, y, na.rm = TRUE) {

  metric <- log2(mean(x, na.rm = na.rm) / mean(y, na.rm = na.rm))
  return(metric)

}

#' Split, Function and Merge
#'
#' Splits the input data frame using a grouping column, applies
#' the provided function over each list element and return the row-binded
#' resulting data frame.
#'
#' @param df Data frame to process.
#' @param splitCol Column to split by.
#' @param fun Function to be applied to each sub-data frame.
#'
#' @return The processed data frame.
#'
#' @importFrom dplyr bind_rows
#'
splitFunMerge <- function(df, splitCol, fun) {

  splitted <- split(x = df, f = df[, splitCol])
  splittedRes <- lapply(splitted, fun)
  binded <- dplyr::bind_rows(splittedRes, .id = splitCol)
  return(binded)

}
