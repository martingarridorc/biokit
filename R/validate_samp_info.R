#' Validate sample information
#'
#' Evaluates integrity of the grouping variable at the sample information data frame.
#'
#' @param sampInfo The sample information data frame.
#' @param groupCol The column containing the grouping variable
#' @param checkNames Check that grouping variable contains syntactically valid names? If TRUE, they are fixed in the output data frame.
#' @param checkSingleSample Check that any group contains only one sample? A warning message appears if true.
#'
#' @return TRUE if no errors are found.
#'
validateSampInfo <- function(sampInfo, groupCol, checkNames = TRUE, checkSingleSample = TRUE) {

  # extract grouping variable
  groupVar <- sampInfo[ , groupCol]
  # check not syntactically valid names
  if(checkNames) {
    goodNames <- make.names(groupVar)
    if(any(goodNames != groupVar)) {
      message("Some levels at the groupìng variable are not syntactically valid names. They will be fixed.")
      sampInfo[ , groupCol] <- goodNames
    }
  }
  # check presence of single sample groups
  if(checkSingleSample) {
    groupN <- table(sampInfo[, groupCol])
    if(any(groupN == 1)) {
      warning("Some sample groups contains only 1 sample. This may produce unexpected results.")
    }
  }
  return(sampInfo)

}
