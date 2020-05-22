#' Validate matrix
#'
#' Evaluates integrity of the input matrix.
#'
#' @param mat The matrix to validate.
#' @param checkRowNames Check rownames are present? Stop execution if TRUE and mat does not contain rownames.
#' @param allowNa Allow NA presence in mat? Stop execution if NAs are found in mat.
#'
#'
#' @return TRUE if no errors are found.
#'
validateMatrix <- function(mat, checkRowNames = TRUE, allowNa = TRUE) {

  # check class
  matClass <- class(mat)
  if(! all(matClass == c("matrix", "array"))) stop("Input is not a matrix.")
  # check rownames
  if(checkRowNames & is.null(rownames(mat))) stop("Matrix does not have rownames.")
  # check NAs
  naCount <- sum(is.na(mat))
  if(naCount != 0 & allowNa) message("Matrix contain NAs. This may produce unexpected results.")
  if(naCount != 0 & (! allowNa)) stop("Matrix contain NAs.")
  return(TRUE)

}

#' Validate sample information
#'
#' Evaluates integrity of the grouping variable at the sample information data frame.
#'
#' @param sampInfo The sample information data frame.
#' @param groupCol The column containing the grouping variable.
#' @param checkNames Check if grouping variable contains syntactically valid names? If TRUE, they are fixed in the output data frame.
#' @param checkSingleSample Check if any group contains only one sample? A warning message appears if true.
#' @param mat Matrix containing data. If provided, it will check number of columns in mat.
#'
#' @return The checked sampInfo data frame
#'
validateSampInfo <- function(sampInfo, groupCol, mat,
                             checkNames = TRUE, checkSingleSample = TRUE) {

  # check not syntactically valid names
  if(checkNames) {
    groupVar <- sampInfo[ , groupCol]
    goodNames <- make.names(groupVar)
    if(any(goodNames != groupVar)) {
      message("Some levels at the grouping variable are not syntactically valid names. They will be fixed.")
      sampInfo[ , groupCol] <- goodNames
    }
  }
  # check number of rows is same than matrix ncol
  if(!missing(mat)) {
    if(ncol(mat) != nrow(sampInfo)) {
      stop("Number of rows for sampInfo is not equal to the number of matrix columns.")
    }
  }
  # check presence of single sample groups
  if(checkSingleSample) {
    groupN <- table(sampInfo[, groupCol])
    if(any(groupN == 1)) {
      message("Some sample groups contains only 1 sample. This may produce unexpected results.")
    }
  }
  return(sampInfo)

}

#' Evaluate and print translator df information
#'
#' @param df A 2 columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#'
messageMappingInfo <- function(df, sourceKey, targetKey) {

  # count source ids
  uniqueSource <- length(unique(df[, sourceKey]))
  # remove NAs
  notMapped <- is.na(df[, targetKey])
  df <- df[!notMapped, ]
  # count multimapping IDs
  sourceMultiMap <- sum(table(df[, sourceKey]) >= 2)
  targetMultiMap <- sum(table(df[, targetKey]) >= 2)
  # count target ids
  uniqueTarget <- length(unique(df[, targetKey]))
  # print info
  m <- paste("------------------------------------------------",
             paste0(sum(notMapped), " of ", uniqueSource, " input ids on the translator data frame could not be mapped."),
             paste0(sourceMultiMap, " of ", uniqueSource, " input ids on the translator data frame were mapped to 2 or more target ids."),
             paste0(targetMultiMap, " of ", uniqueTarget, " target ids on the translator data frame were mapped to 2 or more input ids."),
             "------------------------------------------------",
             paste0("Input keys were finally mapped to ", uniqueTarget, " target ids."),
             "------------------------------------------------",
             sep = "\n")
  message(m)

}

#' Validate translator data frame
#'
#' Evaluates integrity of the translator data frame, subsetting it
#' to selected columns and removing duplicated rows.
#'
#' @param df The translator data frame to validate.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#'
#' @return The validated and ammended translator data frame if no problems are found.
#'
validateTranslatorDf <- function(df, sourceKey, targetKey) {

  # subset to desired columns
  translatorDf <- df[ , c(sourceKey, targetKey)]
  # remove duplicated rows in source to target df
  uniqueDf <- unique(translatorDf)
  if(nrow(uniqueDf) != nrow(translatorDf)) {
    message("Translator data frame contains duplicated rows and will be removed.")
  }
  # print mapping info
  messageMappingInfo(uniqueDf, sourceKey, targetKey)
  # remove not mapped ids
  mapped <- !is.na(uniqueDf[, targetKey])
  uniqueDf <- uniqueDf[mapped, ]
  return(uniqueDf)

}
