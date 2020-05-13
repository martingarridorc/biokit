#' Evaluate and print translator df information
#'
#' @param df A 2-columns translator data frame with source to target ids.
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
             paste0(sum(notMapped), " of ", uniqueSource, " input ids could not be mapped."),
             paste0(sourceMultiMap, " of ", uniqueSource, " input ids were mapped to 2 or more target ids."),
             paste0(targetMultiMap, " of ", uniqueTarget, " target ids were mapped to 2 or more input ids."),
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
#' @return The validated and ammended translator
#' data frame if no problems are found.
#'
validateTranslatorDf <- function(df, sourceKey, targetKey) {

  # check data frame class
  dfClass <- class(df)
  if(dfClass != "data.frame") {
    stop("Input is not a data frame.")
  }
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
