#' Translate the input matrix rownames to the desired DB ids
#'
#' Uses a database like org.Hs.eg.db (https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) or
#' org.Mm.eg.db (https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) to translate the rownames of the input matrix
#' into the desired ids. When input ids maps to several target ids, uses the input function to resolve conflicts.
#'
#' @param x The input matrix. Should have rownames.
#' @param db The database to use for the translation.
#' @param sourceKey The source key in the DB.
#' @param targetKey The target key in the DB.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#'
#' @importFrom AnnotationDbi select
#' @export
#'
translateMatrix <- function(x, db, sourceKey, targetKey, summariseFun) {

  # get mapping df
  srcTarget <- AnnotationDbi::select(db, keys = rownames(x), keytype = sourceKey, columns = targetKey)
  # translate matrix
  translatedMat <- summariseMatrix(x, df = srcTarget, sourceKey = sourceKey, targetKey = targetKey, summariseFun = summariseFun)
  return(translatedMat)

}

#' Translate input matrix using a translator df
#'
#' Translates the rownames of the input matrix
#' into the desired ids. When input ids maps to several
#' target ids, uses the input function to resolve conflicts.
#'
#' @param x The input matrix. Should have rownames.
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#' @export
#'
#' @importFrom dplyr %>% group_by summarise_all
#' @importFrom tibble rownames_to_column column_to_rownames
#'
summariseMatrix <- function(x, df, sourceKey, targetKey, summariseFun) {

  # print mapping info
  messageMappingInfo(df, sourceKey, targetKey)
  # remove not mapped ids
  mapped <- !is.na(df[,targetKey])
  df <- df[mapped,]
  # transform mat into data frame
  xDf <- as.data.frame(x) %>% tibble::rownames_to_column("id")
  # merge with data
  mergedDf <- merge(x = df, y = xDf, by.x = sourceKey, by.y = "id")
  # remove source id, group and summarize
  mergedDf[,sourceKey] <- NULL
  # tidy eval for targetKey (https://tidyeval.tidyverse.org/introduction.html)
  outMat <- dplyr::group_by(mergedDf, !!sym(targetKey)) %>%
    dplyr::summarise_all(.funs = summariseFun) %>%
    tibble::column_to_rownames(var = targetKey) %>%
    as.matrix()
  return(outMat)

}


#' Evaluate and print translator df information
#'
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#'
messageMappingInfo <- function(df, sourceKey, targetKey) {

  # remove NAs
  notMapped <- is.na(df[,targetKey])
  df <- df[!notMapped,]
  # count multimapping IDs
  sourceMultiMap <- sum(table(df[,sourceKey]) >= 2)
  targetMultiMap <- sum(table(df[,targetKey]) >= 2)
  # count unique IDs
  uniqueSource <- length(unique(df[,sourceKey]))
  uniqueTarget <- length(unique(df[,targetKey]))
  # print info
  message(paste0(sum(notMapped), " of ", uniqueSource, " input ids could not be mapped."))
  message(paste0(sourceMultiMap, " of ", uniqueSource, " input ids were mapped to 2 or more target ids."))
  message(paste0(targetMultiMap, " of ", uniqueTarget, " target ids were mapped to 2 or more input ids."))
  message("------------------------------------------------")
  message(paste0( "Input keys were finally mapped to ", uniqueTarget, " target ids."))

}




