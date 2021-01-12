#' Translate matrix rownames
#'
#' Translates the rownames of the input matrix into the
#' desired ids using a translator data frame. When
#' input ids maps to several target ids, uses the summarise
#' function to resolve conflicts.
#'
#' @param mat Input matrix. Should have rownames.
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#'
#' @export
#'
#' @importFrom dplyr %>% group_by summarise_all sym select
#' @importFrom tibble rownames_to_column column_to_rownames
#'
translateMatrix <- function(mat, df, sourceKey, targetKey, summariseFun) {

  # if not function is passed, use first match as default function
  if(missing(summariseFun)) {
    message("------------------------------------------------")
    message("No input summarise function detected, using first match on multi-mapping situations.")
    summariseFun <- function(x) return(x[1])
  }
  # check matrix and translator df, printing mapping info
  inputMatrixCheck <- validateMatrix(mat)
  df <- validateTranslatorDf(df, sourceKey, targetKey)
  # transform mat into data frame and rownames into column
  xDf <- as.data.frame(mat) %>%
    tibble::rownames_to_column("mergingVariable")
  # merge with data matrix
  mergedDf <- merge(x = df, y = xDf, by.x = sourceKey, by.y = "mergingVariable")
  # get unique sources and keys
  uniqueSources <- df[, sourceKey][df[, sourceKey] %in% names(table(df[, sourceKey]))[table(df[, sourceKey]) == 1]]
  uniqueTargets <- df[, targetKey][df[, targetKey] %in% names(table(df[, targetKey]))[table(df[, targetKey]) == 1]]
  uniqueMatches <- mergedDf[, sourceKey] %in% uniqueSources & mergedDf[,targetKey] %in% uniqueTargets
  # remove source id, group and split matrixes
  if(sum(uniqueMatches) != 0) {
    mergedDf[, sourceKey] <- NULL
    uniqueMat <- mergedDf[uniqueMatches, ] %>%
      tibble::rownames_to_column(var = "toDiscard") %>%
      dplyr::select(-toDiscard) %>%
      column_to_rownames(var = targetKey) %>%
      as.matrix()
  } else {
    uniqueMat <- NULL
  }
  # tidy eval for targetKey (https://tidyeval.tidyverse.org/introduction.html)
  if(sum(!uniqueMatches) != 0) {
    duplicatedDf <- mergedDf[ !uniqueMatches, ]
    ids <- unique(duplicatedDf[, targetKey])
    duplicatedMat <- sapply(ids, function(id) {
      outVec <- subset(duplicatedDf, duplicatedDf[, targetKey] == id) %>%
        dplyr::select( -any_of(targetKey)) %>%
        as.matrix() %>%
        apply(., 2, FUN = summariseFun)
      return(outVec)
    }) %>%
      # transpose
      t()
    # and set new rownames
    rownames(duplicatedMat) <- ids
  } else {
    duplicatedMat <- NULL
  }
  # bind cols and return matrix
  outMat <- rbind(uniqueMat, duplicatedMat)
  return(outMat)

}

#' Translate matrix rownames using annotation package
#'
#' Uses an annotation package like
#' \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db} or
#' \href{https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db}
#' to translate the rownames of the input matrix using \link[biokit]{translateMatrix}.
#'
#' @param mat Input matrix.
#' @param db Annotation package object to use for the translation.
#' @param sourceKey Source key in the DB.
#' @param targetKey Target key in the DB.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#'
#' @importFrom AnnotationDbi select
#'
#' @export
#'
translateMatrixWithDb <- function(mat, db, sourceKey, targetKey, summariseFun) {

  # get mapping df
  srcTarget <- AnnotationDbi::select(db, keys = rownames(mat), keytype = sourceKey, columns = targetKey)
  # translate matrix
  translatedMat <- translateMatrix(mat = mat, df = srcTarget, sourceKey = sourceKey,
                                   targetKey = targetKey, summariseFun = summariseFun)
  return(translatedMat)

}


