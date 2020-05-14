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
#' @importFrom dplyr %>% group_by summarise_all sym
#' @importFrom tibble rownames_to_column column_to_rownames
#'
translateMatrix <- function(mat, df, sourceKey, targetKey, summariseFun) {

  # if not function is passed, use first match as default function
  if(missing(summariseFun)) {
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
  # remove source id, group and summarize
  mergedDf[, sourceKey] <- NULL
  # tidy eval for targetKey (https://tidyeval.tidyverse.org/introduction.html)
  outMat <- dplyr::group_by(mergedDf, !!dplyr::sym(targetKey)) %>%
    dplyr::summarise_all(.funs = summariseFun) %>%
    tibble::column_to_rownames(var = targetKey) %>%
    as.matrix()
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


