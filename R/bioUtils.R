#' Translate the input matrix rownames to the desired DB ids
#'
#' Uses a database like \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db} or
#' \href{https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db}
#' to translate the rownames of the input matrix
#' into the desired ids. When input ids maps to several target ids, uses the input function to resolve conflicts.
#'
#' @param x Iinput matrix. Should have rownames.
#' @param db Database to use for the translation.
#' @param sourceKey Source key in the DB.
#' @param targetKey Target key in the DB.
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
    translatedMat <- summariseMatrix(x, df = srcTarget, sourceKey = sourceKey,
                                     targetKey = targetKey, summariseFun = summariseFun)
    return(translatedMat)

}

#' Translate input matrix using a translator df
#'
#' Translates the rownames of the input matrix
#' into the desired ids. When input ids maps to several
#' target ids, uses the input function to resolve conflicts.
#'
#' @param x Input matrix. Should have rownames.
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#' @export
#'
#' @importFrom dplyr %>% group_by summarise_all sym distinct
#' @importFrom tibble rownames_to_column column_to_rownames
#'
summariseMatrix <- function(x, df, sourceKey, targetKey, summariseFun) {

    # remove duplicated rows in source to target df
    df <- dplyr::distinct(df)
    # print mapping info
    messageMappingInfo(df, sourceKey, targetKey)
    # remove not mapped ids
    mapped <- !is.na(df[, targetKey])
    df <- df[mapped, ]
    # transform mat into data frame
    xDf <- as.data.frame(x) %>% tibble::rownames_to_column("id")
    # merge with data
    mergedDf <- merge(x = df, y = xDf, by.x = sourceKey, by.y = "id")
    # remove source id, group and summarize
    mergedDf[, sourceKey] <- NULL
    # tidy eval for targetKey (https://tidyeval.tidyverse.org/introduction.html)
    outMat <- dplyr::group_by(mergedDf, !!dplyr::sym(targetKey)) %>%
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
    notMapped <- is.na(df[, targetKey])
    df <- df[!notMapped, ]
    # count multimapping IDs
    sourceMultiMap <- sum(table(df[, sourceKey]) >= 2)
    targetMultiMap <- sum(table(df[, targetKey]) >= 2)
    # count unique IDs
    uniqueSource <- length(unique(df[, sourceKey]))
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

#' Over-representation analysis over a list of features
#'
#' Uses a over-representation strategy (Fisher Exact T-Test) to
#' evaluate the over-represented categories in a list of features. Uses
#' the \link[clusterProfiler]{enricher} function to perform the analyses.
#'
#' @param x The list of features to be analysed.
#' @param funCategories Two-columns data frame containing the functional categories info.
#' First column should correspond to categories and the second to features.
#' Passed to \link[clusterProfiler]{enricher} TERM2GENE argument.
#' @param minFeatures Minimum features to evaluate.
#' The function will filter all list elements with a lower number of features than this value.
#' @param ... Rest of arguments passed to \link[clusterProfiler]{enricher}
#'
#' @return A list of \link[clusterProfiler]{enricher} results, one per list element.
#'
#' @export
#'
#' @importFrom clusterProfiler enricher
#'
oraFromList <- function(x, funCategories, minFeatures = 10, ...) {

    # remove list elements lower than minFeatures
    l <- unlist(lapply(x, length))
    keep <- unlist(l) >= minFeatures
    x <- x[keep]
    # perform enrichment analysis in all elements of the list
    result <- lapply(x, function(y) enricher(gene = y, TERM2GENE = funCategories, ...))
    return(result)

}

#' Transform list of clusterProfiler results into tidy data frame
#'
#' Binds the rows of a list of clusterProfiler results annotating the
#' resulting data frame with the list names.
#'
#' @param x The list of clusterProfiler results
#' @param splitName The name for the new id column.
#'
#' @return A data frame with the clusterProfiler results in a tidy format.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#'
cpResultsToDf <- function(x, splitName = "comparison") {

    dfList <- lapply(x, function(y)y@result)
    df <- dplyr::bind_rows(dfList, .id = splitName)
    return(df)

}


#' Perform Gene Set Enrichment Analysis from list of ranked vectors
#'
#' Uses \link[clusterProfiler]{GSEA} to perform the gene set enrichment analysis
#' over a list of ranked vectors.
#'
#' @param x The list of ranked vectors.
#' @param funCategories Two-columns data frame containing the functional categories info.
#' First column should correspond to categories and the second to features.
#' Passed to \link[clusterProfiler]{GSEA} TERM2GENE argument.
#' @param ... Rest of arguments passed to \link[clusterProfiler]{GSEA}.
#'
#' @return A list of \link[clusterProfiler]{GSEA} results.
#'
#' @export
#'
#' @importFrom clusterProfiler GSEA
#'
gseaFromList <- function(x, funCategories, ...) {

    result <- lapply(x, function(y) GSEA(geneList = y, TERM2GENE = funCategories, ...))
    return(result)

}
