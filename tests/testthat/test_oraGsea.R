library(org.Mm.eg.db)
data("blmSubset")
data("hallmarks")
blmCounts <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL", summariseFun = sum)
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
statDf <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 0, useFilterByExpr = FALSE)
statDf <- annotateMultiComparison(x = statDf, useCutoff = FALSE, n = 100, sortCol = "pAdj")

# add noise to avoid ties in GSEA
statDf$logFc <- jitter(statDf$logFc)

test_that("Split features and ORA", {

  feats <- splitFeatures(x = statDf)
  feats2 <- splitFeatures(x = statDf, splitCol = "comparison")

  t <- oraFromList(x = feats, funCategories = mouseHallmarks, pvalueCutoff = 1, qvalueCutoff = 1)
  t2 <- oraFromList(x = feats2, funCategories = mouseHallmarks, pvalueCutoff = 1, qvalueCutoff = 1)
  t3 <- cpResultsToDf(t2)

  expect_equal(length(t), 2)
  expect_equal(length(t2), 12)
  expect_equal(class(t3), "data.frame")

})
test_that("Ranked lists and GSEA", {

  ranked <- getRankedVector(x = statDf)
  rankedList <- getRankedVectorList(x = statDf)

  t <- gseaFromList(x = rankedList, funCategories = mouseHallmarks, pvalueCutoff = 1)
  # get results in df format
  t2 <- cpResultsToDf(t)

  expect_equal(length(t), 6)
  expect_equal(class(t2), "data.frame")

})

