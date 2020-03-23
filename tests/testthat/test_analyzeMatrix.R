library(org.Mm.eg.db)
data("blmSubset")
data("hallmarks")
blmCounts <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, "ENTREZID", "SYMBOL", sum)
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
blmNorm <- SummarizedExperiment::SummarizedExperiment(assays = jitter(countsToTmm(blmCounts)), colData = blmSampInfo)

test_that("Default analysis from table", {

  t <- analyzeCounts(se = blmSe, groupCol = "group", funCategories = mouseHallmarks)
  t2 <- analyzeMatrix(se = blmNorm, groupCol = "group", funCategories = mouseHallmarks, useLimma = TRUE)
  t3 <- analyzeMatrix(se = blmNorm, groupCol = "group", funCategories = mouseHallmarks, useLimma = FALSE, useCutoff = FALSE)

  expect_equal(class(t), "list")
  expect_equal(class(t2), "list")
  expect_equal(class(t3), "list")

})
