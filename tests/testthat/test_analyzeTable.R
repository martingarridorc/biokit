library(org.Mm.eg.db)
data("blmSubset")
data("hallmarks")
blmCounts <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, "ENTREZID", "SYMBOL", sum)
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
blmDeTable <- autoEdgeR(se = blmSe, groupColumn = "group")
blmSingle <- subset(blmDeTable, comparison == "blm-control")

test_that("Default analysis from table", {

  t <- analyzeDeTable(blmDeTable, funCategories = mouseHallmarks)
  t2 <- analyzeDeTable(blmSingle, funCategories = mouseHallmarks)
  t3 <- analyzeDeTable(blmDeTable, funCategories = mouseHallmarks, useCutoff = FALSE)
  t4 <- analyzeDeTable(blmSingle, funCategories = mouseHallmarks, useCutoff = FALSE)

  expect_equal(class(t), "list")
  expect_equal(class(t2), "list")
  expect_equal(class(t3), "list")
  expect_equal(class(t4), "list")

})
