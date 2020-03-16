library(org.Mm.eg.db)
data("blmSubset")
data("hallmarks")
blmCounts <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, "ENTREZID", "SYMBOL", sum)
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
blmDeTable <- autoEdgeR(se = blmSe, groupColumn = "group")
blmSingle <- subset(blmDeTable, comparison == "blm-control")

test_that("Default analysis from table", {

  t <- defaultFromTable(blmDeTable, funCategories = mouseHallmarks, gseaPCutoff = 1, oraPCutoff = 1)
  t2 <- defaultFromTable(blmSingle, funCategories = mouseHallmarks, gseaPCutoff = 1, oraPCutoff = 1)
  t3 <- defaultFromTable(blmDeTable, funCategories = mouseHallmarks, useCutoff = FALSE, gseaPCutoff = 1, oraPCutoff = 1)
  t4 <- defaultFromTable(blmSingle, funCategories = mouseHallmarks, useCutoff = FALSE, gseaPCutoff = 1, oraPCutoff = 1)

  expect_equal(class(t), "list")
  expect_equal(class(t2), "list")
  expect_equal(class(t3), "list")
  expect_equal(class(t4), "list")

})
