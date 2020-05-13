library(org.Mm.eg.db)
library(dplyr)
data("blmSubset")
data("hallmarks")
blmCounts <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL", summariseFun = sum)
blmTmms <- countsToTmm(blmCounts)
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
statDf <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 0, useFilterByExpr = FALSE) %>%
  annotateMultiComparison(x = ., useCutoff = FALSE, n = 100, sortCol = "pAdj") %>%
  mutate(logFc = jitter(logFc))
oraDf <- splitFeatures(statDf) %>% oraFromList(x = ., funCategories = mouseHallmarks) %>% cpResultsToDf()
gseaDf <- getRankedVectorList(statDf) %>% gseaFromList(x = ., funCategories = mouseHallmarks,pvalueCutoff = 1 ) %>%  cpResultsToDf()

test_that("Default pca plot", {

  res <- pcaToList(blmTmms)
  p <- defaultPcaPlot(x = res, sampInfoDf = blmSampInfo, groupCol = "group")
  expect_equal(class(p), c("gg", "ggplot"))
  expect_error(defaultPcaPlot(x = res, sampInfoDf = mutate(blmSampInfo, group = letters[1:12])))

})


test_that("Default volcano plot", {

  p <- defaultVolcano(x = statDf, compCol = "comparison")
  expect_equal(class(p), c("gg", "ggplot"))

})


test_that("Default ora plot", {

  expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = TRUE)), c("gg", "ggplot"))
  expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = FALSE)), c("gg", "ggplot"))
  expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = TRUE, pCutoff = 1)), c("gg", "ggplot"))
  expect_error(defaultOraPlot(x = oraDf, splitStatus = TRUE, pCutoff = -1))

})

test_that("Defaults gsea plot", {

  expect_equal(class(defaultGseaDotPlot(x = gseaDf, splitById = FALSE)), c("gg", "ggplot"))
  expect_equal(class(defaultGseaDotPlot(x = gseaDf, splitById = TRUE)), c("gg", "ggplot"))
  expect_equal(class(defaultGseaDotPlot(x = gseaDf, splitById = TRUE, pCutoff = 1)), c("gg", "ggplot"))
  expect_error(defaultGseaDotPlot(x = gseaDf, splitById = TRUE, pCutoff = -1))

})
