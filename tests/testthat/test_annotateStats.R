library(dplyr)
data("blmSubset")
blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
statDf <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 0, useFilterByExpr = FALSE)

p <- case_when(statDf$logFc <= -1 & statDf$pAdj <= 0.05 ~ "Down",
               statDf$logFc >= 1 & statDf$pAdj <= 0.05 ~ "Up",
               TRUE ~ "No change")

test_that("Annotate by cutoff", {

  t <- annotateByCutoff(statDf)
  t2 <- annotateByCutoff(statDf, sigCutoff =  NULL)
  t3 <- annotateByCutoff(statDf, splitUpDown = FALSE)
  t4 <- annotateByCutoff(statDf, splitUpDown = FALSE, sigCutoff = NULL)
  t5  <- annotateByCutoff(statDf, splitUpDown = FALSE, metricCutoff =  NULL)

  expect_equal(t$status, p)
  expect_equal(t2$status == "Up", statDf$logFc >= 1)
  expect_equal(t3$status == "Significant", statDf$pAdj <= 0.05 & abs(statDf$logFc) >= 1)
  expect_equal(t4$status == "Significant", abs(statDf$logFc) >= 1)
  expect_equal(t5$status == "Significant", abs(statDf$pAdj) <= 0.05)

})

test_that("Annotate top N features", {

  byPadj <- arrange(statDf, pAdj)
  byLogFc <- arrange(statDf, desc(logFc))

  t <- annotateTopN(statDf, n = 50, sortCol = "pAdj")
  t2 <- annotateTopN(statDf, n = 50, sortCol = "pAdj", splitUpDown = FALSE)
  t3 <- annotateTopN(statDf, n = 50, sortCol = "logFc", splitUpDown = TRUE, twoSides = TRUE, decreasing = TRUE)

  expect_equal(sum(t$status == "Up" | t$status == "Down"), 100)
  expect_equal(t2$feature[t2$status == "Significant"], byPadj$feature[1:50])
  expect_equal(t3$feature[t3$status == "Up"], byLogFc$feature[1:50])

})

test_that("Annotate multi comparison", {

  t <- annotateMultiComparison(statDf, useCutoff = TRUE, compCol = "comparison")
  t2 <- annotateMultiComparison(statDf, useCutoff = FALSE, compCol = "comparison", n = 50, sortCol = "pAdj")

  expect_equal(t$status, p)
  expect_equal(sum(t2$status == "Up" | t2$status == "Down"), 600)

})
