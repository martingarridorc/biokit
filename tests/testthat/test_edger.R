data("blmSubset")
data("cbdSubset")

test_that("edgeR from contrasts and design matrix", {

  dge <- edgeR::DGEList(counts = blmCounts, group = blmSampInfo$group)
  des <- designFromSampInfo(x = blmSampInfo, column = "group")
  con <- contrastsFromDesign(des)
  t <- edgeRFromContrasts(object = dge, desMat = des, conMat = con)

  expect_equal(class(t), "data.frame")
  # for the 6000 possible comparisons of 1000 * 6 contrasts
  expect_equal(nrow(t), 6000)

})

test_that("edgeR from SummarizedExperiment", {

  blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
  cbdSe <- SummarizedExperiment::SummarizedExperiment(assays = cbdCounts, colData = cbdSampInfo)
  t <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 0, useFilterByExpr = FALSE)
  t2 <- autoEdgeR(se = cbdSe, groupColumn = "group",  minCounts = 0, useFilterByExpr = FALSE)
  t3 <- autoEdgeR(se = cbdSe, groupColumn = "group",  minCounts = 100, useFilterByExpr = TRUE)

  expect_equal(nrow(t), 6000)
  expect_equal(nrow(t2), 1000)
  expect_true(nrow(t3) != 6000)

})


