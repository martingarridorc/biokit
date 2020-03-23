data("blmSubset")

test_that("edgeR from contrasts and design matrix", {

  dge <- edgeR::DGEList(counts = blmCounts, group = blmSampInfo$group)
  des <- designFromSampInfo(x = blmSampInfo, column = "group")
  con <- contrastsFromDesign(des)
  t <- edgeRFromContrasts(object = dge, desMat = des, conMat = con)

  expect_equal(class(t), "data.frame")
  # for the 6000 possible comparisons of 1000 * 6 contrasts
  expect_equal(nrow(t), 6000)

})

test_that("AutoEdgeR from SummarizedExperiment", {

  blmSe <- SummarizedExperiment::SummarizedExperiment(assays = blmCounts, colData = blmSampInfo)
  t <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 0, useFilterByExpr = FALSE)
  t2 <- autoEdgeR(se = blmSe, groupColumn = "group", minCounts = 100, useFilterByExpr = TRUE)

  expect_equal(nrow(t), 6000)
  expect_true(nrow(t2) != 6000)

})


