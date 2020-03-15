data("blmSubset")
testMat <- matrix(rnorm(12000), ncol = 12, nrow = 1000)
rownames(testMat) <- rownames(blmCounts)
colnames(testMat) <- colnames(blmCounts)
des <- designFromSampInfo(x = blmSampInfo, column = "group")
con <- contrastsFromDesign(des)
testSe <- SummarizedExperiment::SummarizedExperiment(assays = testMat, colData = blmSampInfo)

test_that("Limma from contrasts and design", {

  t <- limmaDfFromContrasts(x = testMat, desMat = des, conMat = con)

  expect_equal(class(t), "data.frame")
  expect_equal(nrow(t), 6000)

})

test_that("AutoLimma from SE", {

  t <- autoLimma(se = testSe, groupCol = "group")
  expect_equal(class(t), "data.frame")
  expect_equal(nrow(t), 6000)

})
