data("blmSubset")
testMat <- countsToTmm(x = blmCounts)
testMat[4,] <- NA
testMat[8,] <- 10
unNamedMat <- testMat
rownames(unNamedMat) <- NULL
badContrasts <- c("this-is-not", "a-real-contrast")
badSampInfo <- blmSampInfo
rownames(badSampInfo) <- letters[1:12]
testSe <- SummarizedExperiment::SummarizedExperiment(assays = testMat, colData = blmSampInfo)

test_that("Single matrix One Sample T-Test", {

  t <- matrixSingleTest(testMat)
  t2 <- matrixSingleTest(unNamedMat)

  expect_equal(t$logFc, unname(rowMeans(testMat, na.rm = TRUE)))
  expect_true(is.na(t$pValue[4]))
  expect_true(is.na(t$t[4]))
  expect_true(is.na(t$pValue[8]))
  expect_true(is.na(t$t[8]))
  expect_equal(t2$id, 1:nrow(unNamedMat))

})

test_that("Pairwise matrix comparison", {

  t <- matrixPairwiseTest(x = testMat, blmSampInfo, "group", numerator = "blm", denominator = "control")
  t2 <- matrixPairwiseTest(x = unNamedMat, blmSampInfo, "group", numerator = "blm", denominator = "control")

  expect_equal(t$logFc[1], log2(mean(testMat[ 1, 4:6])/mean(testMat[ 1, 1:3])))
  expect_true(is.na(t$pValue[4]))
  expect_true(is.na(t$t[4]))
  expect_equal(t2$id, 1:nrow(unNamedMat))
  expect_error(matrixPairwiseTest(x = testMat, badSampInfo, "group", numerator = "blm", denominator = "control"))

})

test_that("Pairwise matrix comparison from contrasts", {

  t <- matrixTestFromContrasts(x = testMat, blmSampInfo, "group", contrasts = c("blm-control", "blm_aja-control"))
  expect_equal(nrow(t), 2000)
  expect_error(matrixTestFromContrasts(x = testMat, blmSampInfo, "group", contrasts = badContrasts))

})

test_that("Auto matrix comp", {

  t <- autoMatrixPairwise(se = testSe, groupCol = "group")
  expect_equal(nrow(t), 6000)

})
