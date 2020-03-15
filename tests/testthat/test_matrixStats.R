testMat <- matrix(rnorm(12000), ncol = 12, nrow = 1000)
testMat[4,] <- NA
testMat[8,] <- 5
namedMat <- testMat
rownames(namedMat) <- paste0("id", 1:1000)

test_that("Single matrix One Sample T-Test", {

  t <- matrixSingleTest(testMat)
  t2 <- matrixSingleTest(namedMat)

  expect_equal(t$logFc, rowMeans(testMat, na.rm = TRUE))
  expect_true(is.na(t$pValue[4]))
  expect_true(is.na(t$t[4]))
  expect_true(is.na(t$pValue[8]))
  expect_true(is.na(t$t[8]))
  expect_equal(t2$id, paste0("id", 1:1000))

})

