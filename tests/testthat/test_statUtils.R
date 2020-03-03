library(biokit)

# create random matrix for testing
set.seed(149)
testMat <- matrix( runif(n = 6000), ncol = 6)
testNaMat <- matrix(NA, ncol = 6, nrow = 1000)

test_that("Pairwise contrasts",{

  expect_equal(pairwiseContrasts(c("A","B","C")), c("A-B", "A-C", "B-C"))

})

test_that("Prcomp results",{

  expect_equal(pcaToList(testMat)$result, prcomp(testMat))

})

test_that("Insensitive T-Test",{

  expect_equal(nsTest(c(1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1)), NA)
  expect_equal(nsTest(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)), NA)

})

test_that("One sample T-Test over log ratio matrix",{

  expect_equal(osTestMatrix(testMat)$logFc, rowMeans(testMat, na.rm = TRUE))
  expect_equal(osTestMatrix(testMat)$pValue[1], t.test(testMat[1,])$p.value)
  expect_true(is.na(osTestMatrix(testNaMat)$logFc)[1])
  expect_true(is.na(osTestMatrix(testNaMat)$pValue)[1])
  expect_true(is.na(osTestMatrix(testNaMat)$pAdj)[1])

})

