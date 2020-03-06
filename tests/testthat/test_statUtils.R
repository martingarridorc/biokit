library(biokit)

# prepare data for testing
testMat <- matrix(rnorm(6000), ncol = 6, nrow = 1000)
testNaMat <- matrix(NA, ncol = 6, nrow = 1000)

test_that("Pairwise contrasts",{

  expect_equal(pairwiseContrasts(c("A","B","C")), c("A-B", "A-C", "B-C"))

})

test_that("Prcomp results",{

  expect_equal(pcaToList(cbdMat)$result, prcomp(cbdMat))

})

test_that("Insensitive T-Test",{

  expect_equal(nsTest(c(1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1)), NA)
  expect_equal(nsTest(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)), NA)

})

test_that("One sample T-Test over log ratio matrix",{

  matRes <- osTestMatrix(testMat)
  naRes <- osTestMatrix(testNaMat)

  expect_equal(matRes$logFc, rowMeans(testMat, na.rm = TRUE))
  expect_equal(osTestMatrix(testMat)$pValue[1], t.test(testMat[1,])$p.value)
  expect_true(is.na(naRes$logFc)[1])
  expect_true(is.na(naRes$pValue)[1])
  expect_true(is.na(naRes$pAdj)[1])

})

