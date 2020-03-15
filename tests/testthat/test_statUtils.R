testMat <- matrix(rnorm(12000), ncol = 12, nrow = 1000)

test_that("Pairwise contrasts", {

  expect_equal(pairwiseContrasts(c("A", "B", "C")), c("A-B", "A-C", "B-C"))

})

test_that("Prcomp results", {

  # transpose test matrix
  expect_equal(pcaToList(testMat)$result, prcomp(t(testMat)))

})
