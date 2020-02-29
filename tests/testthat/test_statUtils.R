library(biokit)

# create random matrix for testing
set.seed(149)
testMat <- matrix( runif(n = 6000), ncol = 6)

test_that("Test stat utils",{

  expect_equal(pairwiseContrasts(c("A","B","C")), c("A-B", "A-C", "B-C"))
  expect_equal(pcaToList(testMat)$result, prcomp(testMat))

})

