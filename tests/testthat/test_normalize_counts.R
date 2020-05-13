library(org.Hs.eg.db)
data("sarsCovData")

testMat <- sarsCovMat[1:20,]
testNaMat <- testMat
testNaMat[1, ] <- NA

test_that("Counts to TMM function works nice",  {

  expect_type(countsToTmm(testMat), "double")
  expect_error(countsToTmm(testNaMat))

})
