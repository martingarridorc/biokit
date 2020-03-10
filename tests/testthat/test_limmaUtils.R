# prepare cbd data for testing
data("blmData")
blmMat <- blmMat[rowSums(blmMat) > 25,]
blmMat <- blmMat[1:50,]
desMat <- designFromSampInfo(blmSampInfo, "group")
conMat <- contrastsFromDesign(desMat)
testSe <- SummarizedExperiment::SummarizedExperiment(blmMat, colData = blmSampInfo)

test_that("Contrast matrix creation", {

  expect_equal(unname(conMat[,6]), c(0,0,1,-1))

})

test_that("Limma Df from contrast and test design", {

  result <- limmaDfFromContrasts(x = log(blmMat), desMat = desMat, conMat = conMat)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 300)

})

test_that("AutoLimma", {

  result <- autoLimma(se = testSe, "group")
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 300)

})
