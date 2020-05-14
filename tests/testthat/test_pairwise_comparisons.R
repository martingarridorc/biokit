data("sarsCovData")

# prepare subset version of data to test
counts <- sarsCovMat[1:20,]
absTmms <- abs(countsToTmm(counts))

genericRes <- autoPairwiseMatrixTest(mat = absTmms, sampInfo = sarsCovSampInfo, groupCol = "group")

test_that("Generic comparison", {

  expect_equal(class(genericRes), "data.frame")

})

genericIndividualRes <- singleMatrixTest(mat = absTmms)

test_that("Generic individual comparison", {

  expect_equal(class(genericIndividualRes), "data.frame")

})

limmaRes <- autoLimmaComparison(mat = absTmms, sampInfo = sarsCovSampInfo, groupCol = "group")

test_that("Limma comparison", {

  expect_equal(class(limmaRes), "data.frame")

})

edgerRes <- autoEdgerComparison(counts = counts, sampInfo = sarsCovSampInfo, groupCol = "group", filterByExpr = TRUE)

test_that("EdgeR comparison", {

  expect_equal(class(edgerRes), "data.frame")

})

