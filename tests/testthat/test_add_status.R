data("sarsCovData")

# prepare subset version of data to test
counts <- sarsCovMat[1:200,]
absTmms <- abs(countsToTmm(counts))
genericRes <- autoPairwiseMatrixTest(mat = absTmms, sampInfo = sarsCovSampInfo, groupCol = "group")

test_that("Add status by cutoff", {

  t <- addStatusByCutoff(resDf = genericRes, splitUpDown = TRUE)
  t2 <- addStatusByCutoff(resDf = genericRes, splitUpDown = FALSE)

  expect_equal(sum(t$status %in% c("Up", "Down")), expected = sum(t2$status == "Significant"))

})

test_that("Add status by topN", {

  t <- addStatusByRank(resDf = genericRes, rankCol = "logFc", topN = 20, splitUpDown = TRUE, absolute = TRUE)
  t2 <- addStatusByRank(resDf = genericRes, rankCol = "logFc", topN = 20, splitUpDown = FALSE, absolute = FALSE)

  expect_equal(sum(t$status %in% c("Up", "Down")), expected = sum(t2$status == "Significant")*2)

})
