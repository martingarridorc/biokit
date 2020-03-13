
testDf <- data.frame(logFc = c(0.9, 0.01, 0.1, 2, -2, -0.6, -1.2, 3), pAdj = c(0.02, 0.01, 0.03, 0.05, 0.04, 0.06, 0.03, 0.01), comparison = rep(c("comparison_A", "comparison_B"),
                                                                                                                                                 each = 4))

test_that("Test annotateByCutoff", {

  t1 <- annotateByCutoff(testDf)
  t2 <- annotateByCutoff(testDf, splitUpDown = FALSE)

  expect_equal(t1[8, "status"], "Up")
  expect_equal(t1[5, "status"], "Down")
  expect_equal(t2[8, "status"], "Significant")
  expect_equal(t2[5, "status"], "Significant")

})

test_that("Test annotateTopN", {

  t1 <- annotateTopN(testDf, n = 3, sortCol = "pAdj")
  t2 <- annotateTopN(testDf, n = 3, twoSides = TRUE, sortCol = "logFc")
  t3 <- annotateTopN(testDf, n = 3, splitUpDown = FALSE, sortCol = "pAdj")

  expect_equal(t1$pAdj, sort(testDf$pAdj))
  expect_equal(t2$status, rep(c("Up", "No change", "Down"), c(3, 2, 3)))
  expect_equal(t3$status, rep(c("Significant", "No change"), c(3, 5)))

})

test_that("Test annotateMultiComp", {

  t1 <- annotateMultiComparison(x = testDf, useCutoff = FALSE, n = 1, sortCol = "logFc", twoSides = TRUE, decreasing = TRUE)
  t2 <- annotateMultiComparison(x = testDf, useCutoff = TRUE)

  expect_equal(t1$status, rep(c("Up", "No change", "No change", "Down"), 2))
  expect_equal(t2$status[4], "Up")

})
