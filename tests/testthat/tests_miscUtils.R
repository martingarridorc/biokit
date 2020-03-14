testDf <- data.frame(logFc = c(0.9, 5, 0.1, 2, -2, -0.6, -1.2, 3),
                     pAdj = c(0.02, 0.07, 0.03, 0.05, 0.04, 0.06, 0.03, 0.01),
                     feature = letters[1:8],
                     comparison = rep(c("comparison_A", "comparison_B"), each = 4))


test_that("Test annotateByCutoff", {

  expect_equal(class(annotateByCutoff(testDf)), "data.frame")
  expect_equal(class(annotateByCutoff(testDf, splitUpDown = FALSE)), "data.frame")
  expect_equal(class(annotateByCutoff(testDf, splitUpDown = TRUE, sigCutoff = NULL)), "data.frame")
  expect_equal(class(annotateByCutoff(testDf, splitUpDown = FALSE, sigCutoff = NULL)), "data.frame")
  expect_equal(class(annotateByCutoff(testDf, splitUpDown = FALSE, metricCutoff = NULL)), "data.frame")

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

test_that("Split by label", {

  annotatedTestDf <- annotateTopN(testDf, n = 3, sortCol = "pAdj")
  t1 <- splitFeatures(annotatedTestDf)
  t2 <- splitFeatures(annotatedTestDf, splitCol = "comparison")
  expect_equal(class(t1), "list")
  expect_equal(class(t2), "list")

})
