# load data for testing
testMat <- matrix(rnorm(9000), ncol = 9, nrow = 1000)
colnames(testMat) <- paste0(rep(c("A","B","C"), each = 3), rep(c(1,2,3), 3))
rownames(testMat) <- paste0("gene_", 1:1000)
testSampInfo <- data.frame(row.names = colnames(testMat), group = rep(c("A","B","C"), each = 3))
testNaMat <- matrix(NA, ncol = 9, nrow = 1000)
testDesign <- model.matrix(~ 0 + group, testSampInfo)
colnames(testDesign) <- c("A","B","C")
testDf <- data.frame(logFc = c(0.9, 0.01, 0.1 ,2 ,-2, -0.6, -1.2, 3),
                     pAdj = c(0.02, 0.01, 0.03, 0.05, 0.04, 0.06, 0.03, 0.01),
                     comparison = rep(c("comparison_A", "comparison_B"), each = 4))

test_that("Pairwise contrasts",{

  expect_equal(pairwiseContrasts(c("A","B","C")), c("A-B", "A-C", "B-C"))

})

test_that("Prcomp results",{

  # transpose test matrix
  expect_equal(pcaToList(testMat)$result, prcomp(t(testMat)))

})

test_that("Insensitive T-Test",{

  expect_equal(nsTest(c(1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1, 1.1, 1.2)), 0.002743489)
  expect_equal(nsTest(c(NA,NA,1)), NA)
  expect_equal(nsTest(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)), NA)

})

test_that("One sample T-Test over log ratio matrix", {

  matRes <- osTestMatrix(testMat)
  naRes <- osTestMatrix(testNaMat)

  expect_equal(matRes$logFc, unname(rowMeans(testMat, na.rm = TRUE)))
  expect_equal(osTestMatrix(testMat)$pValue[1], t.test(testMat[1,])$p.value)
  expect_true(is.na(naRes$logFc)[1])
  expect_true(is.na(naRes$pValue)[1])
  expect_true(is.na(naRes$pAdj)[1])

})

test_that("Creation of the design matrix", {

  desFunResult <- designFromSampInfo(testSampInfo, "group")
  expect_equal(desFunResult, testDesign)

})

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

  expect_equal(t1$pAdj, sort(testDf$pAdj))
  expect_equal(t2$status, rep(c("Up","No change","Down"), c(3,2,3)))

})

test_that("Test annotateMultiComp", {

  t1 <- annotateMultiComparison(x = testDf, useCutoff = FALSE, n = 1, sortCol = "logFc", twoSides = TRUE, decreasing = TRUE)

  expect_equal(t1$status, rep(c("Up","No change", "No change","Down"), 2))

})
