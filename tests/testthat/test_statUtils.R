# load data for testing
testMat <- matrix(rnorm(9000), ncol = 9, nrow = 1000)
colnames(testMat) <- paste0(rep(c("A","B","C"), each = 3), rep(c(1,2,3), 3))
rownames(testMat) <- paste0("gene_", 1:1000)
testSampInfo <- data.frame(row.names = colnames(testMat), group = rep(c("A","B","C"), each = 3))
testSe <- SummarizedExperiment::SummarizedExperiment(assay = testMat, colData = testSampInfo)
testNaMat <- matrix(NA, ncol = 9, nrow = 1000)
testDesign <- model.matrix(~ 0 + group, testSampInfo)
colnames(testDesign) <- c("A","B","C")


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



