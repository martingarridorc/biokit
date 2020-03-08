library(testthat)
library(biokit)
library(edgeR)
library(dplyr)

# prepare cbd data for testing
data("cbdData")
cbdMat <- cbdMat[rowSums(cbdMat) > 10,]
cbdSe <- SummarizedExperiment::SummarizedExperiment(assay = log(cbdMat + 1), colData = cbdSampInfo)

# prepare random data for testing
testMat <- matrix(rnorm(9000), ncol = 9, nrow = 1000)
colnames(testMat) <- paste0(rep(c("A","B","C"), each = 3), rep(c(1,2,3), 3))
rownames(testMat) <- paste0("gene_", 1:1000)
testSampInfo <- data.frame(row.names = colnames(testMat), group = rep(c("A","B","C"), each = 3))
testSe <- SummarizedExperiment::SummarizedExperiment(assay = testMat, colData = testSampInfo)


testNaMat <- matrix(NA, ncol = 9, nrow = 1000)


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

test_that("Creation of the design and contrast matrix", {

  testDesign <- model.matrix(~ 0 + group, testSampInfo)
  colnames(testDesign) <- c("A","B","C")

  desFunResult <- designFromSampInfo(testSampInfo, "group")
  conFunResult <- contrastsFromDesign(testDesign)

  expect_equal(desFunResult, testDesign)
  expect_equal(colnames(conFunResult), c("A-B", "A-C", "B-C"))
  expect_equal(rownames(conFunResult), c("A", "B", "C"))
  expect_equal(conFunResult["C", "A-C"] , -1)

})

test_that("Integrative test for limma comparisons", {

  testDesign <- designFromSampInfo(testSampInfo, "group")
  testCon <- contrastsFromDesign(testDesign)
  limmaRes <- limmaDfFromContrasts(testMat, testDesign, testCon)
  cbdAutoRes <- autoLimma(se = cbdSe, groupColumn = "group")

  expect_equal(nrow(limmaRes), 3000)
  expect_equal(colnames(limmaRes), c("comparison", "feature", "logFc", "AveExpr", "t", "pValue", "pAdj", "B"))
  expect_equal(unique(limmaRes$comparison), c("A-B", "A-C", "B-C"))
  expect_equal(colnames(cbdAutoRes), c("comparison", "feature", "logFc", "AveExpr", "t", "pValue", "pAdj", "B"))
  expect_equal(unique(cbdAutoRes$comparison), "cbd-control")

})

test_that("TMM norm", {

  cbdTmms <- DGEList(counts = cbdMat) %>% calcNormFactors() %>% cpm(y = ., log = TRUE, prior.count = 3)

  expect_equal(countsToTmm(x = cbdMat), cbdTmms)

})


