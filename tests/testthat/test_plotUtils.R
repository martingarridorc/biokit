# prepare data for testing
testDf <- data.frame(logFc = c(0.9, 0.01, 0.1, 2, -2, -0.6, -1.2, 3),
                     pAdj = c(0.02, 0.01, 0.03, 0.05, 0.04, 0.06, 0.03, 0.01),
                     feature = letters[1:8],
                     comparison = rep(c("comparison_A", "comparison_B"), each = 4))
testDf <- annotateByCutoff(testDf)
testMat <- matrix(rnorm(9000), ncol = 9, nrow = 1000)
colnames(testMat) <- paste0(rep(c("A", "B", "C"), each = 3), rep(c(1, 2, 3), 3))
testSampInfo <- data.frame(row.names = colnames(testMat), group = rep(c("A", "B", "C"), each = 3))
badSampInfo <-  data.frame(row.names = letters[1:9], group = rep(c("A", "B", "C"), each = 3))
data("blmData")
blmMat <- blmMat[1:1000,]
blmMat <- translateMatrix(x = blmMat, db = org.Mm.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL",summariseFun = sum)
testSe <- SummarizedExperiment::SummarizedExperiment(blmMat, colData = blmSampInfo)
res <- autoEdgeR(testSe, groupColumn = "group")
res <- annotateMultiComparison(x = res, useCutoff = FALSE, n = 50, sortCol = "pAdj")
data("msigdbHallmarks")
featList <- biokit::splitFeatures(x = res, splitCol = "comparison")
oraRes <- oraFromList(x = featList, funCategories = hallmarks)
oraDf <- cpResultsToDf(oraRes)

test_that("Default volcano plot", {

    p <- defaultVolcano(x = testDf, compCol = "comparison")
    expect_equal(class(p), c("gg", "ggplot"))

})

test_that("Defaults pca plot", {

    res <- pcaToList(testMat)
    p <- defaultPcaPlot(x = res, sampInfoDf = testSampInfo, groupCol = "group")
    expect_equal(class(p), c("gg", "ggplot"))
    expect_error(defaultPcaPlot(x = res, sampInfoDf = badSampInfo))

})

test_that("Defaults ora plot", {

    expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = TRUE)), c("gg", "ggplot"))
    expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = FALSE)), c("gg", "ggplot"))
    expect_equal(class(defaultOraPlot(x = oraDf, splitStatus = TRUE, pCutoff = 0.05)), c("gg", "ggplot"))

})


