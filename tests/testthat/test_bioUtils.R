library(org.Hs.eg.db)
library(org.Mm.eg.db)
# prepare data for testing
data("cbdData")
cbdMat <- cbdMat[1:100, ]
intMat <- matrix(1:60, ncol = 6, nrow = 10)
rownames(intMat) <- letters[1:10]
translatorDf <- data.frame(from = rownames(intMat), to = rep(c(NA, "A", "B", "C", "D", "E"), c(2, 2, 1, 3, 1, 1)), stringsAsFactors = FALSE)
data("blmData")
blmMat <- blmMat[1:1000,]
blmMat <- translateMatrix(x = blmMat, db = org.Mm.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL",summariseFun = sum)
testSe <- SummarizedExperiment::SummarizedExperiment(blmMat, colData = blmSampInfo)
res <- autoEdgeR(testSe, groupColumn = "group")
res <- annotateMultiComparison(x = res, useCutoff = FALSE, n = 50, sortCol = "pAdj")
data("msigdbHallmarks")

test_that("Message outputs nice", {

    m <- evaluate_promise(biokit:::messageMappingInfo(translatorDf, "from", "to"))$messages
    expect_equal(substr(m, 245, 290), "Input keys were finally mapped to 5 target ids")

})

test_that("Summarise matrix", {

    res <- summariseMatrix(x = intMat, df = translatorDf, sourceKey = "from", targetKey = "to", summariseFun = sum)
    expect_type(res, "integer")
    expect_equal(class(res), "matrix")
    expect_true(all(rownames(res) %in% translatorDf$to))
    expect_equal(as.numeric(rowSums(res)[3]), 576)

})

test_that("Translate matrix", {

    t <- translateMatrix(x = cbdMat, db = org.Hs.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL", summariseFun = sum)
    expect_equal(class(t), "matrix")

})

test_that("Ora from list", {

    featList <- biokit::splitFeatures(x = res, splitCol = "comparison")
    t <- oraFromList(x = featList, funCategories = hallmarks)
    df <- cpResultsToDf(t)
    expect_equal(class(t), "list")
    expect_equal(class(df), "data.frame")

})

test_that("GSEA from list", {

    rList <- biokit::getRankedVectorList(x = res)
    t <- gseaFromList(x = rList, funCategories = hallmarks)
    df <- cpResultsToDf(t)
    expect_equal(class(t), "list")
    expect_equal(class(df), "data.frame")

})


