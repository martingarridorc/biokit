library(org.Hs.eg.db)
data("sarsCovData")

# prepare matrix to test
testMat <- sarsCovMat[1:100,]
translatorDf <- data.frame(from = rownames(testMat), to = as.character(c(1:90, 1:10)))

t <- translateMatrix(mat = testMat, df = translatorDf, sourceKey = "from", targetKey = "to")
t2 <- translateMatrixWithDb(testMat, db = org.Hs.eg.db, sourceKey = "SYMBOL", targetKey = "ENTREZID")

test_that("Translator functions work nice", {

  expect_equal(nrow(t), 90)
  expect_equal(class(t2), c("matrix", "array"))

})
