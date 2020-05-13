library(org.Hs.eg.db)
data("sarsCovData")

testMat <- sarsCovMat[1:20,]
# prepare two translator dfs with some problems (NAs, Dups, Multi-mapping)
translatorDf1 <- data.frame(from = rownames(testMat), to = rep(letters[1:4], each = 5))
translatorDf2 <- data.frame(from = rep(rownames(testMat)[1:5], 4), to = c(letters[1:19], NA))
translatorDf2[c(21,22,23), ] <- translatorDf2[c(5,7,10),]
# prepare mat without rownames and with nas
noRowNames <- testMat
rownames(noRowNames) <- NULL
naMat <- testMat
naMat[3,3] <- NA

test_that("Matrix and data frame validators works nice",  {

  # matrix validator
  expect_true(validateMatrix(testMat))
  expect_error(validateMatrix("A"))
  expect_error(validateMatrix(translatorDf1))
  expect_error(validateMatrix(noRowNames))
  expect_message(validateMatrix(naMat))

  # data frame validator
  expect_error(validateTranslatorDf(df = 5, sourceKey = "a", targetKey = "b"))
  expect_error(validateTranslatorDf(df = "a", sourceKey = "a", targetKey = "b"))

})


test_that("Summarise matrix works nice", {

  t <- translateMatrix(mat = testMat, df = translatorDf1, sourceKey = "from", targetKey = "to")
  expect_equal(nrow(t), 4)
  t2 <-  translateMatrix(mat = testMat, df = translatorDf2, sourceKey = "from", targetKey = "to", summariseFun = sum)
  expect_equal(nrow(t2), 19)

})

test_that("Test translate matrix with annotation package", {

  t <- translateMatrixWithDb(mat = sarsCovMat, db = org.Hs.eg.db, sourceKey = "SYMBOL", targetKey = "ENTREZID", summariseFun = sum)

})
