library(org.Mm.eg.db)
data("blmSubset")
testMat <- matrix(1:60, ncol = 6, nrow = 10)
rownames(testMat) <- letters[1:10]
translatorDf <- data.frame(from = rownames(testMat), to = rep(c(NA, "A", "B", "C", "D", "E"), c(2, 2, 1, 3, 1, 1)), stringsAsFactors = FALSE)

test_that("Summarise matrix works nice", {

  t <- summariseMatrix(x = testMat, df = translatorDf, sourceKey = "from", targetKey = "to", summariseFun = sum)
  expect_equal(as.numeric(rowSums(t)[3]), 576)
  expect_equal(class(t), "matrix")
  expect_equal(typeof(t), "integer")

})

test_that("Test translate matrix with db", {

  t <- translateMatrix(x = blmCounts, db = org.Mm.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL", summariseFun = sum)
  expect_equal(class(t), "matrix")
  expect_equal(typeof(t), "integer")

})
