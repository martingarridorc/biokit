# prepare data for testing
intMat <- matrix(1:60, ncol = 6, nrow = 10)
rownames(intMat) <- letters[1:10]
testDf <- data.frame(from = rownames(intMat), to = rep(c(NA, "A", "B", "C", "D", "E"), c(2, 2, 1, 3, 1, 1)), stringsAsFactors = FALSE)

test_that("Message outputs nice", {
    
    m <- evaluate_promise(biokit:::messageMappingInfo(testDf, "from", "to"))$messages
    expect_equal(substr(m, 245, 290), "Input keys were finally mapped to 5 target ids")
    
})

test_that("Summarise matrix", {
    
    res <- summariseMatrix(x = intMat, df = testDf, sourceKey = "from", targetKey = "to", summariseFun = sum)
    expect_type(res, "integer")
    expect_equal(class(res), "matrix")
    expect_true(all(rownames(res) %in% testDf$to))
    expect_equal(as.numeric(rowSums(res)[3]), 576)
    
})
