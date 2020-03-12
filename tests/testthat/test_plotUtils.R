# prepare data for testing
testDf <- data.frame(logFc = c(0.9, 0.01, 0.1, 2, -2, -0.6, -1.2, 3), pAdj = c(0.02, 0.01, 0.03, 0.05, 0.04, 0.06, 0.03, 0.01),
                     comparison = rep(c("comparison_A", "comparison_B"), each = 4),
                     feature = letters[1:8])
testDf <- annotateByCutoff(testDf)

test_that("Defaults volcano plot", {

    p <- defaultVolcano(x = testDf)
    expect_equal(class(p), c("gg","ggplot"))

})
