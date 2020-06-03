data("sarsCovData")

# prepare testing matrix
testMat <- sarsCovMat[1:1000,]
tmm <- countsToTmm(testMat)

test_that("Normalized output is a matrix", expect_equal(class(tmm),  c("matrix", "array")))
