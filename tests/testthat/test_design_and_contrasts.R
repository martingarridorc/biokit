data("sarsCovData")

desMat <- designFromSampInfo(sarsCovSampInfo, groupCol = "group")

test_that("Design creator" , {

  expect_equal(class(desMat), c("matrix", "array"))

})

conMat <- contrastsFromDesign(desMat)

test_that("Contrast creator from design" , {

  expect_equal(class(desMat), c("matrix", "array"))

})

levs <- c("A","B","C", "C", NA)
contrasts <- createPairwiseContrasts(levs)

test_that("Contrast creator from design" , {

  expect_equal(contrasts, c("C-B","C-A","B-A"))

})

conMat <- contrastsFromDesign(desMat)

test_that("Contrast creator from design" , {

  expect_equal(class(desMat), c("matrix", "array"))

})
