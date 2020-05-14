data("sarsCovData")

# subset to first 20 columns to validate
testMat <- countsToTmm(sarsCovMat)[1:20,]

test_that("Matrix validator", {

  expect_error(validateMatrix(data.frame(a = "a")))
  expect_true(validateMatrix(testMat))
  noRnames <- testMat
  rownames(noRnames) <- NULL
  expect_error(validateMatrix(noRnames, checkRowNames = TRUE))
  naMat <- testMat
  naMat[4,] <- NA
  expect_true(validateMatrix(naMat, allowNa = TRUE))
  expect_error(validateMatrix(naMat, allowNa = FALSE))

})

sarsCovSampInfo[ ,"wrongNames"] <- rep(c("123", "this-is-wrong", "this-is-alone"), c(6,5,1))

test_that("Sample information validator", {

  expect_equal(class(validateSampInfo(sarsCovSampInfo, groupCol = "group")), "data.frame")
  expect_warning(validateSampInfo(sarsCovSampInfo, groupCol = "wrongNames"))
  expect_equal(class(validateSampInfo(sampInfo = sarsCovSampInfo, groupCol = "group", mat = sarsCovMat)), "data.frame")
  expect_error(class(validateSampInfo(sampInfo = sarsCovSampInfo[1:3,], groupCol = "group", mat = sarsCovMat)))

})

translatorDf <- data.frame(from = letters[1:7], to = rep(c("A","B"), c(3,4)))
translatorDf[7, ] <- translatorDf[3, ]

test_that("Translator df validator", {

  expect_message(validateTranslatorDf(translatorDf, "from", "to"))

})
