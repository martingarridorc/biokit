data("sarsCovData")

# prepare testing matrix
testMat <- sarsCovMat[1:20,]

# prepare matrix without rownames and with NAs
noRowNames <- testMat
rownames(noRowNames) <- NULL
naMat <- testMat
naMat[1:3,] <- NA

# prepare sample information and a column with wrong groups
sampInfo <- sarsCovSampInfo
sampInfo[ , "wrongGroups"] <- rep(c(" -A", "12", "heyIAmAlone", NA),  c(5,4,1,2))

# prepare translator df
translatorDf <- data.frame(from = letters[1:8], to = letters[9:16])
wrongTranslatorDf <- data.frame(from = letters[1:8], to = c(letters[9:14], NA, NA))
wrongTranslatorDf[9,] <- wrongTranslatorDf[1,]

test_that("Validators work nice", {

  expect_true(validateMatrix(testMat))
  expect_error(validateMatrix(data.frame()))
  expect_error(validateMatrix(noRowNames))
  expect_message(validateMatrix(naMat, allowNa = TRUE))
  expect_error(validateMatrix(naMat, allowNa = FALSE))

  expect_equal(class(validateSampInfo(sampInfo, "group")), "data.frame")
  expect_error(validateSampInfo(sampInfo, "group", testMat[,1:3]))
  expect_message(validateSampInfo(sampInfo, "wrongGroups"))

  expect_equal(class(validateTranslatorDf(translatorDf, "from", "to")), "data.frame")
  expect_equal(class(validateTranslatorDf(wrongTranslatorDf, "from", "to")), "data.frame")

  expect_message(createPairwiseContrasts(c(NA, NA, 'a', 'b', 'b')))

})

