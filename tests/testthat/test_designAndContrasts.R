data("blmSubset")

test_that("Create pairwise contrasts and contrasts from design", {

  d <- designFromSampInfo(blmSampInfo, column = "group")
  c <- contrastsFromDesign(x = d)
  p <- pairwiseContrasts(x = levels(blmSampInfo$group))

  # contrasts with combn
  t <- ncol(combn(levels(blmSampInfo$group), m = 2))

  # equal number of contrasts
  expect_equal(ncol(c), t)
  expect_equal(length(p), t)

})
