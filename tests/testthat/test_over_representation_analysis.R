data("humanHallmarks")
f <- sample(humanHallmarks$gene_symbol, 500)
fList <- split(humanHallmarks, humanHallmarks$gs_name)
fList <- lapply(fList, function(x) as.character(x$gene_symbol))
res <- overRepresentationAnalysis(features = f, funCatList = fList)

test_that("Over-representation analysis", {

  expect_equal(class(res), "data.frame")

})

