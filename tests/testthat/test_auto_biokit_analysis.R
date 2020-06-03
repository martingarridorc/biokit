data("sarsCovData")
data("humanHallmarks")

# prepare testing matrix
testMat <- sarsCovMat[1:1000,]
# remove genes not expressed
notExpr <- apply(testMat, 1, function(x) any(x == 0))
testMat <- testMat[!notExpr,]
# subset to genes with high sd to avoid ties
tmm <- countsToTmm(testMat)

res1 <- autoBiokitAnalysis(mat = testMat, sampInfo = sarsCovSampInfo, groupCol = "group", funCatList = humanHallmarks, statMode = "edgeR", annotationMode = "byCutoff", filterByExpr = TRUE)
res2 <- autoBiokitAnalysis(mat = countsToTmm(testMat), sampInfo = sarsCovSampInfo, groupCol = "group", funCatList = humanHallmarks, statMode = "limma", annotationMode = "byRank", rankCol = "logFc", topN = 50)
res3 <- autoBiokitAnalysis(mat = countsToTmm(testMat), sampInfo = sarsCovSampInfo, groupCol = "group", funCatList = humanHallmarks, statMode = "tTest", annotationMode = "byCutoff")



