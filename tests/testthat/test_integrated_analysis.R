data("sarsCovData")
data("humanHallmarks")

oraRes <- overRepresentationAnalysis(features = rownames(sarsCovMat)[1:200], funCatList = humanHallmarks)
sarsCovMat <- sarsCovMat[rowSums(sarsCovMat) >= 15, ]
tmmMat <- countsToTmm(sarsCovMat)
violinPlot(tmmMat)
pcaPlot(mat = tmmMat, sampInfo = sarsCovSampInfo, groupCol = "group")
heatmapPlot(mat = tmmMat, sampInfo = sarsCovSampInfo, groupCol = "group",scaleBy = "row",  nTop = 25)
diffRes <- biokit::autoLimmaComparison(mat = tmmMat, sampInfo = sarsCovSampInfo, groupCol = "group")
diffEdgeRes <- autoEdgerComparison(counts = sarsCovMat, sampInfo = sarsCovSampInfo, groupCol = "group")
diffGeneric <- autoPairwiseMatrixTest(mat = abs(tmmMat)[1:20,], sampInfo = sarsCovSampInfo, groupCol = "group")
statusRes <- addStatusByCutoff(resDf = diffRes, splitUpDown = TRUE)
statusRes <- addStatusByRank(resDf = diffRes, rankCol = "pAdj", topN = 10, splitUpDown = TRUE)
volcanoPlot(diffRes)
gseaResults <- gseaFromStats(df = diffRes, funCatList = humanHallmarks, rankCol = "logFc", splitCol = "comparison")
gseaPlot(gseaResults)

