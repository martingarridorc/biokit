# prepare cbd data for testing
data("blmData")
blmMat <- blmMat[rowSums(blmMat) > 25, ]
blmMat <- blmMat[1:50, ]
desMat <- designFromSampInfo(blmSampInfo, "group")
conMat <- contrastsFromDesign(desMat)
testDge <- edgeR::DGEList(counts = blmMat, group = blmSampInfo$group)
testSe <- SummarizedExperiment::SummarizedExperiment(blmMat, colData = blmSampInfo)

test_that("TMM normalization", {
    
    tmmRes <- countsToTmm(x = blmMat)
    expect_equal(class(tmmRes), "matrix")
    expect_type(tmmRes, "double")
    
})

test_that("edgeR Df from contrast and test design", {
    
    result <- edgeRFromContrasts(object = testDge, desMat = desMat, conMat = conMat)
    expect_equal(class(result), "data.frame")
    expect_equal(nrow(result), 300)
    
})

test_that("AutoLimma", {
    
    result <- autoEdgeR(se = testSe, "group")
    resultNoFilter <- autoEdgeR(se = testSe, "group", minCounts = 0, useFilterByExpr = FALSE)
    resultFilterByCounts <- autoEdgeR(se = testSe, "group", minCounts = 30, useFilterByExpr = FALSE)
    
    expect_equal(class(result), "data.frame")
    expect_equal(nrow(result), 276)
    expect_equal(nrow(resultNoFilter), 300)
    expect_equal(nrow(resultFilterByCounts), 294)
    
})

