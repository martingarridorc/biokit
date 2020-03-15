#### CBD DATA
cbdCountsUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131565/suppl/GSE131565_raw_counts.txt.gz"
# download file and read table
tempFile <- tempfile()
download.file(cbdCountsUrl, tempFile)
cbdCounts <- read.table(file = tempFile, sep ="\t", header = TRUE, row.names = 1)
cbdCounts <- as.matrix(cbdCounts[, 6:11])
# subset to random 1000 genes
set.seed(149)
selGenes <- sample(rownames(cbdCounts), size = 1000)
cbdCounts <- cbdCounts[selGenes, ]
# create sample Info
cbdSampInfo <- data.frame(row.names = colnames(cbdCounts), group = rep(c("control","cbd"), each = 3))
# save results
save(cbdCounts, cbdSampInfo, file = "data/cbdSubset.rda")

#### BLM DATA
blmCountsUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121659/suppl/GSE121659_raw_counts_merged.txt.gz"
# download file and read table
tempFile <- tempfile()
download.file(blmCountsUrl, tempFile)
blmCounts <- read.table(file = tempFile, sep ="\t", header = TRUE, row.names = 1)
blmCounts <- blmCounts[, 6:17]
# subset to random 1000 genes
set.seed(149)
selGenes <- sample(rownames(blmCounts), size = 1000)
blmCounts <- blmCounts[selGenes, ]
# create sample Info
blmSampinfo <- data.frame(row.names = colnames(blmCounts), group = rep(c("control","blm", "blm_ehp", "blm_aja"), each = 3))
# save results
save(blmCounts, blmSampinfo, file = "data/blmSubset.rda")

