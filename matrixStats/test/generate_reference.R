library(matrixStats)
stopifnot(as.character(packageVersion("matrixStats")) == "1.5.0")
x <- matrix(c(1,4,-2,9,2,4,-1,8,3,8,0,7,4,8,1,6,5,9,2,5),4,5)
xn <- x; xn[1,2] <- NA; xn[3,4] <- NaN
v <- c(4,1,7,1,9,2); vn <- v; vn[2] <- NA
w <- c(1,2,0,3,4,2); wr <- c(1,0,2,3,1); wc <- c(1,2,0,3)
output <- commandArgs(trailingOnly=TRUE)[1]
con <- file(output, "w")
write_ref <- function(x) {
 x <- as.numeric(x)
 text <- vapply(x, function(v) if (is.na(v)) "NaN" else if (v == Inf) "Inf" else if (v == -Inf) "-Inf" else sprintf("%.17g",v), "")
 writeLines(paste(text,collapse=","), con)
}
# rowSums2
write_ref(rowSums2(x))
# rowSums2 selected
write_ref(rowSums2(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowSums2 missing
write_ref(rowSums2(xn))
# rowMeans2
write_ref(rowMeans2(x))
# rowMeans2 selected
write_ref(rowMeans2(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowMeans2 missing
write_ref(rowMeans2(xn))
# rowVars
write_ref(rowVars(x))
# rowVars selected
write_ref(rowVars(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowVars missing
write_ref(rowVars(xn))
# rowSds
write_ref(rowSds(x))
# rowSds selected
write_ref(rowSds(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowSds missing
write_ref(rowSds(xn))
# rowMedians
write_ref(rowMedians(x))
# rowMedians selected
write_ref(rowMedians(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowMedians missing
write_ref(rowMedians(xn))
# rowMads
write_ref(rowMads(x))
# rowMads selected
write_ref(rowMads(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowMads missing
write_ref(rowMads(xn))
# rowIQRs
write_ref(rowIQRs(x))
# rowIQRs selected
write_ref(rowIQRs(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowIQRs missing
write_ref(rowIQRs(xn))
# rowMins
write_ref(rowMins(x))
# rowMins selected
write_ref(rowMins(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowMins missing
write_ref(rowMins(xn))
# rowMaxs
write_ref(rowMaxs(x))
# rowMaxs selected
write_ref(rowMaxs(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowMaxs missing
write_ref(rowMaxs(xn))
# rowLogSumExps
write_ref(rowLogSumExps(x))
# rowLogSumExps selected
write_ref(rowLogSumExps(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowLogSumExps missing
write_ref(rowLogSumExps(xn))
# rowCounts
write_ref(rowCounts(x))
# rowCounts selected
write_ref(rowCounts(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowCounts missing
write_ref(rowCounts(xn))
# rowAlls
write_ref(rowAlls(x))
# rowAlls selected
write_ref(rowAlls(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowAlls missing
write_ref(rowAlls(xn))
# rowAnys
write_ref(rowAnys(x))
# rowAnys selected
write_ref(rowAnys(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowAnys missing
write_ref(rowAnys(xn))
# rowProds
write_ref(rowProds(x))
# rowProds selected
write_ref(rowProds(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# rowProds missing
write_ref(rowProds(xn))
# rowRanges
write_ref(rowRanges(x))
# rowRanges missing
write_ref(rowRanges(xn))
# rowQuantiles
write_ref(rowQuantiles(x,probs=c(0,.25,.5,.75,1)))
# rowRanks
write_ref(rowRanks(x,ties.method="average"))
# rowRanks missing
write_ref(rowRanks(xn,ties.method="average"))
# rowCumsums
write_ref(rowCumsums(x))
# rowCumsums missing
write_ref(rowCumsums(xn))
# rowCumprods
write_ref(rowCumprods(x))
# rowCumprods missing
write_ref(rowCumprods(xn))
# rowCummins
write_ref(rowCummins(x))
# rowCummins missing
write_ref(rowCummins(xn))
# rowCummaxs
write_ref(rowCummaxs(x))
# rowCummaxs missing
write_ref(rowCummaxs(xn))
# rowDiffs
write_ref(rowDiffs(x))
# rowDiffs missing
write_ref(rowDiffs(xn))
# rowOrderStats
write_ref(rowOrderStats(x,which=2))
# rowAnyMissings
write_ref(rowAnyMissings(xn))
# rowAnyNAs
write_ref(rowAnyNAs(xn))
# rowVarDiffs
write_ref(rowVarDiffs(x,diff=2))
# rowVarDiffs missing
write_ref(rowVarDiffs(xn,na.rm=TRUE,trim=.2))
# rowSdDiffs
write_ref(rowSdDiffs(x,diff=2))
# rowSdDiffs missing
write_ref(rowSdDiffs(xn,na.rm=TRUE,trim=.2))
# rowMadDiffs
write_ref(rowMadDiffs(x,diff=2))
# rowMadDiffs missing
write_ref(rowMadDiffs(xn,na.rm=TRUE,trim=.2))
# rowIQRDiffs
write_ref(rowIQRDiffs(x,diff=2))
# rowIQRDiffs missing
write_ref(rowIQRDiffs(xn,na.rm=TRUE,trim=.2))
# rowWeightedMeans no weights
write_ref(rowWeightedMeans(x))
# rowWeightedMeans
write_ref(rowWeightedMeans(x,w=wr))
# rowWeightedMeans missing
write_ref(rowWeightedMeans(xn,w=wr,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# rowWeightedMedians no weights
write_ref(rowWeightedMedians(x))
# rowWeightedMedians
write_ref(rowWeightedMedians(x,w=wr))
# rowWeightedMedians missing
write_ref(rowWeightedMedians(xn,w=wr,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# rowWeightedVars no weights
write_ref(rowWeightedVars(x))
# rowWeightedVars
write_ref(rowWeightedVars(x,w=wr))
# rowWeightedVars missing
write_ref(rowWeightedVars(xn,w=wr,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# rowWeightedSds no weights
write_ref(rowWeightedSds(x))
# rowWeightedSds
write_ref(rowWeightedSds(x,w=wr))
# rowWeightedSds missing
write_ref(rowWeightedSds(xn,w=wr,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# rowWeightedMads no weights
write_ref(rowWeightedMads(x))
# rowWeightedMads
write_ref(rowWeightedMads(x,w=wr))
# rowWeightedMads missing
write_ref(rowWeightedMads(xn,w=wr,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# rowTabulates
write_ref(rowTabulates(matrix(as.integer(xn),4,5),values=c(1,4,NA)))
# rowCollapse
write_ref(rowCollapse(x,idxs=c(3,1,2)))
# colSums2
write_ref(colSums2(x))
# colSums2 selected
write_ref(colSums2(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colSums2 missing
write_ref(colSums2(xn))
# colMeans2
write_ref(colMeans2(x))
# colMeans2 selected
write_ref(colMeans2(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colMeans2 missing
write_ref(colMeans2(xn))
# colVars
write_ref(colVars(x))
# colVars selected
write_ref(colVars(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colVars missing
write_ref(colVars(xn))
# colSds
write_ref(colSds(x))
# colSds selected
write_ref(colSds(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colSds missing
write_ref(colSds(xn))
# colMedians
write_ref(colMedians(x))
# colMedians selected
write_ref(colMedians(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colMedians missing
write_ref(colMedians(xn))
# colMads
write_ref(colMads(x))
# colMads selected
write_ref(colMads(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colMads missing
write_ref(colMads(xn))
# colIQRs
write_ref(colIQRs(x))
# colIQRs selected
write_ref(colIQRs(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colIQRs missing
write_ref(colIQRs(xn))
# colMins
write_ref(colMins(x))
# colMins selected
write_ref(colMins(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colMins missing
write_ref(colMins(xn))
# colMaxs
write_ref(colMaxs(x))
# colMaxs selected
write_ref(colMaxs(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colMaxs missing
write_ref(colMaxs(xn))
# colLogSumExps
write_ref(colLogSumExps(x))
# colLogSumExps selected
write_ref(colLogSumExps(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colLogSumExps missing
write_ref(colLogSumExps(xn))
# colCounts
write_ref(colCounts(x))
# colCounts selected
write_ref(colCounts(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colCounts missing
write_ref(colCounts(xn))
# colAlls
write_ref(colAlls(x))
# colAlls selected
write_ref(colAlls(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colAlls missing
write_ref(colAlls(xn))
# colAnys
write_ref(colAnys(x))
# colAnys selected
write_ref(colAnys(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colAnys missing
write_ref(colAnys(xn))
# colProds
write_ref(colProds(x))
# colProds selected
write_ref(colProds(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE))
# colProds missing
write_ref(colProds(xn))
# colRanges
write_ref(colRanges(x))
# colRanges missing
write_ref(colRanges(xn))
# colQuantiles
write_ref(colQuantiles(x,probs=c(0,.25,.5,.75,1)))
# colRanks
write_ref(colRanks(x,ties.method="average"))
# colRanks missing
write_ref(colRanks(xn,ties.method="average"))
# colCumsums
write_ref(colCumsums(x))
# colCumsums missing
write_ref(colCumsums(xn))
# colCumprods
write_ref(colCumprods(x))
# colCumprods missing
write_ref(colCumprods(xn))
# colCummins
write_ref(colCummins(x))
# colCummins missing
write_ref(colCummins(xn))
# colCummaxs
write_ref(colCummaxs(x))
# colCummaxs missing
write_ref(colCummaxs(xn))
# colDiffs
write_ref(colDiffs(x))
# colDiffs missing
write_ref(colDiffs(xn))
# colOrderStats
write_ref(colOrderStats(x,which=2))
# colAnyMissings
write_ref(colAnyMissings(xn))
# colAnyNAs
write_ref(colAnyNAs(xn))
# colVarDiffs
write_ref(colVarDiffs(x,diff=2))
# colVarDiffs missing
write_ref(colVarDiffs(xn,na.rm=TRUE,trim=.2))
# colSdDiffs
write_ref(colSdDiffs(x,diff=2))
# colSdDiffs missing
write_ref(colSdDiffs(xn,na.rm=TRUE,trim=.2))
# colMadDiffs
write_ref(colMadDiffs(x,diff=2))
# colMadDiffs missing
write_ref(colMadDiffs(xn,na.rm=TRUE,trim=.2))
# colIQRDiffs
write_ref(colIQRDiffs(x,diff=2))
# colIQRDiffs missing
write_ref(colIQRDiffs(xn,na.rm=TRUE,trim=.2))
# colWeightedMeans no weights
write_ref(colWeightedMeans(x))
# colWeightedMeans
write_ref(colWeightedMeans(x,w=wc))
# colWeightedMeans missing
write_ref(colWeightedMeans(xn,w=wc,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# colWeightedMedians no weights
write_ref(colWeightedMedians(x))
# colWeightedMedians
write_ref(colWeightedMedians(x,w=wc))
# colWeightedMedians missing
write_ref(colWeightedMedians(xn,w=wc,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# colWeightedVars no weights
write_ref(colWeightedVars(x))
# colWeightedVars
write_ref(colWeightedVars(x,w=wc))
# colWeightedVars missing
write_ref(colWeightedVars(xn,w=wc,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# colWeightedSds no weights
write_ref(colWeightedSds(x))
# colWeightedSds
write_ref(colWeightedSds(x,w=wc))
# colWeightedSds missing
write_ref(colWeightedSds(xn,w=wc,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# colWeightedMads no weights
write_ref(colWeightedMads(x))
# colWeightedMads
write_ref(colWeightedMads(x,w=wc))
# colWeightedMads missing
write_ref(colWeightedMads(xn,w=wc,na.rm=TRUE,rows=c(4,2),cols=c(5,2,1)))
# colTabulates
write_ref(colTabulates(matrix(as.integer(xn),4,5),values=c(1,4,NA)))
# colCollapse
write_ref(colCollapse(x,idxs=c(3,1,2)))
# sum2
write_ref(sum2(v))
# sum2 subset
write_ref(sum2(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# mean2
write_ref(mean2(v))
# mean2 subset
write_ref(mean2(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# product
write_ref(product(v))
# product subset
write_ref(product(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# allValue
write_ref(allValue(v))
# allValue subset
write_ref(allValue(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# anyValue
write_ref(anyValue(v))
# anyValue subset
write_ref(anyValue(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# count
write_ref(count(v))
# count subset
write_ref(count(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# logSumExp
write_ref(logSumExp(v))
# logSumExp subset
write_ref(logSumExp(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# iqr
write_ref(iqr(v))
# iqr subset
write_ref(iqr(vn,idxs=c(6,1,4,2),na.rm=TRUE))
# anyMissing
write_ref(anyMissing(vn))
# diff2
write_ref(diff2(v,lag=2,differences=2))
# varDiff
write_ref(varDiff(v,diff=2,trim=.2))
# sdDiff
write_ref(sdDiff(v,diff=2,trim=.2))
# madDiff
write_ref(madDiff(v,diff=2,trim=.2))
# iqrDiff
write_ref(iqrDiff(v,diff=2,trim=.2))
# weightedMean
write_ref(weightedMean(v,w))
# weightedMean missing subset
write_ref(weightedMean(vn,w,idxs=c(6,3,2,1),na.rm=TRUE))
# weightedMedian
write_ref(weightedMedian(v,w))
# weightedMedian missing subset
write_ref(weightedMedian(vn,w,idxs=c(6,3,2,1),na.rm=TRUE))
# weightedVar
write_ref(weightedVar(v,w))
# weightedVar missing subset
write_ref(weightedVar(vn,w,idxs=c(6,3,2,1),na.rm=TRUE))
# weightedSd
write_ref(weightedSd(v,w))
# weightedSd missing subset
write_ref(weightedSd(vn,w,idxs=c(6,3,2,1),na.rm=TRUE))
# weightedMad
write_ref(weightedMad(v,w))
# weightedMad missing subset
write_ref(weightedMad(vn,w,idxs=c(6,3,2,1),na.rm=TRUE))
# binCounts
write_ref(binCounts(v,bx=c(0,2,4,8,10)))
# binCounts right
write_ref(binCounts(v,bx=c(0,2,4,8,10),right=TRUE))
# binMeans
write_ref(binMeans(vn,v,bx=c(0,2,4,8,10)))
# signTabulate
write_ref(signTabulate(c(-Inf,-2,0,NA,NaN,Inf,1)))
# indexByRow
write_ref(indexByRow(c(4,5),idxs=c(1,2,5,6,20)))
# rowAvgsPerColSet
write_ref(rowAvgsPerColSet(x,S=matrix(c(1,3,5,2,4,NA),3,2)))
# colAvgsPerRowSet
write_ref(colAvgsPerRowSet(x,S=matrix(c(1,3,4,2,4,NA),3,2)))
# x_OP_y+
write_ref(x_OP_y(xn,c(2,NA,4),OP="+",commute=TRUE,na.rm=TRUE))
# x_OP_y-
write_ref(x_OP_y(xn,c(2,NA,4),OP="-",commute=TRUE,na.rm=TRUE))
# x_OP_y*
write_ref(x_OP_y(xn,c(2,NA,4),OP="*",commute=TRUE,na.rm=TRUE))
# x_OP_y/
write_ref(x_OP_y(xn,c(2,NA,4),OP="/",commute=TRUE,na.rm=TRUE))
# t_tx_OP_y+
write_ref(t_tx_OP_y(xn,c(2,NA,4),OP="+",commute=TRUE,na.rm=TRUE))
# t_tx_OP_y-
write_ref(t_tx_OP_y(xn,c(2,NA,4),OP="-",commute=TRUE,na.rm=TRUE))
# t_tx_OP_y*
write_ref(t_tx_OP_y(xn,c(2,NA,4),OP="*",commute=TRUE,na.rm=TRUE))
# t_tx_OP_y/
write_ref(t_tx_OP_y(xn,c(2,NA,4),OP="/",commute=TRUE,na.rm=TRUE))
# rowMeans centered
write_ref(rowVars(x,center=rowMeans(x)))
# rowMAD centered
write_ref(rowMads(x,center=rowMedians(x),constant=2))
# rownegative subset
write_ref(rowSums2(x,rows=c(-1,-3),cols=c(0,2,5)))
# colMeans centered
write_ref(colVars(x,center=colMeans(x)))
# colMAD centered
write_ref(colMads(x,center=colMedians(x),constant=2))
# colnegative subset
write_ref(colSums2(x,rows=c(-1,-3),cols=c(0,2,5)))
# sum2 empty
write_ref(sum2(numeric()))
# sum2 infinities
write_ref(sum2(c(-Inf,Inf)))
# mean2 empty
write_ref(mean2(numeric()))
# mean2 infinities
write_ref(mean2(c(-Inf,Inf)))
# product empty
write_ref(product(numeric()))
# product infinities
write_ref(product(c(-Inf,Inf)))
# logSumExp empty
write_ref(logSumExp(numeric()))
# logSumExp infinities
write_ref(logSumExp(c(-Inf,Inf)))
# iqr empty
write_ref(iqr(numeric()))
# iqr infinities
write_ref(iqr(c(-Inf,Inf)))
# all missing decisive
write_ref(allValue(c(NA,0),value=1))
# any missing decisive
write_ref(anyValue(c(NA,1),value=1))
# count NA
write_ref(count(c(NA,1,NaN),value=NA_real_))
# colRanks preserveShape
write_ref(colRanks(x,ties.method="dense",preserveShape=TRUE))
close(con)
