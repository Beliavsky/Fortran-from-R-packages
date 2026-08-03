require(CLA)

b64 <- .Machine$sizeof.pointer == 8
cat(sprintf("%d bit platform type '%s'\n", if(b64) 64 else 32, .Platform$OS.type))
(nonWindows <- .Platform$OS.type != "windows")
arch <- Sys.info()[["machine"]]
.M <- .Machine; str(.M[grep("^sizeof", names(.M))]) ## differentiate long-double..
## Do we have 64bit but no-long-double ?
b64nLD <- (arch == "x86_64" && .M$sizeof.longdouble != 16)
if(b64nLD) arch <- paste0(arch, "--no-long-double")
arch

data(muS.sp500)## <<-- working from this data ---------------------

##-->  ../man/findMu.Rd      -- has this example --
set.seed(2016)
iS <- sample.int(length(muS.sp500$mu), 17)
cov17 <- muS.sp500$covar[iS, iS]
CLsp.17 <- CLA(muS.sp500$mu[iS], covar=cov17, lB=0, uB = 1/2)
CLsp.17 # 16 turning points
tpS <- CLsp.17$MS_weights[,"Sig"]
str(s0 <- c(min(tpS), seq(0.0188, 0.039, by = 0.0001), max(tpS)))
mu.. <- findMu(s0, result=CLsp.17, covar=cov17)
str(mu..)
stopifnot(dim(mu..$weight) == c(17, length(s0)))
## and then plots


##-->  ../man/findSig.Rd     -- has this example
set.seed(2018)
iS <- sample.int(length(muS.sp500$mu), 21)
cov21 <- muS.sp500$covar[iS, iS]
CLsp.21 <- CLA(muS.sp500$mu[iS], covar=cov21, lB=0, uB = 1/2)
CLsp.21 # had 14, now (R-devel new sample()) 21 turning points
tpM <- CLsp.21$MS_weights[,"Mu"]
str(m0 <- c(min(tpM), seq(0.0028, 0.0052, by = 0.00005), max(tpM)))
sig. <- findSig(m0, result=CLsp.21, covar=cov21)
str(sig.)
stopifnot(dim(sig.$weight) == c(21, length(m0)))
## and then plots
