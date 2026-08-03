require(CLA)

source(system.file("xtraR", "platform-sessionInfo.R", # <<- ../inst/xtraR/platform-sessionInfo.R
		   package = "CLA", mustWork=TRUE)) # moreSessionInfo(), withAutoprint() ..
mS <- moreSessionInfo(print. = TRUE)

data(muS.sp500)# not 500, but 476 assets <--> ../man/muS.sp500.Rd <<<

if(requireNamespace("FRAPO")) {
    data(SP500, package = "FRAPO")
    print(dim(SP500)) # 265 * 476
}

(n <- length(muS.sp500$mu)) # 476
system.time(# ~ 9 sec (64-bit); 13.8 sec (32-b florence); seen 27.44 sec on Winb.32
    CLs5c.0.120 <- CLA(muS.sp500$mu, muS.sp500$covar, lB=0, uB=1/20)
)
CLs5c.0.120 # -> print() method
ws <- CLs5c.0.120$weights_set
table(rowSums(ws) > 0)
## FALSE  TRUE
##   397    79  == we check below that indeed  79 assets have some non-0 wts


uncall <- function(x) `$<-`(x, call, NULL)

doExtras <- TRUE  # for experiments, not normally
doExtras <- FALSE

if(doExtras) {
  system.time({
    tols <- 10^-c(1,3,5:9,11,14)
    names(tols) <- paste0("10^", round(log10(tols)))
    CLs5c.ls <- lapply(tols, function(tol)
        CLA(muS.sp500$mu, muS.sp500$covar, lB=0, uB=1/20, tol.lambda = tol))
  }) #  78.101 elapsed [nb-mm4] ; 46.108 [lynne 2018-10]

    identical(uncall(CLs5c.ls[["10^-7"]]), uncall(CLs5c.0.120))
    for(i in seq_along(tols)[-1]) {
        cat("--=--=--=--=--\n", (n1 <- names(tols[i-1])), " vs. ", (n2 <- names(tols[i])), ": ")
        ae <- all.equal(uncall(CLs5c.ls[[i-1]]),
                        uncall(CLs5c.ls[[ i ]]))
        if(isTRUE(ae)) cat(" are all.equal()\n")
        else {
            CLA.i. <- uncall(CLs5c.ls[[i-1]]) ; wgt <- CLA.i.$weights_set
            cat("are different [all.equal()]: dim(..[[",n1,"]]$weights_set) =",
                dim(wgt)[1],"x", dim(wgt)[2],"\n")
        }
    }
}
## 2018-10 lynne, 64b Fedora 28
##  10^-1  vs.  10^-3 : are different [all.equal()]: dim(..[[ 10^-1 ]]$weights_set) = 476 x 47
##  10^-3  vs.  10^-5 : are different [all.equal()]: dim(..[[ 10^-3 ]]$weights_set) = 476 x 156
##  10^-5  vs.  10^-6 :  are all.equal()
##   ................     "   "    "

op <- options(width = max(500, getOption("width"))) # then it actually fits

## Mostly keep the 79 assets with non-zero weights:
## --->  wts.non0 { 79 x 161 } matrix
if(require(Matrix)) withAutoprint(local = FALSE, {
    ## visualize how weights change "along turning points"
    spWts <- Matrix(CLs5c.0.120$weights_set, sparse=TRUE)
    image(spWts, xlab = "turning point", ylab = "asset number")
    ##
    dim(spWts.non0 <- spWts[rowSums(spWts) > 0 , ])
    round(1000 * spWts.non0) ##-> e.g. ./Solaris_wts_non0.txt ==> readWn0()
    ##    ~~~~~~~~~~~~~~~~~
    image(spWts.non0, xlab = "turning point", ylab = "asset number")
    wts.non0 <- as(spWts.non0, "matrix")
}) else {
    warning("'Matrix' package not available -- should not happen!")
    wts.non0 <- CLs5c.0.120$weights_set[rowSums(CLs5c.0.120$weights_set) > 0 , ]
    if(is.null(colnames(wts.non0))) ## empty column names for nice printing:
        colnames(wts.non0) <- rep("", ncol(wts.non0))
    print.table(round(1000 * wts.non0), zero.print=".")
}

##' To read the output from the above round(1000*spWts.non) into a (non-sparse) matrix:
readWn0 <- function(txtfile) {
    stopifnot(file.exists(txtfile))
    lns <- strsplit((ln <- readLines(txtfile)), split = " +")
    n <- unique(lengths(lns))
    if(length(n) != 1)
        stop("number of 'words' in each line differ:", paste(n, collapse=", "))
    nms <- lapply(lns, `[[`, 1)
    rnls <- lapply(lns, `[`, -1L) # the non-names, i.e. round(1000 * wts[])
    w1000 <- setNames(rnls, nms)
    nr_non0 <- vapply(w1000, function(ch) sum(ch != "."), 1L)
    list(nr_non0=nr_non0, w1000=w1000)
}
if(FALSE) { # private file
    rr <- readWn0(txtfile = "./SP5c_Solaris_wts_non0.txt")
    rr <- readWn0(txtfile = "./SP5c_Win32_wts_non0.txt")
    rr <- readWn0(txtfile = "./SP5c_M1mac_wts_non0.txt") # this is really special
    nr_non0 <- rr$nr_non0
    ## compare with b64.n0, to construct  'S10.b32.n0' :
    ## same length (== 79)  *and* same names :
    stopifnot(identical(names(nr_non0), names(b64.n0)))
    d <- nr_non0 - b64.n0
    d[d != 0] # Solaris: the first, ADSK is 2, all others are 1 ..
    all((d1 <- d[d != 0][-1]) == 1) # TRUE  (FALSE for M1mac)
    dput(names(d1)) # the names for the '1' ..
    ## Win 32: 1 2 -1 1  ... 1 2 1 1 -- win 32 is really "big" differing:
    ##     table(d) .. only 23 do *not* differ
    ## -1  0  1  2
    ##  5 23 49  2     majority (49 of 79): "+1"
}

table(w.pos <- wts.non0 > 0)
## FALSE  TRUE
##  7808  4911
summary(log10(c(wts.non0[w.pos])))
## Lnx F34:
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## -18.363  -1.816  -1.405  -1.688  -1.301  -1.301

plot(ecdf(l10w <- log10(c(wts.non0[w.pos]))))
abline(v = log10(1/20), col=2, lwd=2, lty=2) # weights upper boundary = 1/20 (here)
rug(l10w) # wow, that's revealing!
mtext(paste(sfsmisc::shortRversion(), osVersion, sep=" -- "))
table(l10w < -10)
## FALSE  TRUE
##  4889    22 <-- get 23 (and different ones), for Lnx(64b) --no-long-double
sum(l10w < -5) # the *same* 22 | 23
round(sort(l10w[l10w < -10]), 1)
## -18.4 ..... Lnx&Win64/32: -15.4 -15.1 -15.1 -14.8 // Lnx_noLD: -15.8 -15.4 -14.8 // M1: .. -16.0 -15.4 -15.4
smll.e_10 <- 0 < wts.non0 & wts.non0 < 1e-10
iSml <- which(smll.e_10, arr.ind = TRUE)
iSml <- iSml[order(iSml[,1], iSml[,2]),]
data.frame(iSml, wt = wts.non0[smll.e_10])
##       row col           wt
## ADSK    2   7 8.881784e-16
## AET     3   4 1.776357e-15
## AMGN    4 159 4.440892e-16
## BCR     9  98 2.220446e-16
## BDX    10 126 8.881784e-16
## BUD    15 155 1.110223e-16
## CI     18  32 5.551115e-17
## CL     19  79 1.110223e-16
## CTL    24 157 1.387779e-17
## ED     29 130 1.387779e-17
## EIX    30  15 5.551115e-17
## FE     33  50 3.469447e-18
## HAL    35 131 1.127570e-17
## MCD    46  73 3.469447e-18
## NEM    53 117 1.734723e-18
## RAI    62  52 6.938894e-18
## SII    66 119 2.775558e-17
## SYK    68  63 6.938894e-18
## WMT    77 156 3.469447e-18
## XTO    79   3 2.602085e-18
## XTO.1  79  12 4.336809e-19
## XTO.2  79  96 8.673617e-19

## NEW:: set all really small wts to zero ==> I expect everything will be platform independent!
wts.non0[smll.e_10] <- 0 # !!
##                    ===

## from the 4911, now have 22 less:
stopifnot(sum(w.pos <- wts.non0 > 0) == 4889)

options(op)
stopifnot(nrow(wts.non0) == 79)
if(FALSE) # once, manually (into tests/ directory
    saveRDS(wts.non0, "wtsn0.rds")
file.info("wtsn0.rds")$size  # 26926
wtsn0.ref <- readRDS("wtsn0.rds")

## see on all platforms what we get -- on OpenBLAS, the dim() differs !
all.equal(target = wtsn0.ref, current = wts.non0, tol=0)
                                        # expect TRUE only on 64bit (Lnx, R's BLAS)
                                        # 3.10416e-15 and 1.366427e-15 on other BLAS
differWts <- ncol(wtsn0.ref) != ncol(wts.non0)
if(differWts) {
    cat("Got",ncol(wts.non0), "weights from CLA() -- different than ref with",
        ncol(wtsn0.ref), "\n")
    strict <- FALSE # !
} else {
    strict <- mS$strictR
    stopifnot(all.equal(target = wtsn0.ref, current = wts.non0,
                        tol = 1e-13))
}
non.0.assets <- Filter(function(.) . > 0, apply(wts.non0, 1, function(c) sum(c > 0)))

## New (after zapping 22 very small weights):

b64.n0 <-
    c(AAPL = 135L, ADSK = 65L, AET = 146L, AMGN = 2L, ATI = 75L,
      AYE = 56L, AZO = 26L, BAX = 95L, BCR = 34L, BDX = 35L, BIIB = 118L,
      BNI = 86L, BRL = 23L, BTU = 27L, BUD = 6L, CCE = 54L, CELG = 128L,
      CI = 68L, CL = 82L, CLX = 53L, CME = 140L, CNX = 16L, COST = 40L,
      CTL = 4L, CVS = 102L, DF = 36L, DGX = 33L, DVN = 14L, ED = 31L,
      EIX = 126L, ESRX = 48L, FCX = 54L, FE = 60L, GILD = 38L, HAL = 30L,
      HES = 40L, HST = 108L, HUM = 71L, INTU = 48L, JNJ = 34L, K = 61L,
      LH = 80L, LLL = 96L, LMT = 83L, LUK = 72L, MCD = 60L, MDT = 43L,
      MMC = 7L, MON = 53L, MRO = 136L, MTW = 66L, MUR = 97L, NEM = 44L,
      NOC = 74L, NUE = 30L, NVDA = 13L, PBG = 72L, PCP = 102L, PDCO = 71L,
      PEP = 69L, PG = 87L, RAI = 109L, RIG = 121L, RRC = 105L, RTN = 90L,
      SII = 26L, SSP = 14L, SYK = 18L, SYMC = 13L, TEX = 36L, TIE = 84L,
      TSO = 115L, TYC = 59L, UST = 127L, WAG = 17L, WFR = 5L, WMT = 5L,
      X = 43L, XTO = 99L)


non.0.TARG <- if(mS$ b64) {
                  if(mS$ b64nLD)
                      if(mS$arch == "arm64" && grepl("darwin", R.version$os)) # since Feb.2021 also for 'M1mac'
                          b64.n0 # M1.b64.n0
                      else
                          b64.n0 # b64nLD.n0
                  else
                      b64.n0
              } else { # 32 bit
                  if(mS$ osVersion == "Solaris 10")
                      b64.n0 # S10.b32.n0
                  else if(.Platform$OS.type == "windows")
                      b64.n0 # win.b32.n0
                  else # notably 32-bit Linux
                      b64.n0 # b32.n0
              }

## see on all platforms what we get;  typically no diff on 64bit *and* using R's BLAS/Lapack
if(all(non.0.assets == non.0.TARG)) { ## show differences:
    cat("Asset results == non.0.TARG\n")
} else { ## (not seen anymore, since "zap < 1e-10"):
    cat("\n'non.0.assets' differing from non.0.TARG:\n")
    cat("+1:\n"); dput(names(b64.n0)[non.0.assets - non.0.TARG == +1])
    cat("-1:\n"); dput(names(b64.n0)[non.0.assets - non.0.TARG == -1])

    if(any(isB <- abs(non.0.assets - non.0.TARG) > 1)) {
        cat("more different (than just +/- 1), showing differences:\n")
        dput((non.0.assets - non.0.TARG)[isB])
    }
}

## They have the same names and only differ by  +/- 1:
stopifnot(exprs = {
    !strict || identical(non.0.assets, non.0.TARG)

    differWts || identical(head(CLs5c.0.120$free_indices, 12),
              list(c(295L, 453L), 453L, c(453L, 472L), c(19L, 453L, 472L),
                   c(19L, 453L), 453L, c(15L, 453L), 15L, c(15L, 320L),
                   c(15L, 105L, 320L), c(105L, 320L), c(105L, 320L, 472L)))
})

## Check some of the 'Env<n>' versions: ---------

##' Transform CLA() result to old style (= Env8 / Env9 results):
claStrip <- function(res) {
    class(res) <- NULL
    res$call <- NULL
    res$emptyLambdas <- NULL
    res$weights_set <- unname(res$weights_set)
    names(res)[[match("MS_weights", names(res))]] <- "MS_weight" # "MS_weights" w/ final "s"
    res
}

rCLA <- claStrip(CLs5c.0.120)

##' Drop "first turning point" from old, pre-0.95,  CLA() result:
claDrop1st <- function(res) {
    res$weights_set <- res$weights_set[, -1L , drop=FALSE] # drop 1st column
    if(is.matrix(res$MS_weight))
        res$MS_weight <- res$MS_weight[ -1L, , drop=FALSE] # drop 1st row
    for(nm in c("free_indices", "gammas", "lambdas"))
        res[[nm]] <- res[[nm]][-1L]
    res
}

## back compatibility to "old" Env8 results:
nsCLA <- asNamespace("CLA")
if(is.environment(e8 <- nsCLA$Env8)) local(withAutoprint({
    system.time(r8 <- e8$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne (2017): 9.6--9.8 sec; 2018: 6.1 sec
    if(ncol(claDrop1st(r8)$weights_set) == ncol(rCLA$weights_set))
        stopifnot(all.equal(claDrop1st(r8), rCLA, tol = 1e-14)) # they are the same!
    else cat("#{columns} differ in r8\n")
}))

if(is.environment(e9 <- nsCLA$Env9)) local(withAutoprint({
    system.time(r9 <- e9$cla.solve(muS.sp500$mu, muS.sp500$covar,
                                   lB = rep(0,n), uB= rep(1/20, n)))
    ## lynne(2017): 10.0 sec;  2018: 6.6 sec
    if(ncol(claDrop1st(r9)$weights_set) == ncol(rCLA$weights_set))
        stopifnot(all.equal(claDrop1st(r9), rCLA, tol = 1e-14)) # they are the same!
    else cat("#{columns} differ in r9\n")
}))
