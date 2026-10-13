# Run against the genuine installed CRAN ESG 1.3 package; not a formula reimplementation.
# Rscript reference/generate_reference.R reference/expected.tsv
suppressPackageStartupMessages(library(ESG))
stopifnot(as.character(packageVersion('ESG')) == '1.3')
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else 'reference/expected.tsv'
ns <- asNamespace('ESG')
env <- parent.env(ns)
old_normal <- get('rnorm', env)
unlockBinding('rnorm', env)
queue <- numeric()
cursor <- 0L
assign('rnorm', function(n, mean=0, sd=1) {
  ix <- cursor + seq_len(n)
  stopifnot(max(ix) <= length(queue))
  cursor <<- cursor + n
  mean + sd * queue[ix]
}, env)
lockBinding('rnorm', env)
reset <- function(x) { queue <<- as.numeric(x); cursor <<- 0L }
n <- 3L; horizon <- 5L
zero <- 0.018 + 0.00003 * seq_len(85) + sin(seq_len(85) * 0.2) * 0.0002
e <- array(sin(seq_len(n*horizon*5) * 0.7) * 0.7, c(n,horizon,5))
obj <- new('Scenarios')
obj <- setParamsBaseScenarios(obj,horizon,n)
obj <- setRiskParamsScenarios(obj, vol=.017,k=.42,volStock=.23,volRealEstate=.12,volDefault=.018,
 alpha=.6,beta=.05,eta=.005,rho=-.35,stock0=1.1,realEstate0=1.4,liquiditySpread0=.011,defaultSpread0=.035)
obj <- setForwardRates(obj,zero,horizon)
obj <- setZCRates(obj,zero,horizon)
rows <- list()
record <- function(name, x) { stopifnot(all(is.finite(x))); rows[[name]] <<- as.numeric(x) }
reset(e)
obj <- customPathsGeneration(obj)
record('forward',getForwardRates(obj))
record('zero',getZCRates(obj))
record('rate',getShortRatePaths(obj))
record('stock',getstockPaths(obj)$stockPaths)
record('estate',getrealEstatePaths(obj)$realEstatePaths)
record('liquidity',getLiquiditySpreadPaths(obj))
record('default',getdefaultSpreadPaths(obj))
reset(e[,,1]); record('rShortRate',rShortRate(horizon,n,zero,.017,.42))
reset(c(e[,,1],e[,,1],e[,,2])); record('rStock',rStock(horizon,n,zero,.017,.42,.23,1.1,-.35)$stockPaths)
reset(c(e[,,1],e[,,1],e[,,3])); record('rRealEstate',rRealEstate(horizon,n,zero,.017,.42,.12,1.4)$realEstatePaths)
reset(e[,,4]); record('rLiquiditySpread',rLiquiditySpread(horizon,n,.005,.011))
reset(e[,,5]); record('rDefaultSpread',rDefaultSpread(horizon,n,.035,.018,.6,.05))
reset(e); allrisk <- rAllRisksFactors(horizon,n,zero,.017,.42,.23,1.1,-.35,.12,1.4,.005,.011,.035,.018,.6,.05)
record('rAllRisksFactors',unlist(allrisk))
# Plotting side effect redirected to a temporary, deleted file; numerical return is the fixture.
plotfile <- tempfile(fileext='.pdf'); pdf(plotfile)
record('martingale',MartingaleTest(obj)); dev.off(); unlink(plotfile)
kinds <- c('Zero-Coupon','Bond','CBond','ConvBond','EuroCall_Stock','EuroPut_Stock','EuroCall_ZC','EuroPut_ZC','CDS')
functions <- c('ZCBond_PriceDistribution','Bond_PriceDistribution','CBond_PriceDistribution','ConvBond_PriceDistribution',
 'EuroCall_Stock_PriceDistribution','EuroPut_Stock_PriceDistribution','EuroCall_ZC_PriceDistribution',
 'EuroPut_ZC_PriceDistribution','CDSPremium_PriceDistribution')
for (t in 0:2) for (j in seq_along(kinds)) {
 a <- list(.Object=obj,t=t,T=3.25)
 if (j %in% c(2,3,4)) a <- c(a,list(nCoupons=2,couponsRate=.025))
 if (j %in% c(3,9)) a <- c(a,list(omega=.4))
 if (j %in% 5:8) a <- c(a,list(Strike=if(j<7) 1.05 else .93))
 if (j %in% 7:8) a <- c(a,list(s=4.5))
 record(paste0('price_',t,'_',j),do.call(get(functions[j]),a))
 a$type <- kinds[j]
 record(paste0('dispatch_',t,'_',j),do.call(Asset_PriceDistribution,a))
}
# End-to-end exported R wrapper with its actual internally selected horizon (ceiling 3.25=4).
reset(e[,1:4,1]); record('rAssetDistribution',rAssetDistribution('Zero-Coupon',1,3.25,.017,.42,zero,nScenarios=n))
record('maturity_stock',EuroCall_Stock_PriceDistribution(obj,2,2,1.05))
record('maturity_cds',CDSPremium_PriceDistribution(obj,2,2,.4))
con <- file(out,'w')
for (key in names(rows)) writeLines(paste(key, seq_along(rows[[key]]), format(rows[[key]],digits=17,scientific=TRUE), sep='\t'),con)
close(con)
unlockBinding('rnorm',env); assign('rnorm',old_normal,env); lockBinding('rnorm',env)
cat('ESG',as.character(packageVersion('ESG')),'R',as.character(getRversion()),'fixtures',length(rows),'values',sum(lengths(rows)),'\n')
