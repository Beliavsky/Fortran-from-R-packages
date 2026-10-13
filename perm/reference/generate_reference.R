# SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
# Run with the official perm 1.0-0.4 package installed.
library(perm)
stopifnot(as.character(packageVersion("perm")) == "1.0.0.4")
args <- commandArgs(trailingOnly=TRUE)
path <- if(length(args)) args[1] else "reference/r_oracle.txt"
con <- file(path,"wt")
emit <- function(name,x) {
  cat(name, paste(sprintf("%.17g",as.numeric(x)),collapse=" "),"\n",file=con)
}
tails <- function(x) x$p.values[c("p.twosided","p.twosidedAbs","p.lte","p.gte","p.equal")]
s <- c(-2.3,-2.3,0,1.7,2.1,5.2)
g <- c(1,1,0,1,0,0)
emit("weighted",tails(calcPvals(c(-4,0,0,1,7),1,12,wgts=c(1,2,3,5,1))))
emit("ce",tails(twosample.exact.ce(s,g)))
emit("network",tails(twosample.exact.network(s,g)))
for(a in c("less","greater","two.sided","two.sidedAbs")) {
  r <- calcPvalsMC(c(1,1,3,6,8,12),5,12,a,p.conf.level=.95)
  emit(paste0("mc_",a),c(tails(r),r$p.conf.int))
}
r <- calcPvalsMC(c(1,1,3,6,8,12),5,12,twosidedTstat=TRUE,p.conf.level=.95)
emit("mc_upper",c(r$p.values,r$p.conf.int))
for(level in c(.95,.99,0,1)) for(hits in c(0,1,5,9,10)) emit("ci",pCI(hits,10,level))
for(hits in c(0,1,5000,9999)) emit("ci_large",pCI(hits,9999,.99))
r <- twosample.pclt(s,g)
emit("ts_pclt",c(r$p.values[c("p.twosided","p.twosidedAbs","p.lte","p.gte")],r$Z))
y <- c(2,8,4,3,9,5)
r <- trend.pclt(s,y)
emit("trend_pclt",c(r$p.values,r$Z))
gk <- c(1,1,2,2,3,3)
r <- ksample.pclt(s,gk)
emit("ks_pclt",c(r$p.values,r$chisq.value,r$df))
s2 <- rep(c(-.25,.5,1.25,2.75),each=10)
g2 <- rep(c(1,0),20)
emit("network_ties",tails(twosample.exact.network(s2,g2)))
emit("choose",as.vector(t(chooseMatrix(5,4))))
emit("getcnt",getcnt(c(1,2,1),c(1,3,4),c(2,1)))
close(con)
