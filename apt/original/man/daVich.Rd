\name{daVich}
\alias{daVich}
\docType{data}
\title{Import prices and values of wooden beds from Vietnam and China}
\description{This data set contains two unit import prices (dollar per piece) and values (million dollars) of wooden beds from Vietnam and China to the United States.
  \tabular{ll}{
  \code{price.vi } \tab Monthly price over Janurary 2002 to Janurary 2010 from Vietnam. \cr
  \code{price.ch } \tab Monthly price over Janurary 2002 to Janurary 2010 from China. \cr
  \code{price.vi } \tab Monthly value over Janurary 2002 to Janurary 2010 from Vietnam. \cr
  \code{price.ch } \tab Monthly value over Janurary 2002 to Janurary 2010 from China. \cr
  }
}
\usage{data(daVich)}

\format{A monthly time series from Janurary 2002 to Janurary 2010 with 97 observations for each of the four series.}

\details{
Under the Harmonized Tariff Schedule (HTS) system, the commodity of wooden beds is classified as HTS 9403.50.9040. The monthly cost-insurance-freight values in dollar and quantities in piece are reported by country from U.S. ITC (2010). The unit price (dollar per piece) is calculated as the ratio of value over quantity by country.
}
\source{
U.S. ITC, 2010. Interactive tariff and trade data web. US International Trade Commision (Assecced on March 1, 2010).
}
\references{
Sun, C. 2011. Price dynamics in the import wooden bed market of the United States. Forest Policy and Economics 13(6): 479-487.
}
\examples{

# The following codes reproduce the main results in Sun (2011 FPE).
# All the codes have been tested and should work.

# 1. Data preparation __________________________________________________________

library(urca); data(daVich)
head(daVich); tail(daVich); str(daVich)
prVi <- y <- daVich[, 1] 
prCh <- x <- daVich[, 2]


# 2. EG cointegration  _________________________________________________________

LR <- lm(y ~ x); summary(LR)
(LR.coef <- round(summary(LR)$coefficients, 3))
(ry <-  ts(residuals(LR), start=start(prCh), end=end(prCh), frequency =12))
summary(eg <- ur.df(ry,  type=c("none"), lags=1)); plot(eg)
(eg4  <- Box.test(eg@res, lag = 4, type="Ljung") )
(eg8  <- Box.test(eg@res, lag = 8, type="Ljung") )
(eg12 <- Box.test(eg@res, lag = 12, type="Ljung"))


\donttest{

# 3. TAR + Cointegration  ______________________________________________________

# best threshold 
t3 <- ciTarThd(y=prVi, x=prCh, model="tar", lag=0) 
(th.tar <- t3$basic); plot(t3)
for (i in 1:12) {                 # 20 seconds
  t3a <- ciTarThd(y=prVi, x=prCh, model="tar", lag=i) 
  th.tar[i+2] <- t3a$basic[,2]
}
th.tar
 
t4 <- ciTarThd(y=prVi, x=prCh, model="mtar", lag=0); (th.mtar <- t4$basic)
plot(t4)
for (i in 1:12) {
  t4a <- ciTarThd(y=prVi, x=prCh, model="mtar", lag=i) 
  th.mtar[i+2] <- t4a$basic[,2]
 }
th.mtar

# The following threshold values are specific to this data. They HAVE TO be  
# revised for another data set. Otherwise, various errors will occur.
  t.tar <- -8.041; t.mtar <- -0.451     # lag = 0 to 4
# t.tar <- -8.701 ; t.mtar <- -0.451    # lag = 5 to 12

# lag selection
mx <- 12
(g1 <-ciTarLag(y=prVi, x=prCh, model="tar", maxlag=mx, thresh= 0));     plot(g1)
(g2 <-ciTarLag(y=prVi, x=prCh, model="mtar",maxlag=mx, thresh= 0));     plot(g2)
(g3 <-ciTarLag(y=prVi, x=prCh, model="tar", maxlag=mx, thresh=t.tar));  plot(g3)
(g4 <-ciTarLag(y=prVi, x=prCh, model="mtar",maxlag=mx, thresh=t.mtar)); plot(g4)


# Table 3 Results of EG and threshold cointegration tests
# Note: Some results in Table 3 in the published paper were incorrect because
#  of a mistake made when the paper was done in 2009. I found the mistake when 
#  the package was build in later 2010. For example, for the consistent MTAR, 
#  the coefficient for the positive term was reported as -0.251 (-2.130) but 
#  it should be -0.106 (-0.764), as cacluated from below codes. 
#  The main conclusion about the research issue is still qualitatively the same.
vv <- 3
(f1 <- ciTarFit(y=prVi, x=prCh, model="tar",  lag=vv, thresh=0     ))
(f2 <- ciTarFit(y=prVi, x=prCh, model="tar",  lag=vv, thresh=t.tar ))
(f3 <- ciTarFit(y=prVi, x=prCh, model="mtar", lag=vv, thresh=0     ))
(f4 <- ciTarFit(y=prVi, x=prCh, model="mtar", lag=vv, thresh=t.mtar))
 
r0 <- cbind(summary(f1)$dia, summary(f2)$dia, summary(f3)$dia, 
  summary(f4)$dia)
diag <- r0[c(1:4, 6:7, 12:14, 8, 9, 11), c(1,2,4,6,8)]
rownames(diag) <- 1:nrow(diag); diag

e1 <- summary(f1)$out; e2 <- summary(f2)$out
e3 <- summary(f3)$out; e4 <- summary(f4)$out; rbind(e1, e2, e3, e4)
ee <- list(e1, e2, e3, e4); vect <- NULL
for (i in 1:4) {
  ef <- data.frame(ee[i])
  vect2 <- c(paste(ef[3, "estimate"], ef[3, "sign"], sep=""), 
             paste("(", ef[3, "t.value"], ")",       sep=""),
             paste(ef[4, "estimate"], ef[4, "sign"], sep=""), 
             paste("(", ef[4, "t.value"], ")",       sep=""))    
  vect <- cbind(vect, vect2)
}
item <- c("pos.coeff","pos.t.value", "neg.coeff","neg.t.value")
ve <- data.frame(cbind(item, vect)); colnames(ve) <- colnames(diag)
( res.CI <- rbind(diag, ve)[c(1:2, 13:16, 3:12), ] )
rownames(res.CI) <- 1:nrow(res.CI)

                        
# 4. APT + ECM   _______________________________________________________________

(sem <- ecmSymFit(y=prVi, x=prCh, lag=4)); names(sem)
aem <- ecmAsyFit(y=prVi, x=prCh,lag=4, model="mtar", split=TRUE, thresh=t.mtar)
aem
(ccc <- summary(aem))
coe <- cbind(as.character(ccc[1:19, 2]), 
  paste(ccc[1:19, "estimate"], ccc$signif[1:19], sep=""), ccc[1:19, "t.value"],
  paste(ccc[20:38,"estimate"], ccc$signif[20:38],sep=""), ccc[20:38,"t.value"]) 
colnames(coe) <- c("item", "China.est", "China.t", "Vietnam.est","Vietnam.t")

(edia <- ecmDiag(aem, 3))
(ed <- edia[c(1,6,7,8,9), ])
ed2 <- cbind(ed[,1:2], "_", ed[,3], "_")
colnames(ed2) <- colnames(coe)

(tes <- ecmAsyTest(aem)$out)
(tes2 <- tes[c(2,3,5,11,12,13,1), -1])
tes3 <- cbind(as.character(tes2[,1]), 
  paste(tes2[,2], tes2[,6], sep=''), paste("[", round(tes2[,4],2), "]", sep=''),
  paste(tes2[,3], tes2[,7], sep=''), paste("[", round(tes2[,5],2), "]", sep=''))
colnames(tes3) <- colnames(coe)

(coe  <- data.frame(apply(coe, 2, as.character), stringsAsFactors=FALSE))
(ed2  <- data.frame(apply(ed2, 2, as.character), stringsAsFactors=FALSE))
(tes3 <- data.frame(apply(tes3,2, as.character), stringsAsFactors=FALSE))
table.4 <- data.frame(rbind(coe, ed2, tes3))
table.4

} 
}
\keyword{datasets}