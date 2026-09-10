\name{makewpstRO}
\alias{makewpstRO}
\title{
Make a wavelet packet regression object from a dependent and independent
time series variable.
}
\description{
The idea here is to try and build facilities to enable
a transfer function model along
the lines of that described by Nason and Sapatinas 2002 in
\emph{Statistics and Computing}. The idea is to turn the
\code{timeseries} variable into a set of nondecimated wavelet packets
which are already pre-selected to have some semblance of relationship
to the \code{response} time series. The function does not actually
perform any regression, in contrast to the related \code{\link{makewpstDO}}
but returns a data frame which the user can use to build their own models.
}
\usage{
makewpstRO(timeseries, response, filter.number = 10,
    family = "DaubExPhase", trans = logabs, percentage = 10)
}
\arguments{
  \item{timeseries}{
	The dependent variable time series. This series is decomposed
	using the \code{\link{wpst}} function into nondecimated
	wavelet packets, need to be a power of two length.
}
  \item{response}{
	The independent or response time series.
}
  \item{filter.number}{
	The type of wavelet used within \code{family}, see
	\code{\link{filter.select}}.
}
  \item{family}{
	The family of wavelet, see \code{\link{filter.select}}
}
  \item{trans}{
	A transform to apply to the nondecimated wavelet packet
	coefficients before any selection
}
  \item{percentage}{
	The top \code{percentage} of nondecimated wavelet packets
	that correlated best with the \code{response} series will
	be preselected.
}
}
\details{The idea behind this methodology is that a \code{response}
time series might not be directly related to the dependent
\code{timeseries} time series, but it might be related to the
nondecimated wavelet packets of the \code{timeseries}, these packets
can pick out various features of the \code{timeseries} including
certain delays, oscillations and others. 

The best packets (the number if controlled by \code{percentage}), those
that correlate best with \code{response} are selected and returned.
The \code{response} and the best nondecimated wavelet packets
are returned in a data frame object and then any convenient form of
statistical modeling can be used to build a model of the \code{response} in
terms of the packet variables. 

Once a model has been built it can be interpreted in the usual way, but
with respect to nondecimated wavelet packets.

Note that nondecimated wavelet packets are essential, as they are all
of the same length as the original response series. If a decimated
wavelet packet algorithm had been used then it is not clear what to
do with the "gaps"!

If new \code{timeseries} data comes along the \code{\link{wpstREGR}}
function can be used to extract the identical packets as the ones
produced by this function (as the result of this function stores the
identities of these packets). Then the statistical modelling that
build the model from the output of this function, can be used to
predict future values of the \code{response} time series from future
values of the \code{timeseries} series.
}
\value{
An object of class \code{wpstRO} containing the following items
\item{df}{A data frame containing the \code{response} time series and
	a number of columns/variables/packets that correlated with
	response series. These are all entitled "Xn" where n is some
	integer} 
\item{ixvec}{A packet index vector. After taking the nondecimated wavelet
	packet transform, all the packets are stored in a matrix. This
	vector indicates those that were preselected}
\item{level}{The original level from which the preselected vectors came from}
\item{pktix}{Another index vector, this time referring to the original
	wavelet packet object, not the matrix in which they subsequently
	got stored}
\item{nlevelsWT}{The number of resolution levels in the original wavelet
	packet object}
\item{cv}{The correlation vector. These are the values of the correlations
	of the packets with the response, then sorted in terms of decreasing
	absolute correlation}
\item{filter}{The wavelet filter details}
\item{trans}{The transformation function actually used}
}
\references{
Nason, G.P. and Sapatinas, T. (2002) Wavelet packet transfer function
	modeling of nonstationary time series. \emph{Statistics and Computing},
	\bold{12}, 45-56.
}
\author{
G P Nason
}
\seealso{
\code{\link{makewpstDO}}, \code{\link{wpst}}, \code{\link{wpstREGR}}}
\examples{
data(BabyECG)
baseseries <- BabyECG[1:256]
#
# Make up a FICTITIOUS response series!
#
response <- BabyECG[6:261]*3+52
#
# Do the modeling
#
BabeModel <- makewpstRO(timeseries=baseseries, response=response)
#Level: 0  ..........
#1  ..........
#2  ..........
#3  ..........
#4  ................
#5  
#6  
#7  
#
#Contains SWP coefficients
#Original time series length:  256 
#Number of bases:  25 
#Some basis selection performed
#       Level Pkt Index Orig Index      Score
#[1,]     5         0        497  0.6729833
#[2,]     4         0        481  0.6120771
#[3,]     6         0        505  0.4550616
#[4,]     3         0        449  0.4309924
#[5,]     7         0        509  0.3779385
#[6,]     1        53        310  0.3275428
#[7,]     2        32        417 -0.3274858
#[8,]     2        59        444 -0.2912863
#[9,]     3        16        465 -0.2649679
#[10,]     1       110        367  0.2605178
#etc. etc.
#
#
# Let's look at the data frame component
#
names(BabeModel$df)
# [1] "response" "X1"       "X2"       "X3"       "X4"       "X5"      
# [7] "X6"       "X7"       "X8"       "X9"       "X10"      "X11"     
#[13] "X12"      "X13"      "X14"      "X15"      "X16"      "X17"     
#[19] "X18"      "X19"      "X20"      "X21"      "X22"      "X23"     
#[25] "X24"      "X25"    
#
# Generate a formula including all of the X's (note we could use the .
# argument, but we later want to be more flexible
#
xnam <- paste("X", 1:25, sep="")
fmla1 <- as.formula(paste("response ~ ", paste(xnam, collapse= "+")))
#
# Now let's fit a linear model, the response on all the Xs
#
Babe.lm1 <- lm(fmla1, data=BabeModel$df)
#
# Do an ANOVA to see what's what
#
anova(Babe.lm1)
#Analysis of Variance Table
#
#Response: response
#	Df Sum Sq Mean Sq  F value    Pr(>F)    
#X1          1 214356  214356 265.7656 < 2.2e-16 ***
#X2          1  21188   21188  26.2701 6.289e-07 ***
#X3          1  30534   30534  37.8565 3.347e-09 ***
#X4          1    312     312   0.3871 0.5344439    
#X5          1   9275    9275  11.4999 0.0008191 ***
#X6          1     35      35   0.0439 0.8343135    
#X7          1    195     195   0.2417 0.6234435    
#X8          1     94      94   0.1171 0.7324600    
#X9          1    331     331   0.4103 0.5224746    
#X10         1      0       0   0.0006 0.9810560    
#X11         1    722     722   0.8952 0.3450597    
#X12         1      0       0   0.0004 0.9850243    
#X13         1     77      77   0.0959 0.7570769    
#X14         1   2770    2770   3.4342 0.0651404 .  
#X15         1      6       6   0.0072 0.9326155   
#X16         1    389     389   0.4821 0.4881649    
#X17         1     44      44   0.0544 0.8157015    
#X18         1     44      44   0.0547 0.8152640    
#X19         1   4639    4639   5.7518 0.0172702 *  
#X20         1    490     490   0.6077 0.4364469    
#X21         1    389     389   0.4823 0.4880660    
#X22         1     85      85   0.1048 0.7463860    
#X23         1   1710    1710   2.1198 0.1467664    
#X24         1     12      12   0.0148 0.9033427    
#X25         1     82      82   0.1019 0.7498804    
#Residuals 230 185509     807                       
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#
# Looks like X1, X2, X3, X5, X14 and X19 are "significant". Also throw in
# X4 as it was a highly ranked preselected variable, and refit
#
fmla2 <- response ~ X1 + X2 + X3 + X4 + X5 + X14 + X19
Babe.lm2 <- lm(fmla2, data=BabeModel$df)
#
# Let's see the ANOVA table for this
#
anova(Babe.lm2)
#Analysis of Variance Table
#
#Response: response
#	Df Sum Sq Mean Sq  F value    Pr(>F)    
#X1          1 214356  214356 279.8073 < 2.2e-16 ***
#X2          1  21188   21188  27.6581 3.128e-07 ***
#X3          1  30534   30534  39.8567 1.252e-09 ***
#X4          1    312     312   0.4076 0.5238034    
#X5          1   9275    9275  12.1075 0.0005931 ***
#X14         1   3095    3095   4.0405 0.0455030 *  
#X19         1   4540    4540   5.9259 0.0156263 *  
#Residuals 248 189989     766                       
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#
# So, let's drop X4, refit, and then do ANOVA
#
Babe.lm3 <- update(Babe.lm2, . ~ . -X4)
anova(Babe.lm3)
#
# After viewing this, drop X14
#
Babe.lm4 <- update(Babe.lm3, . ~ . -X14)
anova(Babe.lm4)
#
# Let's plot the original series, and the "fitted" one
#
\dontrun{ts.plot(BabeModel$df[["response"]])}
\dontrun{lines(fitted(Babe.lm4), col=2)}
#
# Let's plot the wavelet packet basis functions associated with the model
#
\dontrun{oldpar <- par(mfrow=c(2,2))}
\dontrun{z <- rep(0, 256)}
\dontrun{zwp <- wp(z, filter.number=BabeModel$filter$filter.number,
    family=BabeModel$filter$family)}
\dontrun{draw(zwp, level=BabeModel$level[1], index=BabeModel$pktix[1], main="", sub="")}
\dontrun{draw(zwp, level=BabeModel$level[2], index=BabeModel$pktix[2], main="", sub="")}
\dontrun{draw(zwp, level=BabeModel$level[3], index=BabeModel$pktix[3], main="", sub="")}
\dontrun{draw(zwp, level=BabeModel$level[5], index=BabeModel$pktix[5], main="", sub="") }
\dontrun{par(oldpar)}
#
# Now let's do some prediction of future values of the response, given
# future values of the baseseries
#
newseries <- BabyECG[257:512]
#
# Get the new data frame
#
newdfinfo <- wpstREGR(newTS = newseries, wpstRO=BabeModel)
#
# Now use the best model (Babe.lm4) with the new data frame (newdfinfo)
# to predict new values of response
#
newresponse <- predict(object=Babe.lm4, newdata=newdfinfo)
#
# What is the "true" response, well we made up a response earlier, so let's
# construct the true response for this future data (in your case you'll
# have a separate genuine response variable)
#
trucfictresponse <- BabyECG[262:517]*3+52
#
# Let's see them plotted on the same plot
#
\dontrun{ts.plot(trucfictresponse)}
\dontrun{lines(newresponse, col=2)}
#
# On my plot they look tolerably close!
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ts}
