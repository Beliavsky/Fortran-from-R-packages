##  PCA for skewed data
##
##  The present robust PCA methods like ROBPCA work best if the
##  non-outlying data have an approximately symmetric distribution.
##  When the original variables are skewed, too many points tend
##  to be flagged as outlying. Hubert et al. (2009) developed a
##  robust method which is also suitable for skewed data. Also, the
##  outlier map is modified to present adequately the PCA outliers.
##
##  In PcaHubert() the version for skewed data is invoked using the
##  argument skew=TRUE. Example with the computer hardware data
##  follows

library(rrcov)

data(machines)
##  The data set contains 209 observations, and 8 variables.

##  First we robustly center and scale the variables by subtracting
##  the median and dividing by the median absolute deviation (MAD).
##  Of course, this could be done also later, when calling the PCA
##  function, by setting the arguments 'center' and 'sacle' to TRUE.
data <- robustbase::doScale(machines, center=median, scale=mad)
X <- data$x

##  Each of the variables is significantly skewed as measured by
##  its medcouple. The corresponding p-values can be computed using
##  the formulas in using the formulas in Brys et al. (2004). These
##  show the variables are all significantly asymmetric.
data.frame(MC=round(apply(X, 2, mc),2))

## Plot a pairwise scaterplot matrix
mcd <- CovMcd(X[,1:6])
plot(mcd, which="pairs")

##  Remove the rownames (too long)
rownames(X) <- NULL

## Start with robust PCA based on MCD (P << n)
(pca1 <- PcaHubert(X, k=3))
plot(pca1, main="ROBPCA-MCD", off=0.03)

## PCA with the projection algoritm of Hubert
(pca2 <- PcaHubert(X, k=3, mcd=FALSE))
plot(pca2, main="ROBPCA-SD", off=0.03)

## PCA with the adjusted for skewness algorithm of Hubert et al (2009)
(pca3 <- PcaHubert(X, k=3, mcd=FALSE, skew=TRUE))
plot(pca3, main="ROBPCA-AO", off=0.03)
