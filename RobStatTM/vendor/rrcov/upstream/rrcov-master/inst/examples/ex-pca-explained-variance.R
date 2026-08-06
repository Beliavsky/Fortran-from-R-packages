##  Showing correctly the percentage explained variance in PCA
##
##  If not all PCA were extracted, which is a great advantage in
##  the case in high dimensional data, the PCA methods in rrcov
##  could not show correctly the percentage of variance explained,
##  Because the total variance explained, i.e. the sum of _all_
##  eigenvalues was not known. Now this is fixed, differently in
##  the different methods.

##  In PcaClassic, PcaCov and PcaLocantore, this is nt a problem
##  because always all eigenvalues are computed.
##
##  In PcaHubert there is a preliminary step in which the classical
##  PCA are computed on the data set without outliers identified by
##  the Stahel-Donoho Outlyingness. All eigenvalues are calculated
##  and used for selecting the number of components and for
##  presenting the percentage of explaned variance.
##
##  In the pure projection purcuit methods PcaGrid() and PcaProj()
##  this cannot be done and a note is written that the proportion
##  of variance and cumulative proportion are not shown because the
##  chosen number of components is smaller than the rank of the data
##  matrix.

library(rrcov)

data(hbk)

## PCA with all variables
(pca1 <- PcaHubert(hbk, k=ncol(hbk), trace=TRUE, mcd=TRUE, skew=FALSE))
summary(pca1)

## PCA with number of components selected by the algorithm
(pca2 <- PcaHubert(hbk, trace=TRUE, mcd=TRUE, skew=FALSE))
summary(pca2)

## PCA with number of components selected by the user
(pca3 <- PcaHubert(hbk, k=2, trace=TRUE, mcd=TRUE, skew=FALSE))
summary(pca3)

## PCA by the projection algorithm with number of components selected by the user.
##  Here we cannot show the proportion of variance and the cumulative proportion
(pca4 <- PcaGrid(hbk, k=2, trace=TRUE))
summary(pca4)

##  The other PCA methods available in rrcov
summary(PcaClassic(hbk, k=2))
summary(PcaCov(hbk, k=2))
summary(PcaLocantore(hbk, k=2))
summary(PcaProj(hbk, k=2))

## Example with the newly added to 'rrcov' data set fruit ========
data(fruit)
# Remove the first variable, the grouping one
(pca <- PcaHubert(fruit[,-1], trace=TRUE))
summary(pca)

screeplot(pca)

(pca <- PcaHubert(fruit[,-1], k=4)
summapry(pca)
plot(pca)
