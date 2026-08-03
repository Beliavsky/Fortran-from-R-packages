library(rrcov)

data(pottery)
x <- pottery[,c("MG", "CA")]
grp <- pottery$origin

##
## Compute robust location and covariance matrix and
## plot the tolerance ellipses
library(rrcov)
(mcd <- CovMcd(x))
col <- c(3,4)
gcol <- ifelse(grp == "Attic", col[1], col[2])
gpch <- ifelse(grp == "Attic", 16, 1)
plot(mcd, which="tolEllipsePlot", class=TRUE, col=gcol, pch=gpch)


##
## Perform classical LDA and plot the data, 0.975 tolerance ellipses
##  and LDA separation line
##
x <- pottery[,c("MG", "CA")]
grp <- pottery$origin
lda <- LdaClassic(x, grp)
lda

## Calculate the ellipses using the package 'ellipse'
require(ellipse)
e1 <- ellipse(x=lda@cov, centre=lda@center[1,], level=0.975)
e2 <- ellipse(x=lda@cov, centre=lda@center[2,], level=0.975)

## Calculate the ellipses using the package 'rrcov'
e3 <- getEllipse(loc=lda@center[1,], cov=lda@cov)
e4 <- getEllipse(loc=lda@center[2,], cov=lda@cov)

## Calculate the ellipses using the package 'cluster'
require(cluster)
e5 <- ellipsoidPoints(lda@cov, qchisq(0.975, 2), lda@center[1,])
e6 <- ellipsoidPoints(lda@cov, qchisq(0.975, 2), lda@center[2,])

plot(CA~MG, data=pottery, col=gcol, pch=gpch,
    xlim=c(min(MG,e1[,1], e2[,1]), max(MG,e1[,1], e2[,1])),
    ylim=c(min(CA,e1[,2], e2[,2]), max(CA,e1[,2], e2[,2])))

ab <- lda@ldf[1,] - lda@ldf[2,]
cc <- lda@ldfconst[1] - lda@ldfconst[2]
abline(a=-cc/ab[2], b=-ab[1]/ab[2], col=2, lwd=2)

lines(e1, type="l", col=1)
lines(e2, type="l", col=1)
lines(e3, type="l", col=2)
lines(e4, type="l", col=2)
lines(e5, type="l", col=3)
lines(e6, type="l", col=3)
