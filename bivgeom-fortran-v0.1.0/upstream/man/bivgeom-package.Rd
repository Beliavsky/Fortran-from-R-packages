\name{bivgeom-package}
\alias{bivgeom-package}
\alias{bivgeom}
\docType{package}
\title{
\packageTitle{bivgeom}
}
\description{
\packageDescription{bivgeom}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{bivgeom}
\packageIndices{bivgeom}
}
\author{
Alessandro Barbiero

Maintainer: Alessandro Barbiero (alessandro.barbiero@unimi.it)
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.

Barbiero, A. (2018) Properties and estimation of a bivariate geometric model with locally constant failure rates, submitted
}

\keyword{package}
\seealso{
\code{\link{dbivgeomRoy}}, \code{\link{rbivgeomRoy}}, \code{\link{estbivgeomRoy}}, \code{\link{FbivgeomRoy}}
}
\examples{
#####################################
#### MONTE CARLO SIMULATION PLAN ####
#####################################
# setting the parameters' values
theta1 <- 0.3
theta2 <- 0.7
theta3 <- 0.6
N <- 20       # number of Monte Carlo runs
n <- 100      # sample size
# arranging the array containig the simulation results
# N runs, 7 methods, 3 estimates
h <- array(0,c(N,7,3))
# setting the seed 
set.seed(12345)
# function for handling missing values 
# when computing the mean and standard deviation of the estimates:
meanrm <- function(x){mean(x,na.rm=TRUE)}
sdrm <- function(x){sd(x,na.rm=TRUE)}
colnames <- c("ML","MMP","MM1","MM2","MM3","MM4","LS")
dimnames(h)[[2]] <- colnames
# Monte Carlo simulation:
for(i in 1:N)
{
d <- rbivgeomRoy(n,theta1,theta2,theta3)
cat("MC run #",i,"\n")
x<-d[,1]
y<-d[,2]
# implementing all the estimation methods 
# and saving the point estimates in the array
h[i,1,] <- estbivgeomRoy(x, y, "ML")
h[i,2,] <- estbivgeomRoy(x, y, "MMP")
h[i,3,] <- estbivgeomRoy(x, y, "MM1")
h[i,4,] <- estbivgeomRoy(x, y, "MM2")
h[i,5,] <- estbivgeomRoy(x, y, "MM3")
h[i,6,] <- estbivgeomRoy(x, y, "MM4")
h[i,7,] <- estbivgeomRoy(x, y, "LS")
}
# printing MC expected values and standard errors
# for each of the proposed estimation methods
cat("hattheta1:","\n")
cbind(mean=apply(h,c(2,3),meanrm)[,1],se=apply(h,c(2,3),sdrm)[,1]) 
cat("hattheta2:","\n")
cbind(mean=apply(h,c(2,3),meanrm)[,2],se=apply(h,c(2,3),sdrm)[,2]) 
cat("hattheta3:","\n")
cbind(mean=apply(h,c(2,3),meanrm)[,3],se=apply(h,c(2,3),sdrm)[,3]) 
# boxplots of MC distribution of the estimators of theta3
boxplot(h[,,3])
abline(h=theta3, lty=3)
}
