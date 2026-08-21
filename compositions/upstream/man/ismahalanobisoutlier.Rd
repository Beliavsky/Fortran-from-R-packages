\name{IsMahalanobisOutlier}
\alias{IsMahalanobisOutlier}
\title{Checking for outliers}
\description{
  Detect outliers with respect to a normal distribution model.
}
\usage{
IsMahalanobisOutlier(X,...,alpha=0.05,goodOnly=NULL,
                 replicates=1000,corrected=TRUE,robust=TRUE,crit=NULL)
}
\arguments{
  \item{X}{a dataset (e.g. given as acomp, rcomp, aplus, rplus or rmult) object
    to which \code{\link{idt}} and \code{MahalanobisDist} apply. }
  \item{\dots}{further arguments to MahalanobisDist/gsi.mahOutlier}
  \item{alpha}{The confidence level for identifying outliers.}
  \item{goodOnly}{an integer vector. Only the specified index of the
    dataset should be used
    for estimation of the outlier criteria. This parameter if only a small
    portion of the dataset is reliable.}
  \item{replicates}{The number of replicates to be used in the Monte
    Carlo simulations for determination of the quantiles. The
    \code{replicates} not given a minimum is computed from the alpha
    level to ensure reasonable precission.}
  \item{corrected}{logical. Literatur often proposed to compare the
    Mahalanobis distances with Chisq-Approximations of there
    distributions. However this does not correct for multiple
    testing. If corrected is true a correction for multiple testing is
    used. In any case we do not use the chisq-approximation, but a
    simulation based procedure to compute confidence bounds. }
  \item{robust}{A robustness description as define in
    \code{\link{robustnessInCompositions}}}
  \item{crit}{The critical value to be used. Typically the routine is
    called mainly for the purpose of finding this value, which it does,
    when crit is NULL, however
    sometimes we might want to specifiy a value used by someone else to
    reproduce the results.}
}
\value{
A logical vector giving for each element the result of the alpha-level
test for beeing an outlier. TRUE corresponds to a significant result.   
}
\details{
  See \link{outliersInCompositions} and \link{robustnessInCompositions}
  for a comprehensive introduction into the outlier treatment in
  compositions.
  
  See \code{\link{OutlierClassifier1}} for a highlevel method to
  classify observations in the context of outliers.  
}
\note{
  For some unkown reasons the computation sometimes produces NaN's. In
  this case a warning is issued and a recomputation is tried.
  
  The package \pkg{robustbase} is required for using the
  robust estimations.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{OutlierClassifier1}} , \code{\link{outlierplot}},
  \code{\link{ClusterFinder1}}
}

\examples{
\dontrun{
data(SimulatedAmounts)

datas <- list(data1=sa.outliers1,data2=sa.outliers2,data3=sa.outliers3,
              data4=sa.outliers4,data5=sa.outliers5,data6=sa.outliers6)

opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))
tmp<-mapply(function(x,y){
plot(x,col=ifelse(IsMahalanobisOutlier(x),"red","gray"))
  title(y)
},datas,names(datas))
}
}
\keyword{multivariate}
