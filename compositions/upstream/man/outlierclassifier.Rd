\name{outlierclassifier}
\alias{OutlierClassifier1}
\alias{OutlierClassifier1.acomp}
\title{Detect and classify compositional outliers.}
\description{
  Detects outliers and classifies them according to different possible
  explanations.
}
\usage{
OutlierClassifier1(X,...)
\method{OutlierClassifier1}{acomp}(X,...,alpha=0.05,
           type=c("best","all","type","outlier","grade"),goodOnly=NULL,
           corrected=TRUE,RedCorrected=FALSE,robust=TRUE) 
}
\arguments{
  \item{X}{the dataset as an \code{acomp} object}
  \item{\dots}{further arguments to MahalanobisDist/gsi.mahOutlier}
  \item{alpha}{The confidence level for identifying outliers.}
  \item{type}{ What type of classification should be used: best: Which
    single component would best explain the outlier. all: Give a binary coding
    specifying all components, which could explain the outlier. type: Is
  it a a normal observation \code{"ok"}, a single componentn outlier \code{"1"} or can it not
  be explained by a single wrong component "?". outlier: All outliers 
  are marked as \code{"outlier"}, others are marked as
  \code{"ok"}. \code{"grade"}: Proven Outliers are marked as "outlier"s,
suspected outliers, detected without correction of the p-value are
reported as "extreme", the rest is reported as "ok".}
\item{goodOnly}{an integer vector. Only the specified index of the
dataset should be used
  for estimation of the outlier criteria. This parameter if only a small
  portion of the dataset is reliable.}
\item{corrected}{logical. Literatur often proposed to compare the
  Mahalanobis distances with Chisq-Approximations of there
  distributions. However this does not correct for multiple
  testing. If corrected is true a correction for multiple testing is
  used. In any case we do not use the chisq-approximation, but a
  simulation based procedure to compute confidence bounds. }
\item{RedCorrected}{logical. If an outlier is detected we can try to
  find out wether a single component would be sufficient to drop the
  outlier under the outlier detection limit. Since in this second case
  we only check a few outliers no second correction step applies as long
as the number of outliers is not very high.}
\item{robust}{A robustness description as define in
  \code{\link{var.acomp}}}
}
\value{
A factor classifying the observations in the dataset as "ok" or some
type of outlier. 
  
}
\details{
  See \link{outliersInCompositions} for a comprehensive introduction into the outlier
  treatment in compositions.
  
  See \code{\link{ClusterFinder1}} for an alternative method to classify
  observations in the context of outliers.  
}
\note{
  The package \pkg{robustbase} is required for using the
  robust estimations.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{outlierplot}}, \code{\link{ClusterFinder1}}
}

\examples{
\dontrun{
tmp<-set.seed(1400)
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A\%*\%t(A))
Mcenter <- acomp(c(1,2,1))
data(SimulatedAmounts)
datas <- list(data1=sa.outliers1,data2=sa.outliers2,data3=sa.outliers3,
              data4=sa.outliers4,data5=sa.outliers5,data6=sa.outliers6)
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
tmp<-mapply(function(x,y) {
outlierplot(x,type="scatter",class.type="grade");
  title(y)
},datas,names(datas))


par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
tmp<-mapply(function(x,y) {
  myCls2 <- OutlierClassifier1(x,alpha=0.05,type="all",corrected=TRUE)
  outlierplot(x,type="scatter",classifier=OutlierClassifier1,class.type="best",
  Legend=legend(1,1,levels(myCls),xjust=1,col=colcode,pch=pchcode),
  pch=as.numeric(myCls2));
  legend(0,1,legend=levels(myCls2),pch=1:length(levels(myCls2)))
  title(y)
},datas,names(datas))

par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="ecdf",main=names(datas)[i])
par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="portion",main=names(datas)[i])
par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="nout",main=names(datas)[i])
par(opar)
}
}
\keyword{multivariate}
