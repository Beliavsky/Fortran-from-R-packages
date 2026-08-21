\name{outliersInCompositions}
\alias{outliersInCompositions}
\title{Analysing outliers in compositions.}
\description{
  The Philosophy behind outlier treatment in library(compositions).
}
\details{
  Outliers are omnipresent in all kinds of data analysis. To avoid
  catastrophic misinterpreations robust statistics has developed some
  methods to avoid the distracting influence of the outliers. The
  introduction of robust methods into the compositions package is
  described in \link{robustnessInCompositions}.

  However sometimes we are
  interested directly in the analysis of outliers. The central
  philosophy of the the outlier classification subsystem in
  compositions is
  that outlier are in most cases not simply erroneous observations, but
  rather products of some systematic anomality. This can e.g. be an
  error in an individual component, a secondary process or a minor
  undetected but different subpopulation. The package provides various
  concepts to investigate possible reasons for outliers in compositional
  datasets. 
  \describe{
    \item{Proven Outliers}{
      The package relies on an additive--lognormal reference distribution
      in the simplex (and the correponding normal distribution in each
      other scale). The central tool for the detection of outliers is the
      Mahalanobis distance of the observation from a robustly estimated
      center based on a robustly estimated covariance. The robust
      estimation can be influenced by the given robust attributes. An
      outlier is considered as proven if its Mahalanobis distance is
      larger that the (1-alpha) quantile of the distribution of the
      maximum Mahalanobis distance of a dataset of the same size with a
      corresponding
      (additive)(log)normal distribution. This relies heavily on the
      presumption that the robust estimation is invariant under linear
      transformation, but make no assumptions about the actually used
      robust estimation method. The corresponding distributions are thus
      only defined with respect to a specific implementation of the robust
      estimation algorithm. See 
      \code{\link{OutlierClassifier1}(...,type="outlier")},
      \code{\link{outlierplot}(...,type=c("scatter","biplot"),class.type="outlier")}, 
      \code{\link{qMaxMahalanobis}(...)}.
      
    }
    \item{Extrem Values / Possible outliers}{
      Some cases of the dataset might have unusually high Mahalanobis
      distances, e.g. such that we would expect the probility of a random
      case to have such a value or higher might be below alpha. In
      Literature these cases are often rendered as outliers, because this
      level is approximated by the correponding chisq-based criterion
      proposed. However we consider these only as extrem values, but
      however provide tools to detect and plot them. See 
      \code{\link{OutlierClassifier1}(...,type="grade")},
      \code{\link{outlierplot}(...,type=c("scatter","biplot"),class.type="grade")}, 
      \code{\link{qEmpiricalMahalanobis}(...)} 
    }
    \item{Single Component Outliers}{
      Some Outliers can be explained by a single component, e.g. because
      this single measurement error was wrong. These sort of outliers is
      detected when we reduce the dataset to a subcomposition with one
      component less and realise that our former outlier is now a fairly
      normal member of the dataset, maybe not even extrem. Thus a outlier
      is considered as as single component outlier, when it does not
      appear extrem in any of the subcompositions with one component
      less. For other outliers we can prove that they are still extrem for
      all subcomposition with one component removed. Thus these have to be
      as multicomponent outliers, that can not be explained by a single
      measurment error. For remaining single component outliers, we can
      ask which component is able to explain the outlying character. See 
      \code{\link{OutlierClassifier1}(...,type=c("best","type","all"))}.  
    }
    \item{Counting hidden outliers}{
      If outliers are not outlying far enough to be detected by the test
      for outlyingness are only at first sight harmless. One outlier is
      within the reasonable bounds of what a normal distribution could
      have delivered should not harm the analysis and might not even
      detectable in any way. However if there is more than one they could
      akt together to disrupt our analysis and more interestingly there
      might be some joint reason, which than might make them an
      interesting object of investigation in themselfs. Thus the package
      provides methods (e.g. \code{\link{outlierplot}(...,type="portions")}),
      to prove the existence of such outliers, to give a lower bound
      for there number and to provide us with suspects, with an associated
      outlyingness
      probability. See \code{\link{outlierplot}(...,type="portions")},
      \code{\link{outlierplot}(...,type="nout")}, \code{\link{pQuantileMahalanobis}(...)} 
    }
    \item{Finding atypical subpopulations}{
      When we assume smaller subpopulation we need a tool finding these
      clusters. However usual cluster analysis tends to ignore the
      subgroups, split the main mass and then associate the subgroups
      prematurely to the next part of the main mass. For this task we have
      developed special tools to find
      clusters of atypical populations clearly inducing secondary modes,
      without ripping apart the central
      nonoutlying mass. See \code{\link{ClusterFinder1}}.
    }
    \item{Identifying multiple distracting processes}{
      Outliers that are not due to a seperate subpopulation or due to a
      single component error, might still belong together for beeing
      influenced by the same secondary process distorting the composition
      to a different degrees. Out proposal is to cluster the direction of
      the outliers from the center, e.g. by a command like:
      \code{take<-OutlierClassifier1(data,type="grade")!="ok"}
      \code{hc<-hclust(dist(normalize(acomp(scale(data)[take,]))),method="compact")} and to plot by a command like:
      \code{plot(hc)} and \code{plot(acomp(data[take,]),col=cutree(hc,1.5))}
    }
  }  
  With these tools we hope to provide a systematic approach to identify
  various types of outliers in a exploratory analysis. 
}
  
\note{
  The package \pkg{robustbase} is required for using the
  robust estimations and the outlier subsystem of compositions. To
  simplify installation it is not listed as required, but it will be
  loaded, whenever any sort of outlierdetection or robust estimation is
used.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \link{compositions-package}, \link{missingsInCompositions},
  \link{robustnessInCompositions}, \link{outliersInCompositions},
  \code{\link{outlierplot}},
  \code{\link{OutlierClassifier1}}, \code{\link{ClusterFinder1}}
}
\references{
K. Gerald van den Boogaart, Raimon Tolosana-Delgado, Matevz-Bren (2009)
Robustness, classification and visualization of outliers in
compositional data, in prep.
}
\examples{
\dontrun{
# To slow
tmp<-set.seed(1400)
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A\%*\%t(A))
Mcenter <- acomp(c(1,2,1))
typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
colnames(typicalData)<-c("A","B","C")
data1 <- acomp(rnorm.acomp(100,Mcenter,Mvar))
data2 <- acomp(rbind(typicalData+rbinom(100,1,p=0.1)*rnorm(100)*acomp(c(4,1,1))))
data3 <- acomp(rbind(typicalData,acomp(c(0.5,1.5,2))))
colnames(data3)<-colnames(typicalData)
tmp<-set.seed(30)
rcauchy.acomp <- function (n, mean, var){
    D <- gsi.getD(mean)-1
    perturbe(ilrInv(matrix(rnorm(n*D)/rep(rnorm(n),D), ncol = D) \%*\% chol(clrvar2ilr(var))), mean)
}
data4 <- acomp(rcauchy.acomp(100,acomp(c(1,2,1)),Mvar/4))
colnames(data4)<-colnames(typicalData)
data5 <- acomp(rbind(unclass(typicalData)+outer(rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2))))
data6 <- acomp(rbind(typicalData,rnorm.acomp(20,acomp(c(4,4,1)),Mvar)))
datas <- list(data1=data1,data2=data2,data3=data3,data4=data4,data5=data5,data6=data6)
tmp <-c()
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

moreData <- acomp(rbind(data3,data5,data6))
take<-OutlierClassifier1(moreData,type="grade")!="ok"
hc<-hclust(dist(normalize(acomp(scale(moreData)[take,]))),method="complete")
plot(hc)
plot(acomp(moreData[take,]),col=cutree(hc,1.5))
}
}
\keyword{robust}
