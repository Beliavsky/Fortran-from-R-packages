\name{outlierplot}
\alias{outlierplot}
\alias{outlierplot.acomp}
\title{Plot various graphics to analyse outliers.}
\description{
  A collection of plots emphasing different aspects of possible outliers.
}
\usage{
outlierplot(X,...)
\method{outlierplot}{acomp}(X,colcode=colorsForOutliers1,
  pchcode=pchForOutliers1,
  type=c("scatter","biplot","dendrogram","ecdf","portion","nout","distdist"),
  legend.position,pch=19,...,clusterMethod="ward",
  myCls=classifier(X,alpha=alpha,type=class.type,corrected=corrected),
  classifier=OutlierClassifier1,
  alpha=0.05,
  class.type="best",
  Legend,pow=1,
  main=paste(deparse(substitute(X))),
  corrected=TRUE,robust=TRUE,princomp.robust=FALSE,
                              mahRange=exp(c(-5,5))^pow,
                              flagColor="red",
                              meanColor="blue",
                              grayColor="gray40",
                              goodColor="green",
                              mahalanobisLabel="Mahalanobis Distance"
                              )
}
\arguments{
  \item{X}{The dataset as an \code{acomp} object}
  \item{colcode}{A color palette for factor given by the \code{myCls},
    or function to create it from the factor. Use \code{colorForOutliers2} if
    \code{class.method="all"} is used. }
  \item{pchcode}{A function to create a plot character palette for the factor
    returned by the \code{myCls} call}
  \item{type}{ The type of plot to be produced. See details for more
    precise definitions.}
  \item{legend.position}{The location of the legend. Must!!! be given to
    draw a classical legend.}
  \item{pch}{A default plotting char}
  \item{\dots}{Further arguments to the used plotting function}
  \item{clusterMethod}{The clustering method for \code{\link{hclust}}
    based outlier grouping.}
  \item{myCls}{A factor presenting the groups of outliers}
  \item{classifier}{The routine to create a factor presenting the groups
    of outliers heuristically. It is only used in the default argument
    to \code{myCls}.}
  \item{alpha}{The confidence level to be used for outlier
    classification tests}
  \item{class.type}{The type of classification that should be generated
    by \code{classifier}}
  \item{Legend}{The content will be substituted and stored as list entry
    legend in the result of the function. It can than be evaluated to
    actually create a seperate legend on another device (e.g. for
    publications).}
  \item{pow}{The power of Mahalanobis distances to be used.}
  \item{main}{The title of the graphic}
  \item{corrected}{Literature typically proposes to compare the
    Mahalanobis distances with the distribution of a random Mahalanobis
    distance. However it would be needed to correct this for (dependent)
    multiple testing, since we always test the whole dataset, which means
    comparing against the distribution of the maximum Mahalanobis
    distance. This argument switches to this second behavior, giving less
    outliers.}
  \item{robust}{A robustness description as define in
    \code{\link{robustnessInCompositions}}}
  \item{princomp.robust}{Either a logical determining wether or not the
    principal component analysis should be done robustly or a principal
    component object for the dataset.}
  \item{mahRange}{The range of Mahalanobis distances displayed. This is
    fixed to make views comparable among datasets. However if the preset
    default is not enough a warning is issued and a red mark is drawn in
    the plot}
  \item{flagColor}{The color to draw critical situations.}
  \item{meanColor}{The color to draw typical curves.}
  \item{goodColor}{The color to draw confidence bounds.}
  \item{grayColor}{The color to draw less important things.}
  \item{mahalanobisLabel}{The axis label to be used for axes displaying
    Mahalanobis distances.}
}
\value{
  a list respresenting the criteria computed to create the plots. The
  content of the list depends on the plotting type selected.
}
\details{
  See \link{outliersInCompositions} for a comprehensive introduction
into the outlier
  treatment in compositions.
  \describe{
    \item{\code{type="scatter"}}{
      Produces an appropriate standard plot such as a tenary diagram with
      the outliers marked by there codes according to the given classifier
      and colorcoding and pch coding.
      \cr
      This shows the actual values of the identified outliers.
    }
    \item{\code{type="biplot"}}{
      Creates a biplot based on a nonrobust principal component analysis
      showing the outliers classified through outliers in the given color
      scheme. We use the nonrobust principal component analyis since it
      rotates according to a good visibility of the extreme values.
      \cr
      This shows the position of the outliers in the usual principal
      components analysis. However note that a \code{\link{coloredBiplot}}
      is used rather than the usual one. 
    }
    %%  \item{\code{type="barplot"}}{
    %%  I don't understand that infl.Mahalanobis stuff.    
    %%  }
    \item{\code{type="dendrogram"}}{
      Shows a dendrogram based on robust Mahalanobis distance
      based hierachical clustering, where the observations are labeled
      with the identified outlier classes. 
      \cr
      This plot can be used to see how good different categories of
      outliers cluster.
    }
    \item{\code{type="ecdf"}}{
      This plot provides a cummulated distribution function of the
      Mahalanobis distances along with an expeced curve and a lower
      confidence limit. The empirical cdf is plotted in the default
      color. The expected cdf is displayed in \code{meanColor}. The
      \code{alpha}-quantile -- i.e. a lower prediction bound -- for the
      cdf is given in goodColor. A line in \code{grayColor} show the
      minium portion of observations above some limit to be
      outliers, based on the portion of observations necessary to move
      down to make the empirical distribution function get above its lower
      prediction limit under the assumption of normality.
      \cr
      This plot shows the basic construction for the minimal number of
      outlier computation done in \code{type="portion"}. 
    }
    \item{\code{type="portion"}}{
      This plot focusses on numbers of outliers. The horizontal axis
      give Mahalanobis distances and the vertical axis number of
      observations. In \code{meanColor} we see a curve of an estimated
      number of outliers above some limit, generated by estimating the
      portion of outliers with a Mahalanobis distance over the given
      limit by max(0,1-ecdf/cdf). The minimum
      number of outliers is computed by replacing cdf by its lower
      confidence limit and displayed in \code{goodColor}. The
      Mahalanobis distances of the individual data points are added as a
      stacked \code{\link{stripchart}}, such that the influence of
      individual observations can be seen.
      \cr
      The true problem of outlier detection is to detect "near"
      outliers. Near outliers are outliers so near to the dataset that
      they could well be extrem observation. These near outliers would
      provide no problem unless they are not many showing up in
      groups. Graphic allows at least to count them and to show there
      probable Mahalanobis distance such, however it still does not
      allow to conclude that an individual observation is an
      outlier. However still the outlier candidates can be identified
      comparing their mahalanobis distance (returned by the plot
      as\code{$mahalanobis}) with a cutoff inferred from this graphic. 
    }
    \item{\code{type="nout"}}{
      This is a simplification of the previous plot simply providing the
      number of outliers over a given limit.
      \cr
      %??? MORE DOCUMENTATION NEEDED ???
    }
    \item{\code{type="distdist"}}{
      Plots a scatterplot of the the classical and robust Mahalanobis
      distance with the given classification for colors and plot
      symbols. Furthermore it plots a horizontal line giving the 0.95-Quantil
      of the distribution of the maximum robust Mahalanobis distance of
      normally distributed dataset.
    }
  }
}
\note{
  The package \pkg{robustbase} is required for using the
  robust estimations.}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{OutlierClassifier1}}, \code{\link{ClusterFinder1}}
}

\examples{
\dontrun{
data(SimulatedAmounts)
outlierplot(acomp(sa.outliers5))

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
# To slow
par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="ecdf",main=names(datas)[i])
par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="portion",main=names(datas)[i])
par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="nout",main=names(datas)[i])
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="distdist",main=names(datas)[i])
par(opar)

}
}
\keyword{multivariate}
