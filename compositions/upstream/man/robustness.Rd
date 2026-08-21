\name{robustnessInCompositions}
\alias{robust}
\alias{robustnessInCompositions}
\title{Handling robustness issues and outliers in compositions.}
\description{
  The seamless transition to robust estimations in library(compositions).
}
\details{
  A statistical method is called nonrobust if an arbitrary contamination of
  a small portion of the dataset can produce results radically different
  from the results without the contamination. In this sense many
  classical procedures relying on distributional models or on moments
  like mean and variance are highly nonrobust. 
  \cr
  We consider robustness as an essential prerequierement of all
  statistical analysis. However in the context of compositional data
  analysis robustness is still in its first years. 
  \cr
  As of Mai 2008 we provide a new approach to robustness in the
  package. The central idea is that robustness should be more or less
  automatic and that there should be no necessity to change the code to
  compare results obtained from robust procedures and results from there
  more efficient nonrobust counterparts. 
  \cr
  To achieve this all routines that rely on distributional models (such
  as e.g. mean,
  variance, principle component analysis, scaling) and routines relying
  on those routines get a new standard argument of the form:
  \cr
  \code{fkt(...,robust=getOption("robust"))}
  \cr
  which defaults to a new option "robust". This option can take several
  values:
  \describe{
    \item{FALSE}{The classical estimators such as arithmetic mean and
      persons product moment variance are used and the results are to be
      considered nonrobust. }
    \item{TRUE}{The default for robust estimation in the package is
      used. At this time this is \code{\link[robustbase]{covMcd}} in the
      \pkg{robustbase}-package. This default might change in future. }
    \item{"pearson"}{This is a synonym for FALSE and explicitly states
      that no robustness should be used.}
    \item{"mcd"}{Minimum Covariance Determinant. This option explicitly
      selects the use of \code{\link[robustbase]{covMcd}} in the
      \pkg{robustbase}-package as the main robustness engine.}
  }
  More options might follow later.
  To control specific parameters of the
  model the string can get an attribute named "control" which contains
  additional options for the robustness engine used. In this moment the
  control attribute of mcd is a control object of
  \code{\link[robustbase]{covMcd}}. The control argument of "pearson" is a list
  containing addition options to the mean, like trim. 
  \cr
  The standard value for getOption("robust") is FALSE to avoid situation
  in which the user thinks he uses a classical technique. Robustness
  must be switched on explicitly. Either by setting the option with
  \code{options(robust=TRUE)} or by giving the argument. This default
  might change later if the authors come to the impression that robust
  estimation is now considered to be the default.
  
  For those not only interested in avoiding the influence of the
  outliers, but in an analysis of the outliers we added a subsystem for
  outlier classification. This subsystem is described in
  \link{outliersInCompositions} and also relies on the
  robust option. However evidently for these routines the factory
  default for the robust option is always TRUE, because it is only
  applicable in an outlieraware context.
  
  We hope that in this way we can provide a seamless transition from
  nonrobust analysis to a robust analysis. 
}
\note{
  IMPORTANT: The robust argument only works with the classes of the
  package. Only your compositional analysis is suddenly robust.
  \cr
  The package \pkg{robustbase} is required for using the
  robust estimations and the outlier subsystem of compositions. To
  simplify installation it is not listed as required, but it will be
  loaded, whenever any sort of outlierdetection or robust estimation is
used.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{var.acomp}}, \code{\link{mean.acomp}},
  \pkg{robustbase}, \link{compositions-package},
\link{missings}, \code{\link{outlierplot}},
  \code{\link{OutlierClassifier1}}, \code{\link{ClusterFinder1}}
}
\examples{
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A\%*\%t(A))
Mcenter <- acomp(c(1,2,1))
typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
colnames(typicalData)<-c("A","B","C")
data5 <- acomp(rbind(unclass(typicalData)+outer(rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2))))

mean(data5)
mean(data5,robust=TRUE)
var(data5)
var(data5,robust=TRUE)
Mvar
biplot(princomp(data5))
biplot(princomp(data5,robust=TRUE))

}
\keyword{robust}
