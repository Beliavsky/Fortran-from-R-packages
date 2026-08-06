\name{soil}
\alias{soil}
\docType{data}
\title{
Exchangable cations in forest soil data set
}
\description{
The forest soil data set contains measurements on 58 soil 
pits in the Hubbard Brook Experimental Forest in north-central
New Hampshire. The excavations were done in 1983 and 1986. The 
soil samples were analyzed for the exchangeable cations of 
aluminium, calcium, magnesium, potassium and sodium. The pit 
locations in both data sets can be classified by the 
type of the forest: 

\itemize{
    \item 1: spruce-fir (11 samples),
    \item 2: high elevation hardwood (23 samples) and
    \item 3: low elevation hardwood (24 samples)).
}

Additionally the degree of logging disturbance can be considered
(all 0 in the 1983 data set):

\itemize{
    \item 0: uncut forest,
    \item 1: cut, undisturbed by machinery and
    \item 2: cut, disturbed.
}

The observations are expressed in grams of exchangeable 
cations per square meter.

}
\usage{data(soil)}
\format{
  A data frame with 116 observations on the following 7 variables.
  \describe{
    \item{\code{F}}{Type of forest}
    \item{\code{D}}{Degree of logging disturbance}
    \item{\code{Al}}{Level of the exchangable cations in Al}
    \item{\code{Ca}}{Level of the exchangable cations in Ca}
    \item{\code{Mg}}{Level of the exchangable cations in Mg}
    \item{\code{K}}{Level of the exchangable cations in K}
    \item{\code{Na}}{Level of the exchangable cations in Na}
  }
}
\source{
    Morrison D.F., 2005, Multivariate Statistical Methods, Thompson
}
\references{
    Vanden Branden K, Hubert M (2005). 
    Robust Classiffication in High Dimensions Based on the SIMCA Method. 
    \emph{Cbemometrics and Intelligent Laboratoty Sysiem}, 79: 10--21.   
}
\examples{
data(soil)
soil1983 <- soil[soil$D == 0, -2]       # only 1983, remove column D (always 0)

(cc <- Linda(F~., data=soil))
(pr <- predict(cc))
pr@classification

}
\keyword{datasets}
