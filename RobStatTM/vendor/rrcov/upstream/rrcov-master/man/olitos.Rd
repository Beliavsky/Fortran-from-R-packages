\name{olitos}
\alias{olitos}
\docType{data}
\title{
Olive Oil Data
}
\description{
This dataset consists of 120 olive oil samples on measurements on
25 chemical compositions (fatty acids, sterols, triterpenic alcohols) 
of olive oils from Tuscany, Italy (Armanino et al. 1989).
There are 4 classes corresponding to different production areas.
Class 1, Class 2, Class 3, and Class 4 contain 50, 25, 34, and 11
observations, respectively. 

}
\usage{data(olitos)}
\format{
  A data frame with 120 observations on the following 26 variables.
  \describe{
    \item{\code{X1}}{Free fatty acids}
    \item{\code{X2}}{Refractive index}
    \item{\code{X3}}{K268}
    \item{\code{X4}}{delta K}
    \item{\code{X5}}{Palmitic acid}
    \item{\code{X6}}{Palmitoleic acid}
    \item{\code{X7}}{a numeric vector}
    \item{\code{X8}}{a numeric vector}
    \item{\code{X9}}{a numeric vector}
    \item{\code{X10}}{a numeric vector}
    \item{\code{X11}}{a numeric vector}
    \item{\code{X12}}{a numeric vector}
    \item{\code{X13}}{a numeric vector}
    \item{\code{X14}}{a numeric vector}
    \item{\code{X15}}{a numeric vector}
    \item{\code{X16}}{a numeric vector}
    \item{\code{X17}}{a numeric vector}
    \item{\code{X18}}{a numeric vector}
    \item{\code{X19}}{a numeric vector}
    \item{\code{X20}}{a numeric vector}
    \item{\code{X21}}{a numeric vector}
    \item{\code{X22}}{a numeric vector}
    \item{\code{X23}}{a numeric vector}
    \item{\code{X24}}{a numeric vector}
    \item{\code{X25}}{a numeric vector}
    \item{\code{grp}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
  }
}
\source{
    Prof. Roberto Todeschini, Milano Chemometrics and QSAR Research Group
    \url{https://michem.unimib.it}
}
\references{
    C. Armanino, R. Leardi, S. Lanteri and G. Modi, 1989. 
    Chemometric analysis of Tuscan olive oils. 
    \emph{Cbemometrics and Intelligent Laboratoty Sysiem}, 5: 343--354.
    
    R. Todeschini, V. Consonni, A. Mauri, M. Pavan (2004) 
    Software for the calculation of molecular descriptors. 
    Pavan M. Talete slr, Milan, Italy, \url{http://www.talete.mi.it}
}
\examples{
data(olitos)
cc <- Linda(grp~., data=olitos, method="mcdC", l1med=TRUE)
cc
pr <- predict(cc)
tt <- mtxconfusion(cc@grp, pr@classification, printit=TRUE)
}
\keyword{datasets}
