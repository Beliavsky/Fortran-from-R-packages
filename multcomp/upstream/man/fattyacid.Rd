\name{fattyacid}
\alias{fattyacid}
\docType{data}
\title{Fatty Acid Content of Bacillus simplex.}
\description{
  Fatty acid content of different putative ecotypes of Bacillus simplex. 
}
\usage{data("fattyacid")}
\format{
  A data frame with 93 observations on the following 2 variables.
  \describe{
    \item{\code{PE}}{a factor with levels \code{PE3}, \code{PE4}, \code{PE5}, \code{PE6}, \code{PE7}, \code{PE9} indicating the putative ecotype (PE).}
    \item{\code{FA}}{a numeric vector indicating the content of fatty acid (FA).}
    }
}
\details{
The data give the fatty acid content for different putative ecotypes of Bacillus simplex. 
Variances of the values of fatty acid are heterogeneous among the putative ecotypes
\bibcitep{multcomp::sikorski2008}.
}
\references{\bibshow{*}}
\examples{

if (require("sandwich")) {
    data("fattyacid")
    ### all-pairwise comparisons of the means of fatty acid content 
    ### FA between different putative ecotypes PE accounting for 
    ### heteroscedasticity by using a heteroscedastic consistent 
    ### covariance estimation
    amod <- aov(FA ~ PE, data = fattyacid)
    amod_glht <- glht(amod, mcp(PE = "Tukey"), vcov. = vcovHC)
    summary(amod_glht)

    ### simultaneous confidence intervals for the differences of 
    ### means of fatty acid content between the putative ecotypes
    confint(amod_glht)
}
}
\keyword{datasets}
