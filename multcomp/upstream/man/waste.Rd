\name{waste}
\alias{waste}
\docType{data}
\title{ Industrial Waste Data Set }
\usage{data("waste")}
\description{
  Industrial waste output in a manufactoring plant.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{temp}{temperature, a factor at three levels: \code{low}, \code{medium}, \code{high}.}
    \item{envir}{environment, a factor at five levels: \code{env1} \dots \code{env5}.}
    \item{waste}{response variable: waste output in a manufacturing plant.}
  }
}
\details{

  The data are from an experiment designed to study the effect of temperature
  (\code{temp}) and environment (\code{envir}) on waste output in a manufactoring plant.
  Two replicate measurements were taken at each temperature / environment combination.
  Data were obtained from \bibcitet{multcomp::westfall1999}, page 177.

}
\references{\bibshow{*}}
\examples{

  ### set up two-way ANOVA with interactions
  amod <- aov(waste ~ temp * envir, data=waste)

  ### comparisons of main effects only
  K <- glht(amod, linfct = mcp(temp = "Tukey"))$linfct
  K
  glht(amod, K)

  ### comparisons of means (by averaging interaction effects)
  low <- grep("low:envi", colnames(K))
  med <- grep("medium:envi", colnames(K))
  K[1, low] <- 1 / (length(low) + 1)
  K[2, med] <- 1 / (length(low) + 1)
  K[3, med] <- 1 / (length(low) + 1)
  K[3, low] <- - 1 / (length(low) + 1)
  K
  confint(glht(amod, K))

  ### same as TukeyHSD
  TukeyHSD(amod, "temp")

  ### set up linear hypotheses for all-pairs of both factors
  wht <- glht(amod, linfct = mcp(temp = "Tukey", envir = "Tukey"))

  ### cf. Westfall et al. (1999, page 181)
  summary(wht, test = adjusted("Shaffer"))

}
\keyword{datasets}
