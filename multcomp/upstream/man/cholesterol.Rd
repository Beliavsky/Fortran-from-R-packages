\name{cholesterol}
\alias{cholesterol}
\docType{data}
\title{ Cholesterol Reduction Data Set }
\usage{data("cholesterol")}
\description{
 Cholesterol reduction for five treatments.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{trt}{treatment groups, a factor at levels \code{1time}, \code{2times},
               \code{4times}, \code{drugD} and \code{drugE}.}
    \item{response}{cholesterol reduction.}
  }
}
\details{

  A clinical study was conducted to assess the effect of three formulations
  of the same drug on reducing cholesterol. The formulations were
  20mg at once (\code{1time}), 10mg twice a day (\code{2times}), and 5mg
  four times a day (\code{4times}). In addition, two competing drugs
  were used as control group (\code{drugD} and \code{drugE}). The purpose of 
  the study was to find which of the formulations, if any, is efficacious and how
  these formulations compare with the existing drugs. The data were
  published by \bibcitet{multcomp::westfall1999}, page 153.

}
\references{\bibshow{*}}
\examples{

  ### adjusted p-values for all-pairwise comparisons in a one-way layout 
  ### set up ANOVA model  
  amod <- aov(response ~ trt, data = cholesterol)

  ### set up multiple comparisons object for all-pair comparisons
  cht <- glht(amod, linfct = mcp(trt = "Tukey"))

  ### cf. Westfall et al. (1999, page 171)
  summary(cht, test = univariate())
  summary(cht, test = adjusted("Shaffer"))
  summary(cht, test = adjusted("Westfall"))

  ### use only a subset of all pairwise hypotheses
  K <- contrMat(table(cholesterol$trt), type="Tukey")
  Ksub <- rbind(K[c(1,2,5),],
                "D - test" = c(-1, -1, -1, 3, 0),
                "E - test" = c(-1, -1, -1, 0, 3))

  ### reproduce results in Westfall et al. (1999, page 172)
  ### note: the ordering of our estimates here is different
  amod <- aov(response ~ trt - 1, data = cholesterol)
  summary(glht(amod, linfct = mcp(trt = Ksub[,5:1])), 
          test = adjusted("Westfall"))
}
\keyword{datasets}
