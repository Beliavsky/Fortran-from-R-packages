\name{icc}
\alias{icc}
\title{Intraclass correlation coefficient (ICC) for oneway and twoway models}
\description{
Computes single score or average score ICCs as an index of interrater reliability of quantitative data. Additionally, F-test and confidence interval are computed.
}
\usage{
icc(ratings, model = c("oneway", "twoway"), 
    type = c("consistency", "agreement"), 
    unit = c("single", "average"), r0 = 0, conf.level = 0.95)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{model}{a character string specifying if a '"oneway"' model (default) with row effects random, or a '"twoway"' model with column and row effects random should be applied. You can specify just the initial letter.}
  \item{type}{a character string specifying if '"consistency"' (default) or '"agreement"' between raters should be estimated. If a '"oneway"' model is used, only '"consistency"' could be computed. You can specify just the initial letter.}
  \item{unit}{a character string specifying the unit of analysis: Must be one of '"single"' (default) or '"average"'. You can specify just the initial letter.}
  \item{r0}{specification of the null hypothesis r = r0. Note that a one sided test (H1: r > r0) is performed.}
  \item{conf.level}{confidence level of the interval.}
}
\details{
Missing data are omitted in a listwise way.\cr
When considering which form of ICC is appropriate for an actual set of data, one has take several decisions (Shrout & Fleiss, 1979):\cr\cr
1. Should only the subjects be considered as random effects ('"oneway"' model) or are subjects and raters randomly chosen from a bigger pool of persons ('"twoway"' model).\cr\cr
2. If differences in judges' mean ratings are of interest, interrater '"agreement"' instead of '"consistency"' should be computed.\cr\cr
3. If the unit of analysis is a mean of several ratings, unit should be changed to '"average"'. In most cases, however, single values (unit='"single"') are regarded.
}
\value{
A list with class '"icclist"' containing the following components:
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$model}{a character string describing the selected model for the analysis.}
  \item{$type}{a character string describing the selected type of interrater reliability.}
  \item{$unit}{a character string describing the unit of analysis.}
  \item{$icc.name}{a character string specifying the name of ICC according to McGraw & Wong (1996).}
  \item{$value}{the intraclass correlation coefficient.}
  \item{$r0}{the specified null hypothesis.}
  \item{$Fvalue}{the value of the F-statistic.}
  \item{$df1}{the numerator degrees of freedom.}
  \item{$df2}{the denominator degrees of freedom.}
  \item{$p.value}{the p-value for a two-sided test.}
  \item{$conf.level}{the confidence level for the interval.}
  \item{$lbound}{the lower bound of the confidence interval.}
  \item{$ubound}{the upper bound of the confidence interval.}
}
\references{
Bartko, J.J. (1966). The intraclass correlation coefficient as a measure of reliability. Psychological Reports, 19, 3-11.\cr\cr
McGraw, K.O., & Wong, S.P. (1996), Forming inferences about some intraclass correlation coefficients. Psychological Methods, 1, 30-46.\cr\cr
Shrout, P.E., & Fleiss, J.L. (1979), Intraclass correlation: uses in assessing rater reliability. Psychological Bulletin, 86, 420-428.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{finn}},
  \code{\link[irr]{meancor}},
  \code{\link[irr]{robinson}}
}
\examples{
data(anxiety)
icc(anxiety, model="twoway", type="agreement")

r1 <- round(rnorm(20, 10, 4))
r2 <- round(r1 + 10 + rnorm(20, 0, 2))
r3 <- round(r1 + 20 + rnorm(20, 0, 2))
icc(cbind(r1, r2, r3), "twoway")              # High consistency
icc(cbind(r1, r2, r3), "twoway", "agreement") # Low agreement
}
\keyword{univar}
