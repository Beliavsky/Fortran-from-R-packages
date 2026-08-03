\name{xclara}
\alias{xclara}
\title{Bivariate Data Set with 3 Clusters}
\description{
  An artificial data set consisting of 3000 points in 3 quite well-separated
  clusters.
}
\usage{data(xclara)}
\format{
  A data frame with 3000 observations on 2 numeric variables (named
  \code{V1} and \code{V2}) giving the
  \eqn{x} and \eqn{y} coordinates of the points, respectively.
}
\source{
  Sample data set accompanying the reference below (file
  \file{xclara.dat} in side \file{clus_examples.tar.gz}).
}
\note{
  Our version of the \code{xclara} is slightly more rounded than the one
  from \code{\link{read.table}("xclara.dat")} and the relative
  difference measured by \code{\link{all.equal}} is \code{1.15e-7} for
  \code{V1} and \code{1.17e-7} for \code{V2} which suggests that our
  version has been the result of a \code{\link{options}(digits = 7)}
  formatting.

  Previously (before May 2017), it was claimed the three cluster were
  each of size 1000, which is clearly wrong.  \code{\link{pam}(*, 3)}
  gives cluster sizes of 899, 1149, and 952, which apart from seven
  \dQuote{outliers} (or \dQuote{mislabellings}) correspond to
  observation indices \eqn{\{1:900\}}{1:900}, \eqn{\{901:2050\}}{901:2050}, and
  \eqn{\{2051:3000\}}{2051:3000}, see the example.
}
\references{
  Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996)
  Clustering in an Object-Oriented Environment.
  \emph{Journal of Statistical Software} \bold{1}.
  \doi{10.18637/jss.v001.i04}
}
\keyword{datasets}
\examples{
## Visualization: Assuming groups are defined as {1:1000}, {1001:2000}, {2001:3000}
plot(xclara, cex = 3/4, col = rep(1:3, each=1000))
p.ID <- c(78, 1411, 2535) ## PAM's medoid indices  == pam(xclara, 3)$id.med
text(xclara[p.ID,], labels = 1:3, cex=2, col=1:3)
\donttest{%% TODO: a clara() call with the _identical_ clustering (but faster!)
 px <- pam(xclara, 3) ## takes ~2 seconds
 cxcl <- px$clustering ; iCl <- split(seq_along(cxcl), cxcl)
 boxplot(iCl, range = 0.7, horizontal=TRUE,
         main = "Indices of the 3 clusters of  pam(xclara, 3)")

 ## Look more closely now:
 bxCl <- boxplot(iCl, range = 0.7, plot=FALSE)
 ## We see 3 + 2 + 2 = 7  clear "outlier"s  or "wrong group" observations:
 with(bxCl, rbind(out, group))
 ## out   1038 1451 1610   30  327  562  770
 ## group    1    1    1    2    2    3    3
 ## Apart from these, what are the robust ranges of indices? -- Robust range:
 t(iR <- bxCl$stats[c(1,5),])
 ##    1  900
 ##  901 2050
 ## 2051 3000
 gc <- adjustcolor("gray20",1/2)
 abline(v = iR, col = gc, lty=3)
 axis(3, at = c(0, iR[2,]), padj = 1.2, col=gc, col.axis=gc)
}% dont
}