\name{plot.CLA}
\alias{plot.CLA}
\title{Plotting CLA() results including Efficient Frontier}
\description{
  A partly experimental \code{\link{plot}()} method for
  \code{\link{CLA}()} objects.

  It draws the efficient frontier in the \eqn{\mu_w,\sigma_w} aka (mean,
  std.dev.) plane.

  Currently, this is quite rudimentary.\cr
  Future improvements would allow
  - to add the/some single asset points,
  - to correctly   (`hyperbolically') interpolate between turning points
  - add text about the number of (unique) critical points
  - add option \code{add = FALSE} which when TRUE would use
    \code{\link{lines}} instead \code{plot}.
}
\usage{%% --> ../R/CLA.R
\method{plot}{CLA}(x, type = "o", main = "Efficient Frontier",
     xlab = expression(sigma(w)),
     ylab = expression(mu(w)),
     col = adjustcolor("blue", alpha.f = 0.5),
     pch = 16, \dots)
}
\arguments{
  \item{x}{a named \code{\link{list}} as resulting from \code{\link{CLA}()}.}
  \item{type}{the \code{\link{lines}}/\code{\link{points}} types used
  for the efficient frontier.  This will become more sophisticated,
  i.e., \emph{may change non-compatibly!!}}
  \item{main}{main \code{\link{title}}.}
  \item{xlab, ylab}{x- and y- axis labels, passed to \code{\link{plot.default}}.}
  \item{col, pch}{color and point type, passed to
  \code{\link{plot.default}}, but with differing defaults in this method.}
  \item{\dots}{potentially further arguments passed to
    \code{\link{plot}}, i.e., \code{\link{plot.default}}.}
}
%% \details{
%% }
%% \value{
%% }
\author{Martin Maechler.}
\seealso{
  \code{\link{CLA}}, \code{\link{plot.default}}.
}
\examples{
## TODO %% Add A. Norring's small 12-asset example see --> ../TODO
## ---- one example is in help(CLA)
}
\keyword{hplot}

