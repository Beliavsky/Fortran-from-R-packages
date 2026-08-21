\name{gsipairs}
\alias{gsi.pairs}
\alias{gsi.add2pairs}
\alias{gsi.plots}
\title{Internal functions of the compositions package}
\description{
  Creates a paired plot like \code{pairs} but allows to add additional panels
  afterwards
}
\usage{
gsi.pairs(x, labels, panel = points, \dots, main = NULL, oma = NULL, 
    font.main = par("font.main"), cex.main = par("cex.main"), 
    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, cex.labels = NULL, 
    font.labels = 1, row1attop = TRUE, gap = 1, add=list(), 
    xlim=apply(x,2,range), ylim=apply(x,2,range), log="",
    onlyPanel=NULL,noplot=FALSE,trimode=FALSE)
gsi.add2pairs(x,panel,...,noplot=FALSE)
gsi.plots
}
\arguments{
  \item{x}{a multivariate dataset}
  \item{labels}{The names of the variables}
  \item{panel}{The function to performe the actual pairwise plots.}
  \item{\dots}{see \code{graphics::pairs}} %\code{\link[graphics]{pairs}}}
  \item{main}{see \code{graphics::pairs}} %\code{\link[graphics]{pairs}}}
  \item{oma}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{font.main}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{cex.main}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{lower.panel}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{upper.panel}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{diag.panel}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{text.panel}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{label.pos}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{cex.labels}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{font.labels}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{row1attop}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{gap}{see \code{graphics::pairs}} % \code{\link[graphics]{pairs}}}
  \item{add}{additional parameter containing a list of additional
    panels}
  \item{xlim}{additional 2x\code{ncol(x)}-matrix parameter giving in
    \code{xlim[,j]} the xlims of the j-th column
  }
  \item{ylim}{additional 2x\code{ncol()}-matrix parameter giving in
    \code{ylim[,i]} the ylims of the j-th column
  }
  \item{log}{additional parameter with possible values like in
    %\code{\link[graphics]{plot}} 
    \code{graphics::plot}
    allowing to log some plots, without a warning}
  \item{noplot}{Logical indicating wether the plotting should be
    suppressed. This is usefull for plotting single page postscripts.}
}
\details{
  \code{gsi.pairs}  essentially copies the functionality of
  %\code{\link[graphics]{pairs}}. 
  \code{graphics::pairs}.
  However it additionally stores its own
  parameters in the dev.cur() position of gsi.plots and allows to modify
  the parameters and re-do a modified plot afterwards. This is mainly
  done by \code{gsi.add2pairs} by modifying the additional \code{add}
  parameter, that specifies more panels. This mechanism should not be
  used directly, since it is planed to replace it by a more rigorous solution soon.
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}, 
}
\keyword{internal}
