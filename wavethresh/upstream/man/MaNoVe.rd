\name{MaNoVe}
\alias{MaNoVe}
\title{Make Node Vector (using Coifman-Wickerhauser best-basis type algorithm)}
\usage{
MaNoVe(\dots)
}
\arguments{
\item{\dots}{Methods may have other arguments}
}
\description{
This generic function chooses a ``best-basis'' using the Coifman-Wickerhauser (1992) algorithm.  This function is generic.  Particular methods exist:
\code{\link{MaNoVe.wp}} and \code{\link{MaNoVe.wst}}.
}
\details{
Description says all.
}
\value{
A node vector, which describes a particular basis specification relevant
to the kind of object that the function was applied to.
}
\seealso{
\code{\link{MaNoVe.wp}},
\code{\link{MaNoVe.wst}},
\code{\link{wp.object}},
\code{\link{wst.object}},
\code{\link{wp}},
\code{\link{wst}}}
\author{G P Nason}
\keyword{smooth}
