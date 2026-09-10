\name{bestm}
\alias{bestm}
\title{
Function called by makewpstRO to identify which packets are individually
good for correlating with a response
}
\description{
This function is used when you have a huge number of packets where you
want to identify which ones are, individually, candidates for the good
prediction of a response
}
\usage{
bestm(w2mobj, y, percentage = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w2mobj}{
	The w2m object that contains the packets you wish to preselect
}
  \item{y}{
	The response time series
}
  \item{percentage}{
	The percentage of the w2m packets that you wish to select
}
}
\details{
This function naively addresses a very common problem. The object
w2mobj contains a huge number of variables which might shed some light
on the response object \code{y}. The problem is that the dimensionality
of \code{w2mobj} is larger than that of the length of the series \code{y}.

The solution here is to choose a large, but not huge, subset of the variables
that might be potentially useful in correlating with \code{y}, discard the
rest, and return the "best" or preselected variables. Then the dimensionality
is reduced and more sophisticated methods can be used to perform better
quality modelling of the response \code{y} on the packets in \code{w2mobj}.
}
\value{
A list of class w2m with the following components:
\item{m}{A matrix containing the select packets (as columns), reordered
	so that the best packets come first}
\item{ixvec}{A vector which indexes the best packets into the original
	supplied matrix}
\item{pktix}{The original wavelet packet indices corresponding to each packet}
\item{level}{As \code{pktix} but for the wavelet packet levels}
\item{nlevelsWT}{The number of resolution levels in the original wavelet
	packet object}
\item{cv}{The ordered correlations}
}
\author{ G P Nason }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{makewpstRO}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
