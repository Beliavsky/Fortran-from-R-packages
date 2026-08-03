\name{LdaPP}
\alias{LdaPP}
\alias{LdaPP.formula}
\alias{LdaPP.default}
\title{ Robust Linear Discriminant Analysis by Projection Pursuit}
\description{
    Performs robust linear discriminant analysis by the projection-pursuit approach -
    proposed by Pires and Branco (2010) - and returns the results as an object of
    class \code{LdaPP} (aka constructor).
}
\usage{
LdaPP(x, ...)
\method{LdaPP}{formula}(formula, data, subset, na.action, \dots)
\method{LdaPP}{default}(x, grouping, prior = proportions, tol = 1.0e-4,
                 method = c("huber", "mad", "sest", "class"),
                 optim = FALSE,
                 trace=FALSE, ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y~x}, it describes the response
  and the predictors. The formula can be more complicated, such as
  \code{y~log(x)+z} etc (see \code{\link{formula}} for more details).
  The response should
  be a factor representing the response variable, or any vector
  that can be coerced to such (such as a logical variable).}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{method}{method}
  \item{optim}{wheather to perform the approximation using the Nelder and Mead simplex method
  (see function \code{optim()} from package \code{stats}). Default is \code{optim = FALSE}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
    Currently the algorithm is implemented only for binary classification
    and in the following will be assumed that only two groups are present.

  The PP algorithm searches for low-dimensional projections of higher-dimensional
  data where a projection index is maximized. Similar to the original Fisher's proposal
  the squared standardized distance between the observations in the two groups is maximized.
  Instead of the sample univariate mean and standard deviation \code{(T,S)} robust
  alternatives are used. These are selected through the argument \code{method} and can be one of
    \describe{
        \item{huber}{the pair \code{(T,S)} are the robust M-estimates of location and scale}
        \item{mad}{\code{(T,S)} are the Median and the Median Absolute Deviation}
        \item{sest}{the pair \code{(T,S)} are the robust S-estimates of location and scale}
        \item{class}{\code{(T,S)} are the mean and the standard deviation.}
    }

    The first approximation \emph{A1} to the solution is obtained by investigating
    a finite number of candidate directions, the unit vectors defined
    by all pairs of points such that one belongs to the first group
    and the other to the second group. The found solution is stored in the slots
    \code{raw.ldf} and \code{raw.ldfconst}.

    The second approximation \emph{A2} (optional) is performed by
    a numerical optimization algorithm using \emph{A1} as initial solution.
    The Nelder and Mead method implemented in the function \code{\link{optim}} is applied.
    Whether this refinement will be used is controlled by the argument \code{optim}.
    If \code{optim=TRUE} the result of the optimization is stored into the slots
    \code{ldf} and \code{ldfconst}. Otherwise these slots are set equal to
    \code{raw.ldf} and \code{raw.ldfconst}.
}
\value{
  Returns an S4 object of class \code{\link[rrcov]{LdaPP-class}}
}

\references{
    Pires, A. M. and A. Branco, J. (2010)
    Projection-pursuit approach to robust linear discriminant analysis
    \emph{Journal Multivariate Analysis}, Academic Press, Inc., \bold{101}, 2464--2485.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} and
    Ana Pires \email{apires@math.ist.utl.pt}
}

\section{Warning }{Still an experimental version! Only binary classification is supported.}

\seealso{  \code{\link[rrcov]{Linda}}, \code{\link[rrcov]{LdaClassic}} }

\examples{

##
## Function to plot a LDA separation line
##
lda.line <- function(lda, ...)
{
    ab <- lda@ldf[1,] - lda@ldf[2,]
    cc <- lda@ldfconst[1] - lda@ldfconst[2]
    abline(a=-cc/ab[2], b=-ab[1]/ab[2],...)
}

data(pottery)
x <- pottery[,c("MG", "CA")]
grp <- pottery$origin
col <- c(3,4)
gcol <- ifelse(grp == "Attic", col[1], col[2])
gpch <- ifelse(grp == "Attic", 16, 1)

##
## Reproduce Fig. 2. from Pires and branco (2010)
##
plot(CA~MG, data=pottery, col=gcol, pch=gpch)

\dontrun{

ppc <- LdaPP(x, grp, method="class", optim=TRUE)
lda.line(ppc, col=1, lwd=2, lty=1)

pph <- LdaPP(x, grp, method="huber",optim=TRUE)
lda.line(pph, col=3, lty=3)

pps <- LdaPP(x, grp, method="sest", optim=TRUE)
lda.line(pps, col=4, lty=4)

ppm <- LdaPP(x, grp, method="mad", optim=TRUE)
lda.line(ppm, col=5, lty=5)

rlda <- Linda(x, grp, method="mcd")
lda.line(rlda, col=6, lty=1)

fsa <- Linda(x, grp, method="fsa")
lda.line(fsa, col=8, lty=6)

## Use the formula interface:
##
LdaPP(origin~MG+CA, data=pottery)       ## use the same two predictors
LdaPP(origin~., data=pottery)           ## use all predictor variables

##
## Predict method
data(pottery)
fit <- LdaPP(origin~., data = pottery)
predict(fit)
}

}
\keyword{robust}
\keyword{multivariate}
