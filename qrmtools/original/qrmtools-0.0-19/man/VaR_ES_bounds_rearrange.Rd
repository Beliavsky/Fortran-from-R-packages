\name{VaR_ES_bounds_rearrange}
\alias{rearrange}
\alias{block_rearrange}
\alias{RA}
\alias{ARA}
\alias{ABRA}
\title{Worst and Best Value-at-Risk and Best Expected Shortfall
       for Given Marginals via Rearrangements}
\description{
  Compute the worst and best Value-at-Risk (VaR) and the best expected
  shortfall (ES) for given marginal distributions via rearrangements.
}
\usage{
## Workhorses
## Column rearrangements
rearrange(X, tol = 0, tol.type = c("relative", "absolute"),
          n.lookback = ncol(X), max.ra = Inf,
          method = c("worst.VaR", "best.VaR", "best.ES"),
	  sample = TRUE, is.sorted = FALSE, trace = FALSE, ...)
## Block rearrangements
block_rearrange(X, tol = 0, tol.type = c("absolute", "relative"),
                n.lookback = ncol(X), max.ra = Inf,
                method = c("worst.VaR", "best.VaR", "best.ES"),
                sample = TRUE, trace = FALSE, ...)

## User interfaces
## Rearrangement Algorithm
RA(level, qF, N, abstol = 0, n.lookback = length(qF), max.ra = Inf,
   method = c("worst.VaR", "best.VaR", "best.ES"), sample = TRUE)
## Adaptive Rearrangement Algorithm
ARA(level, qF, N.exp = seq(8, 19, by = 1), reltol = c(0, 0.01),
    n.lookback = length(qF), max.ra = 10*length(qF),
    method = c("worst.VaR", "best.VaR", "best.ES"),
    sample = TRUE)
## Adaptive Block Rearrangement Algorithm
ABRA(level, qF, N.exp = seq(8, 19, by = 1), absreltol = c(0, 0.01),
     n.lookback = NULL, max.ra = Inf,
     method = c("worst.VaR", "best.VaR", "best.ES"),
     sample = TRUE)
}
\arguments{
  \item{X}{(\code{N}, \code{d})-matrix of quantiles (to be
    rearranged). If \code{is.sorted} it is assumed that the columns of
    \code{X} are sorted in \emph{increasing} order.}
  \item{tol}{(absolute or relative) tolerance to determine (the
	individual) convergence.  This should normally be a number
	greater than or equal to 0, but \code{rearrange()} also allows
	for \code{tol = NULL} which means that columns are rearranged
	until each column is oppositely ordered to the sum of all other
	columns.}
  \item{tol.type}{\code{\link{character}} string indicating the
    type of convergence tolerance function to be used (\code{"relative"}
    for relative tolerance and \code{"absolute"} for absolute tolerance).}
  \item{n.lookback}{number of rearrangements to look back for deciding
    about numerical convergence. Use this option with care.}
  \item{max.ra}{maximal number of (considered) column rearrangements
    of the underlying matrix of quantiles (can be set to \code{Inf}).}
  \item{method}{\code{\link{character}} string indicating whether bounds
    for the worst/best VaR or the best ES should be computed.
    These bounds are termed \eqn{\underline{s}_N} and \eqn{\overline{s}_N}
    in the literature (and below) and are theoretically not guaranteed
    bounds of worst/best VaR or best ES;
    however, they are treated as such in practice and are typically in
    line with results from \code{\link{VaR_bounds_hom}()} in the homogeneous
    case, for example.}
  \item{sample}{\code{\link{logical}} indicating whether each column of
    the two underlying matrices of quantiles (see Step 3 of the Rearrangement
    Algorithm in Embrechts et al. (2013))
    are randomly permuted before the rearrangements begin. This typically has
    quite a positive effect on run time (as most of the time is spent
    (oppositely) ordering columns (for \code{rearrange()}) or blocks
    (for \code{block_rearrange()})).}
  \item{is.sorted}{\code{\link{logical}} indicating whether the columns of
    \code{X} are sorted in increasing order.}
  \item{trace}{\code{\link{logical}} indicating whether the underlying
    matrix is printed after each rearrangement step. See
    \code{vignette("VaR_bounds", package = "qrmtools")} for how to interpret
    the output.}
  \item{level}{confidence level \eqn{\alpha}{alpha} for VaR and ES (e.g., 0.99).}
  \item{qF}{\code{d}-list containing the marginal quantile functions.}
  \item{N}{number of discretization points.}
  \item{abstol}{absolute convergence tolerance \eqn{\epsilon}{epsilon}
    to determine the individual convergence, i.e., the change in the computed
    minimal row sums (for \code{method = "worst.VaR"}) or maximal
    row sums (for \code{method = "best.VaR"}) or expected shortfalls (for
    \code{method = "best.ES"}) for the lower bound \eqn{\underline{s}_N}
    and the upper bound \eqn{\overline{s}_N}. \code{abstol} is typically
    \eqn{\ge0}{>= 0}; it can also be \code{\link{NULL}}, see \code{tol}
    above.}
  \item{N.exp}{exponents of the number of discretization points
    (a \code{\link{vector}}) over which the algorithm iterates to find
    the smallest number of discretization points for which the desired
    accuracy (specified by \code{abstol} and \code{reltol}) is attained;
    for each number of discretization points, at most \code{max.ra}-many
    column rearrangements are of the underlying matrix of quantiles
    are considered.}
  \item{reltol}{\code{\link{vector}} of length two containing the
    individual (first component; used to determine convergence of the
    minimal row sums (for \code{method = "worst.VaR"}) or maximal row
    sums (for \code{method = "best.VaR"}) or expected shortfalls (for
    \code{method = "best.ES"}) for \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N}) and the joint (second component; relative
    tolerance between the computed \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N} with respect to \eqn{\overline{s}_N}) relative
    convergence tolerances. \code{reltol} can also be of length one
    in which case it denotes the joint relative tolerance; the
    individual relative tolerance is taken as \code{NULL}
    (see \code{tol} above) in this case.}
  \item{absreltol}{\code{\link{vector}} of length two containing the
    individual (first component; used to determine convergence of the
    minimal row sums (for \code{method = "worst.VaR"}) or maximal row
    sums (for \code{method = "best.VaR"}) or expected shortfalls (for
    \code{method = "best.ES"}) for \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N}) absolute and the joint (second component; relative
    tolerance between the computed \eqn{\underline{s}_N} and
    \eqn{\overline{s}_N} with respect to \eqn{\overline{s}_N}) relative
    convergence tolerances. \code{absreltol} can also be of length one
    in which case it denotes the joint relative tolerance; the
    individual absolute tolerance is taken as 0 in this case.}
  \item{\dots}{additional arguments passed to the underlying
    optimization function. Currently, this is only used if
    \code{method = "best.ES"} in which case the required confidence
    level \eqn{\alpha}{alpha} must be provided as argument \code{level}.}
}
\value{
  \code{rearrange()} and \code{block_rearrange()} return a
  \code{\link{list}} containing
  \describe{
    \item{\code{bound}:}{computed \eqn{\underline{s}_N}
      or \eqn{\overline{s}_N}.}
    \item{\code{tol}:}{reached tolerance (i.e., the (absolute or
      relative) change of the minimal row sum (for
      \code{method = "worst.VaR"}) or maximal row sum
      (for \code{method = "best.VaR"}) or expected shortfall (for
      \code{method = "best.ES"}) after the last rearrangement).}
    \item{\code{converged}:}{\code{\link{logical}} indicating whether
      the desired (absolute or relative) tolerance \code{tol} has been
      reached.}
    \item{\code{opt.row.sums}:}{\code{\link{vector}} containing the
      computed optima (minima for \code{method = "worst.VaR"}; maxima
      for \code{method = "best.VaR"}; expected shortfalls for
      \code{method = "best.ES"}) for the row sums after each (considered)
      rearrangement.}
    \item{\code{X.rearranged}:}{(\code{N}, \code{d})-\code{\link{matrix}}
      containing the rearranged \code{X}.}
    \item{\code{X.rearranged.opt.row}:}{\code{\link{vector}} containing
      the row of \code{X.rearranged} which leads to the final optimal
      sum. If there is more than one such row, the columnwise averaged
      row is returned.}
  }

  \code{RA()} returns a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{bivariate vector containing the computed
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} (the so-called
      rearrangement range) which are typically treated as bounds for
      worst/best VaR or best ES; see also above.}
    \item{\code{rel.ra.gap}:}{reached relative tolerance (also known as
      relative rearrangement gap) between
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N} computed with
      respect to \eqn{\overline{s}_N}.}
    \item{\code{ind.abs.tol}:}{bivariate \code{\link{vector}} containing
      the reached individual absolute tolerances (i.e., the absolute change
      of the minimal row sums
      (for \code{method = "worst.VaR"}) or maximal row sums
      (for \code{method = "best.VaR"}) or expected shortfalls
      (for \code{mehtod = "best.ES"}) for computing \eqn{\underline{s}_N}
      and \eqn{\overline{s}_N};
      see also \code{tol} returned by \code{rearrange()} above).}
    \item{\code{converged}:}{bivariate \code{\link{logical}} vector
      indicating convergence of the computed \eqn{\underline{s}_N} and
      \eqn{\overline{s}_N} (i.e., whether the desired tolerances were
      reached).}
    \item{\code{num.ra}:}{bivariate vector containing the number
      of column rearrangments of the underlying matrices
      of quantiles for
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{opt.row.sums}:}{\code{\link{list}} of length two containing
      the computed optima (minima for \code{method = "worst.VaR"}; maxima
      for \code{method = "best.VaR"}; expected shortfalls for
      \code{method = "best.ES"}) for the row sums after each
      (considered) column rearrangement for the computed
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}; see also
      \code{rearrange()}.}
    \item{\code{X}:}{initially constructed (\code{N}, \code{d})-matrices
      of quantiles for computing
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{X.rearranged}:}{rearranged matrices \code{X} for
      \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{X.rearranged.opt.row}:}{rows corresponding to optimal
      row sum (see \code{X.rearranged.opt.row} as returned by
      \code{rearrange()}) for \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
  }

  \code{ARA()} and \code{ABRA()} return a \code{\link{list}} containing
  \describe{
    \item{\code{bounds}:}{see \code{RA()}.}
    \item{\code{rel.ra.gap}:}{see \code{RA()}.}
    \item{\code{tol}:}{trivariate \code{\link{vector}} containing
      the reached individual (relative for \code{ARA()}; absolute for
      \code{ABRA()}) tolerances and the reached joint
      relative tolerance (computed with respect to \eqn{\overline{s}_N}).}
    \item{\code{converged}:}{trivariate \code{\link{logical}}
      \code{\link{vector}} indicating individual convergence of the computed
      \eqn{\underline{s}_N} (first entry) and \eqn{\overline{s}_N}
      (second entry) and indicating joint convergence of the two bounds
      according to the attained joint relative tolerance (third entry).}
    \item{\code{N.used}:}{actual \code{N} used for computing
      the (final) \eqn{\underline{s}_N} and \eqn{\overline{s}_N}.}
    \item{\code{num.ra}:}{see \code{RA()}; computed for \code{N.used}.}
    \item{\code{opt.row.sums}:}{see \code{RA()}; computed for \code{N.used}.}
    \item{\code{X}:}{see \code{RA()}; computed for \code{N.used}.}
    \item{\code{X.rearranged}:}{see \code{RA()}; computed for \code{N.used}.}
    \item{\code{X.rearranged.opt.row}:}{see \code{RA()}; computed for \code{N.used}.}
  }
}
\details{
  \code{rearrange()} is an auxiliary
  function (workhorse). It is called by \code{RA()} and \code{ARA()}.
  After a column rearrangement of \code{X}, the tolerance between the
  minimal row sum (for the worst VaR) or maximal row sum (for the best
  VaR) or expected shortfall (obtained from the row sums; for the best
  ES) after this rearrangement and the one of \eqn{n.lookback}
  rearrangement steps before is computed and convergence determined.
  For performance reasons, no input checking is done for
  \code{rearrange()} and it can change in future versions to (futher)
  improve run time. Overall it should only be used by experts.

  \code{block_rearrange()}, the workhorse underlying \code{ABRA()},
  is similar to \code{rearrange()} in that it
  checks whether convergence has occurred after every rearrangement by
  comparing the change to the row sum variance from \code{n.lookback}
  rearrangement steps back. \code{block_rearrange()} differs from
  \code{rearrange} in the following ways. First, instead of single columns,
  whole (randomly chosen) blocks (two at a time) are chosen and
  oppositely ordered. Since some of the ideas for improving the speed of
  \code{rearrange()} do not carry over to \code{block_rearrange()}, the
  latter should in general not be as fast as the former.
  Second, instead of using minimal or maximal row
  sums or expected shortfall to determine numerical convergence,
  \code{block_rearrange()} uses the variance of the vector of row sums
  to determine numerical convergence. By default, it targets a variance
  of 0 (which is also why the default \code{tol.type} is \code{"absolute"}).

  For the Rearrangement Algorithm \code{RA()}, convergence of
  \eqn{\underline{s}_N} and \eqn{\overline{s}_N} is determined if the
  minimal row sum (for the worst VaR) or maximal row sum (for the best
  VaR) or expected shortfall (obtained from the row sums; for the best ES)
  satisfies the specified \code{abstol} (so \eqn{\le\epsilon}{<= eps})
  after at most \code{max.ra}-many column rearrangements. This is different
  from Embrechts et al. (2013) who use \eqn{<\epsilon}{< eps} and
  only check for convergence after an iteration through all
  columns of the underlying matrix of quantiles has been completed.

  For the Adaptive Rearrangement Algorithm \code{ARA()}
  and the Adaptive Block Rearrangement Algorithm \code{ABRA()},
  convergence of \eqn{\underline{s}_N} and \eqn{\overline{s}_N}
  is determined if, after at most \code{max.ra}-many column
  rearrangements, the (the individual relative tolerance)
  \code{reltol[1]} is satisfied \emph{and} the
  relative (joint) tolerance between both bounds is at most \code{reltol[2]}.

  Note that \code{RA()}, \code{ARA()} and \code{ABRA()} need to evalute the
  0-quantile (for the lower bound for the best VaR) and
  the 1-quantile (for the upper bound for the
  worst VaR). As the algorithms, due to performance reasons, can only
  handle finite values, the 0-quantile and the 1-quantile need to be
  adjusted if infinite. Instead of the 0-quantile,
  the \eqn{\alpha/(2N)}{alpha/(2N)}-quantile is
  computed and instead of the 1-quantile the
  \eqn{\alpha+(1-\alpha)(1-1/(2N))}{alpha+(1-alpha)(1-1/(2N))}-quantile
  is computed for such margins (if the 0-quantile or the 1-quantile is
  finite, no adjustment is made).

  \code{rearrange()}, \code{block_rearrange()}, \code{RA()}, \code{ARA()}
  and \code{ABRA()} compute \eqn{\underline{s}_N} and
  \eqn{\overline{s}_N} which are, from a practical
  point of view, treated as bounds for the worst (i.e., largest) or the
  best (i.e., smallest) VaR or the best (i.e., smallest ES), but which are
  not known to be such bounds from a theoretical point of view; see also above.
  Calling them \dQuote{bounds} for worst/best VaR or best ES is thus
  theoretically not correct (unless proven) but \dQuote{practical}.
  The literature thus speaks of \eqn{(\underline{s}_N, \overline{s}_N)} as
  the rearrangement gap.
}
\author{Marius Hofert}
\references{
  Embrechts, P., Puccetti, G., \enc{Rüschendorf}{Rueschendorf}, L.,
  Wang, R. and Beleraj, A. (2014).  An Academic Response to Basel
  3.5. \emph{Risks} \bold{2}(1), 25--48.

  Embrechts, P., Puccetti, G. and \enc{Rüschendorf}{Rueschendorf}, L. (2013).
  Model uncertainty and VaR aggregation. \emph{Journal of Banking &
    Finance} \bold{37}, 2750--2764.

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Hofert, M., Memartoluie, A., Saunders, D. and Wirjanto, T. (2017).
  Improved Algorithms for Computing Worst
  Value-at-Risk. \emph{Statistics & Risk Modeling}
  or, for an earlier version, \url{https://arxiv.org/abs/1505.02281}.

  Bernard, C., \enc{Rüschendorf}{Rueschendorf}, L. and Vanduffel,
  S. (2013). Value-at-Risk bounds with variance constraints.
  See \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2342068}.

  Bernard, C. and McLeish, D. (2014). Algorithms for Finding Copulas
  Minimizing Convex Functions of Sums.
  See \url{https://arxiv.org/abs/1502.02130v3}.
}
\seealso{
  \code{\link{VaR_bounds_hom}()} for an ``analytical'' approach for
  computing best and worst Value-at-Risk in the homogeneous casse.

  \code{vignette("VaR_bounds", package = "qrmtools")}
  for more example calls, numerical challenges
  encoutered and a comparison of the different methods for computing
  the worst (i.e., largest) Value-at-Risk.
}
\examples{
### 1 Reproducing selected examples of McNeil et al. (2015; Table 8.1) #########

## Setup
alpha <- 0.95
d <- 8
theta <- 3
qF <- rep(list(function(p) qPar(p, shape = theta)), d)

## Worst VaR
N <- 5e4
set.seed(271)
system.time(RA.worst.VaR <- RA(alpha, qF = qF, N = N, method = "worst.VaR"))
RA.worst.VaR$bounds
stopifnot(RA.worst.VaR$converged,
          all.equal(RA.worst.VaR$bounds[["low"]],
                    RA.worst.VaR$bounds[["up"]], tol = 1e-4))

## Best VaR
N <- 5e4
set.seed(271)
system.time(RA.best.VaR <- RA(alpha, qF = qF, N = N, method = "best.VaR"))
RA.best.VaR$bounds
stopifnot(RA.best.VaR$converged,
          all.equal(RA.best.VaR$bounds[["low"]],
                    RA.best.VaR$bounds[["up"]], tol = 1e-4))

## Best ES
N <- 5e4 # actually, we need a (much larger) N here (but that's time consuming)
set.seed(271)
system.time(RA.best.ES <- RA(alpha, qF = qF, N = N, method = "best.ES"))
RA.best.ES$bounds
stopifnot(RA.best.ES$converged,
          all.equal(RA.best.ES$bounds[["low"]],
                    RA.best.ES$bounds[["up"]], tol = 5e-1))


### 2 More Pareto examples (d = 2, d = 8; hom./inhom. case; explicit/RA/ARA) ###

alpha <- 0.99 # VaR confidence level
th <- 2 # Pareto parameter theta
qF <- function(p, theta = th) qPar(p, shape = theta) # Pareto quantile function
pF <- function(q, theta = th) pPar(q, shape = theta) # Pareto distribution function


### 2.1 The case d = 2 #########################################################

d <- 2 # dimension

## ``Analytical''
VaRbounds <- VaR_bounds_hom(alpha, d = d, qF = qF) # (best VaR, worst VaR)

## Adaptive Rearrangement Algorithm (ARA)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, qF = rep(list(qF), d), method = "best.VaR")
ARAworst <- ARA(alpha, qF = rep(list(qF), d))

## Rearrangement Algorithm (RA) with N as in ARA()
RAbest  <- RA(alpha, qF = rep(list(qF), d), N = ARAbest$N.used, method = "best.VaR")
RAworst <- RA(alpha, qF = rep(list(qF), d), N = ARAworst$N.used)

## Compare
stopifnot(all.equal(c(ARAbest$bounds[1], ARAbest$bounds[2],
                       RAbest$bounds[1],  RAbest$bounds[2]),
                    rep(VaRbounds[1], 4), tolerance = 0.004, check.names = FALSE))
stopifnot(all.equal(c(ARAworst$bounds[1], ARAworst$bounds[2],
                       RAworst$bounds[1],  RAworst$bounds[2]),
                    rep(VaRbounds[2], 4), tolerance = 0.003, check.names = FALSE))


### 2.2 The case d = 8 #########################################################

d <- 8 # dimension

## ``Analytical''
I <- crude_VaR_bounds(alpha, qF = qF, d = d) # crude bound
VaR.W     <- VaR_bounds_hom(alpha, d = d, method = "Wang", qF = qF)
VaR.W.Par <- VaR_bounds_hom(alpha, d = d, method = "Wang.Par", shape = th)
VaR.dual  <- VaR_bounds_hom(alpha, d = d, method = "dual", interval = I, pF = pF)

## Adaptive Rearrangement Algorithm (ARA) (with different relative tolerances)
set.seed(271) # set seed (for reproducibility)
ARAbest  <- ARA(alpha, qF = rep(list(qF), d), reltol = c(0.001, 0.01), method = "best.VaR")
ARAworst <- ARA(alpha, qF = rep(list(qF), d), reltol = c(0.001, 0.01))

## Rearrangement Algorithm (RA) with N as in ARA and abstol (roughly) chosen as in ARA
RAbest  <- RA(alpha, qF = rep(list(qF), d), N = ARAbest$N.used,
              abstol = mean(tail(abs(diff(ARAbest$opt.row.sums$low)), n = 1),
                            tail(abs(diff(ARAbest$opt.row.sums$up)), n = 1)),
              method = "best.VaR")
RAworst <- RA(alpha, qF = rep(list(qF), d), N = ARAworst$N.used,
              abstol = mean(tail(abs(diff(ARAworst$opt.row.sums$low)), n = 1),
                            tail(abs(diff(ARAworst$opt.row.sums$up)), n = 1)))

## Compare
stopifnot(all.equal(c(VaR.W[1], ARAbest$bounds, RAbest$bounds),
                    rep(VaR.W.Par[1],5), tolerance = 0.004, check.names = FALSE))
stopifnot(all.equal(c(VaR.W[2], VaR.dual[2], ARAworst$bounds, RAworst$bounds),
                    rep(VaR.W.Par[2],6), tolerance = 0.003, check.names = FALSE))

## Using (some of) the additional results computed by (A)RA()
xlim <- c(1, max(sapply(RAworst$opt.row.sums, length)))
ylim <- range(RAworst$opt.row.sums)
plot(RAworst$opt.row.sums[[2]], type = "l", xlim = xlim, ylim = ylim,
     xlab = "Number or rearranged columns",
     ylab = paste0("Minimal row sum per rearranged column"),
     main = substitute("Worst VaR minimal row sums ("*alpha==a.*","~d==d.*" and Par("*
                       th.*"))", list(a. = alpha, d. = d, th. = th)))
lines(1:length(RAworst$opt.row.sums[[1]]), RAworst$opt.row.sums[[1]], col = "royalblue3")
legend("bottomright", bty = "n", lty = rep(1,2),
       col = c("black", "royalblue3"), legend = c("upper bound", "lower bound"))
## => One should use ARA() instead of RA()


### 3 "Reproducing" examples from Embrechts et al. (2013) ######################

### 3.1 "Reproducing" Table 1 (but seed and eps are unknown) ###################

## Left-hand side of Table 1
N <- 50
d <- 3
qPar <- rep(list(qF), d)
p <- alpha + (1-alpha)*(0:(N-1))/N # for 'worst' (= largest) VaR
X <- sapply(qPar, function(qF) qF(p))
cbind(X, rowSums(X))

## Right-hand side of Table 1
set.seed(271)
res <- RA(alpha, qF = qPar, N = N)
row.sum <- rowSums(res$X.rearranged$low)
cbind(res$X.rearranged$low, row.sum)[order(row.sum),]


### 3.2 "Reproducing" Table 3 for alpha = 0.99 #################################

## Note: The seed for obtaining the exact results as in Table 3 is unknown
N <- 2e4 # we use a smaller N here to save run time
eps <- 0.1 # absolute tolerance
xi <- c(1.19, 1.17, 1.01, 1.39, 1.23, 1.22, 0.85, 0.98)
beta <- c(774, 254, 233, 412, 107, 243, 314, 124)
qF.lst <- lapply(1:8, function(j){ function(p) qGPD(p, shape = xi[j], scale = beta[j])})
set.seed(271)
res.best <- RA(0.99, qF = qF.lst, N = N, abstol = eps, method = "best.VaR")
print(format(res.best$bounds, scientific = TRUE), quote = FALSE) # close to first value of 1st row
res.worst <- RA(0.99, qF = qF.lst, N = N, abstol = eps)
print(format(res.worst$bounds, scientific = TRUE), quote = FALSE) # close to last value of 1st row


### 4 Further checks ###########################################################

## Calling the workhorses directly
set.seed(271)
ra <- rearrange(X)
bra <- block_rearrange(X)
stopifnot(ra$converged, bra$converged,
          all.equal(ra$bound, bra$bound, tolerance = 6e-3))

## Checking ABRA against ARA
set.seed(271)
ara  <- ARA (alpha, qF = qPar)
abra <- ABRA(alpha, qF = qPar)
stopifnot(ara$converged, abra$converged,
          all.equal(ara$bound[["low"]], abra$bound[["low"]], tolerance = 2e-3),
          all.equal(ara$bound[["up"]],  abra$bound[["up"]],  tolerance = 6e-3))
}
\keyword{programming}