#' The Bernoulli Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Bernoulli distribution with probability equals to \samp{prob}.
#'
#' These functions use the corresponding functions from the \code{\link{binomial}}
#' distribution with argument \samp{size = 1}. Thus, 1 is for success, 0 is for failure.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \samp{length(n) > 1}, the length is taken to be the number required.
#' @param prob vector of probabilities of success of each trial.
#' @param log,log.p logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.
#' @param lower.tail logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.
#' @return \samp{dbern} gives the density, \samp{pbern} gives the distribution function,
#' \samp{qbern} gives the quantile function, and \samp{rbern} generates random deviates.
#' @seealso \code{\link{Binomial}}
#' @keywords distribution
#' @name bernoulli
#' @examples
#' rbern(n = 10, prob = .5)
#' rbern(n = 3, prob = c(0, .5, 1))
#' @export
dbern <- function(x,prob=.5,log=FALSE)
{
  return(dbinom(x, 1, prob=prob, log = log))}

#' @rdname bernoulli
#' @export
pbern <- function(q,prob=.5,lower.tail = TRUE, log.p = FALSE)
{
  return(pbinom(q, 1, prob=prob, lower.tail = lower.tail, log.p = log.p))
  }

#' @rdname bernoulli
#' @export
qbern <- function(p,prob=.5,lower.tail=TRUE, log.p = FALSE)
{
  return(qbinom(p, 1, prob=prob, lower.tail = lower.tail, log.p = log.p))
  }

#' @rdname bernoulli
#' @export
rbern <- function(n,prob=.5)
{
  return(rbinom(n, 1, prob=prob))
  }
