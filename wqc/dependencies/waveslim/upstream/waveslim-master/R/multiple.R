#' Rotated Cumulative Variance
#' 
#' Provides the normalized cumulative sums of squares from a sequence of
#' coefficients with the diagonal line removed.
#' 
#' The rotated cumulative variance, when plotted, provides a qualitative way to
#' study the time dependence of the variance of a series.  If the variance is
#' stationary over time, then only small deviations from zero should be
#' present.  If on the other hand the variance is non-stationary, then large
#' departures may exist.  Formal hypothesis testing may be performed based on
#' boundary crossings of Brownian bridge processes.
#' 
#' @param x vector of coefficients to be cumulatively summed (missing values
#' excluded)
#' @return Vector of coefficients that are the sumulative sum of squared input
#' coefficients.
#' @author B. Whitcher
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @export rotcumvar
rotcumvar <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  plus <- 1:n/(n-1) - cumsum(x^2)/sum(x^2)
  minus <- cumsum(x^2)/sum(x^2) - 0:(n-1)/(n-1)
  pmax(abs(plus), abs(minus))
}



#' Testing for Homogeneity of Variance
#' 
#' A recursive algorithm for detecting and locating multiple variance change
#' points in a sequence of random variables with long-range dependence.
#' 
#' For details see Section 9.6 of Percival and Walden (2000) or Section 7.3 in
#' Gencay, Selcuk and Whitcher (2001).
#' 
#' @param x Sequence of observations from a (long memory) time series.
#' @param wf Name of the wavelet filter to use in the decomposition.
#' @param J Specifies the depth of the decomposition.  This must be a number
#' less than or equal to \eqn{\log(\mbox{length}(x),2)}{log(length(x),2)}.
#' @param min.coef Minimum number of wavelet coefficients for testing purposes.
#' Empirical results suggest that 128 is a reasonable number in order to apply
#' asymptotic critical values.
#' @param debug Boolean variable: if set to \code{TRUE}, actions taken by the
#' algorithm are printed to the screen.
#' @return Matrix whose columns include (1) the level of the wavelet transform
#' where the variance change occurs, (2) the value of the test statistic, (3)
#' the DWT coefficient where the change point is located, (4) the MODWT
#' coefficient where the change point is located.  Note, there is currently no
#' checking that the MODWT is contained within the associated support of the
#' DWT coefficient.  This could lead to incorrect estimates of the location of
#' the variance change.
#' @author B. Whitcher
#' @seealso \code{\link{dwt}}, \code{\link{modwt}}, \code{\link{rotcumvar}},
#' \code{\link{mult.loc}}.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @export testing.hov
testing.hov <- function(x, wf, J, min.coef=128, debug=FALSE) {
  n <- length(x)
  change.points <- NULL

  x.dwt <- dwt(x, wf, J)
  x.dwt.bw <- brick.wall(x.dwt, wf, method="dwt") 
  x.modwt <- modwt(x, wf, J)
  x.modwt.bw <- brick.wall(x.modwt, wf)

  for(j in 1:J) {
    cat("##### Level ", j, " #####", fill=TRUE)
    Nj <- n/2^j
    dwt.list <- list(dwt = (x.dwt.bw[[j]])[!is.na(x.dwt.bw[[j]])],
                     left = min((1:Nj)[!is.na(x.dwt.bw[[j]])]) + 1,
                     right = sum(!is.na(x.dwt.bw[[j]])))
    modwt.list <- list(modwt = (x.modwt.bw[[j]])[!is.na(x.modwt.bw[[j]])],
                       left = min((1:n)[!is.na(x.modwt.bw[[j]])]) + 1,
                       right = sum(!is.na(x.modwt.bw[[j]])))
    if(debug) cat("Starting recursion; using", dwt.list$left,
                  "to", dwt.list$right - 1, "...  ")
    change.points <-
      rbind(change.points,
            mult.loc(dwt.list, modwt.list, wf, j, min.coef, debug))
  }
  dimnames(change.points) <-
    list(NULL, c("level", "crit.value", "loc.dwt", "loc.modwt"))
  return(change.points)
}



#' Wavelet-based Testing and Locating for Variance Change Points
#' 
#' This is the major subroutine for \code{\link{testing.hov}}, providing the
#' workhorse algorithm to recursively test and locate multiple variance changes
#' in so-called long memory processes.
#' 
#' For details see Section 9.6 of Percival and Walden (2000) or Section 7.3 in
#' Gencay, Selcuk and Whitcher (2001).
#' 
#' @param dwt.list List of wavelet vector coefficients from the \code{dwt}.
#' @param modwt.list List of wavelet vector coefficients from the \code{modwt}.
#' @param wf Name of the wavelet filter to use in the decomposition.
#' @param level Specifies the depth of the decomposition.
#' @param min.coef Minimum number of wavelet coefficients for testing purposes.
#' @param debug Boolean variable: if set to \code{TRUE}, actions taken by the
#' algorithm are printed to the screen.
#' @return Matrix.
#' @author B. Whitcher
#' @seealso \code{\link{rotcumvar}}, \code{\link{testing.hov}}.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @export mult.loc
mult.loc <- function(dwt.list, modwt.list, wf, level, min.coef, debug)
{
  Nj <- length(dwt.list$dwt)
  N <- length(modwt.list$modwt)
  crit <- 1.358
  change.points <- NULL
  
  if(Nj > min.coef) {
    ## test statistic using the DWT
    P <- cumsum(dwt.list$dwt^2) / sum(dwt.list$dwt^2)
    test.stat <- pmax((1:Nj) / (Nj-1) - P, P - (1:Nj - 1) / (Nj-1))
    loc.dwt <- (1:Nj)[max(test.stat) == test.stat]
    test.stat <- max(test.stat)

    ## location using the MODWT
    P <- cumsum(modwt.list$modwt^2) / sum(modwt.list$modwt^2)
    loc.stat <- pmax((1:N) / (N-1) - P, P - (1:N - 1) / (N-1))
    loc.modwt <- (1:N)[max(loc.stat) == loc.stat]

    if(test.stat > sqrt(2) * crit / sqrt(Nj)) {
      if(debug) cat("Accepted!", fill=TRUE)
      ## Left
      if(debug) cat("Going left; using", dwt.list$left,
                    "to", loc.dwt + dwt.list$left - 1, "...  ")
      temp.dwt.list <- list(dwt = dwt.list$dwt[1:(loc.dwt-1)],
                            left = dwt.list$left,
                            right = loc.dwt + dwt.list$left - 1)
      temp.modwt.list <- list(modwt = modwt.list$modwt[1:(loc.modwt-1)],
                              left = modwt.list$left,
                              right = loc.modwt + modwt.list$left - 1)
      change.points <-
        rbind(c(level, test.stat, loc.dwt + dwt.list$left,
                loc.modwt + modwt.list$left),
              Recall(temp.dwt.list, temp.modwt.list, wf, level, min.coef, debug))
      ## Right
      if(debug) cat("Going right; using", loc.dwt + dwt.list$left + 1,
                    "to", dwt.list$right, "...  ")
      temp.dwt.list <- list(dwt = dwt.list$dwt[(loc.dwt+1):Nj],
                            left = loc.dwt + dwt.list$left + 1,
                            right = dwt.list$right)
      temp.modwt.list <- list(modwt = modwt.list$modwt[(loc.modwt+1):N],
                              left = loc.modwt + modwt.list$left + 1,
                              right = modwt.list$right)
      change.points <-
        rbind(change.points,
              Recall(temp.dwt.list, temp.modwt.list, wf, level, min.coef, debug))
    }
    else
      if(debug) cat("Rejected!", fill=TRUE)
  }
  else
    if(debug) cat("Sample size does not exceed ", min.coef, "!",
                  sep="", fill=TRUE)

  return(change.points)
}

