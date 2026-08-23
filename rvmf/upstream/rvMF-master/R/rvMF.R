#' von Mises--Fisher Distributed Pseudo-Random Vector Generator
#'
#' `rvMF()` generates von Mises--Fisher distributed pseudo-random vectors,
#' without resorting to the rejection-based sampling method proposed by Wood
#' (1994). See Kang and Oh (2024) for details. This function partly uses the
#' code from the function [Rfast::rvmf()] and the article Marsaglia et al.
#' (2004).
#'
#' @param n number of pseudo-random vectors to generate.
#' @param mu mean direction.
#' @param k concentration parameter. \ifelse{html}{\code{k} \out{&ge; 0}}{\eqn{k\ge 0}}.
#' @returns matrix where each row independently follows the specified von
#'   Mises-Fisher distribution. The number of columns equals the length of `mu`,
#'   and the number of rows equals `n` for `rvMF`.
#' @seealso [rvMFangle()], [dvMFangle()], [Rfast::rvmf()].
#' @examples
#' rvMF(10, c(0,0,1), 10)
#' rvMF(10, c(1,1)/sqrt(2), 0)
#' @references
#' S. Kang and H.-S. Oh. Novel sampling method for the von Mises--Fisher
#' distribution. \emph{Statistics and Computing}, 34(3):106, 2024.
#'
#' K. V. Mardia and P. E. Jupp. \emph{Directional Statistics}, volume 494. John
#' Wiley & Sons, Chichester, 1999.
#'
#' G. Marsaglia, W. W. Tsang, and J. Wang. Fast generation of discrete random
#' variables. \emph{Journal of Statistical Software}, 11(3):1–11, 2004.
#'
#' M. Papadakis, M. Tsagris, M. Dimitriadis, S. Fafalios, I. Tsamardinos, M.
#' Fasiolo, G. Borboudakis, J. Burkardt, C. Zou, K. Lakiotaki, and C.
#' Chatzipantsiou. \emph{Rfast: A Collection of Efficient and Extremely Fast R
#' Functions}, 2022. <https://CRAN.R-project.org/package=Rfast>. R package
#' version 2.0.6.
#'
#' A. T. Wood. Simulation of the von Mises Fisher distribution. \emph{Communications
#' in Statistics– Simulation and Computation}, 23(1):157–164, 1994.
#' @export
rvMF <- function (n, mu, k)
{
  rotation <- function(a, b) {
    p <- length(a)
    ab <- sum(a * b)
    ca <- a - b * ab
    ca <- ca/sqrt(sum(ca^2))
    A <- b %*% t(ca)
    A <- A - t(A)
    theta <- acos(ab)
    diag(p) + sin(theta) * A + (cos(theta) - 1) * (b %*% t(b) + ca %*% t(ca))
  }
  d <- length(mu)
  if (k > 0) {
    mu <- mu/sqrt(sum(mu^2))
    ini <- c(numeric(d - 1), 1)
    d1 <- d - 1
    v1 <- Rfast::matrnorm(n, d1)
    v <- v1/sqrt(Rfast::rowsums(v1^2))
    # w <- .rvMF64(n,d,k,scModels::chf_1F1(2*k,d1/2,d1))
    w <- .rvMF64(n,d,k,.chf(k, d1))
    S <- cbind(sqrt(1 - w^2) * v, w)
    if (isTRUE(all.equal(ini, mu, check.attributes = FALSE))) {
      x <- S
    }
    else if (isTRUE(all.equal(-ini, mu, check.attributes = FALSE))) {
      x <- -S
    }
    else {
      A <- rotation(ini, mu)
      x <- tcrossprod(S, A)
    }
  }
  else {
    x1 <- Rfast::matrnorm(n, d)
    x <- x1/sqrt(Rfast::rowsums(x1^2))
  }
  colnames(x) <- names(mu)
  x
}

#' @name vMFangle
#' @title Inner Product of von Mises--Fisher Random Vector and Mean Direction
#'
#' @description
#' These functions provide information about the distribution of an inner
#' product between von Mises--Fisher random vector and its mean direction.
#' Specifically, if \eqn{X} follows a von Mises--Fisher distribution with mean
#' direction \eqn{\mu}, the inner product \eqn{X'\mu} will be a random variable
#' following some distribution. See page 170 of Mardia and Jupp (1999).
#' `rvMFangle()` generates random variates using the algorithm proposed in Kang
#' and Oh (2024), and `dvMFangle` gives the density from this distribution. This
#' function partly uses the code from the article Marsaglia et al. (2004).
#'
#' @param n number of random vectors to generate.
#' @param r vector of quantiles. \ifelse{html}{\out{-1 &le;} \code{r} \out{&le; 1}}{\eqn{-1\le r\le 1}}.
#' @param p dimension of the sphere. i.e.,
#'   \ifelse{html}{\out{S<sup>p-1</sup>}}{\eqn{S^{p-1}}}, \ifelse{html}{\code{p} \out{&ge; 2}}{\eqn{p\ge 2}}.
#' @param kappa concentration parameter. \code{kappa > 0}.
#'   Setting `kappa = 0` may cause errors.
#' @returns
#' * `rvMFangle()` returns a vector whose components independently follow the
#' aforementioned distribution. The length of the result is determined by `n`
#' for `rvMFangle()`.
#'
#' * `dvMFangle()` returns a vector of density function value. The length of the
#' result is determined by the length of `r` for `dvMFangle()`.
#' @examples
#' rvMFangle(10, 2, 10)
#' rvMFangle(10, 3, 0.1)
#' dvMFangle(seq(-1,1,by=0.01), 2, 10)
#' dvMFangle(seq(0,1,by=0.01), 3, 0.1)
#' @seealso [rvMF()] wrapper of `rvMFangle()`.
#' @references
#' S. Kang and H.-S. Oh. Novel sampling method for the von Mises--Fisher
#' distribution. \emph{Statistics and Computing}, 34(3):106, 2024.
#'
#' K. V. Mardia and P. E. Jupp. \emph{Directional Statistics}, volume 494. John Wiley
#' & Sons, Chichester, 1999.
#'
#' G. Marsaglia, W. W. Tsang, and J. Wang. Fast generation of discrete random
#' variables. \emph{Journal of Statistical Software}, 11(3):1–11, 2004.
#' @export
# rvMFangle <- function(n,p,kappa) .rvMF64(n,p,kappa,scModels::chf_1F1(2*kappa,(p-1)/2,p-1))
rvMFangle <- function(n,p,kappa) .rvMF64(n,p,kappa,.chf(kappa,p-1))

#' @name vMFangle
#' @export
dvMFangle <- function(r,p,kappa){
  nu <- p/2-1
  if(kappa == 0){
    return(exp((2-p)*log(2) + lgamma(p-1) - lgamma((2*nu+1)/2)*2 + (p-3)/2*log(1-r^2)))
  } else{
    if(Re(Bessel::BesselI(kappa,nu)) != Inf & Re(Bessel::BesselI(kappa,nu)) != 0){
      return(exp(nu*log(kappa/2) - 1/2*log(pi) - log(Re(Bessel::BesselI(kappa,nu))) - lgamma((2*nu+1)/2) +
                   (p-3)/2*log(1-r^2) + kappa*r))
    } else{ # more precision
      return(exp(nu*log(kappa/2) - 1/2*log(pi) - Bessel::besselI.nuAsym(kappa,nu, k.max=5, expon.scaled = T, log = T) - kappa -
                   lgamma((2*nu+1)/2) + (p-3)/2*log(1-r^2) + kappa*r))
    }
  }
}
