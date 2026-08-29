#' @title Display Integrand Information for Tweedie Fourier inversion
#' @name tweedie_integrand
#' @description Plots the integrand for Fourier inversion and the real and imaginary parts separately.
#'
#' @usage
#' tweedie_integrand(y, power, mu, phi, t = seq(0, 5, length = 200), 
#'                   type = "PDF", whichPlots = 1:4, 
#'                   plot_args  = list())
#'
#' @param y vector of quantiles.
#' @param power a synonym for \eqn{\xi}{xi}; the Tweedie power-index on the variance.
#' @param mu the mean parameter \eqn{\mu}{mu}.
#' @param phi the dispersion parameter \eqn{\phi}{phi}.
#' @param t the values of the variable over which to integrate; the default is \code{t = seq(0, 5, length = 200)}.
#' @param type either \code{"PDF"} (the default) for the (probability) density function, or \code{"CDF"} for the (cumulative) distribution function.
#' @param whichPlots which combination of the four plots (described below) are produced; by default, all four are produced (i.e., \code{whichPlots = 1:4}).
#' @param plot_args A named list of arguments controlling the main plot.
#'   These are passed to \code{\link[graphics]{plot}}.
#'
#' @details The Tweedie family of distributions belong to the class of exponential dispersion models (\acronym{edm}s), famous for their role in generalized linear models. 
#' The Tweedie distributions are the \acronym{edm}s with a variance of the form \eqn{\mbox{var}[Y] = \phi\mu^p}{var(Y) = phi * mu^power} where \eqn{p}{power} is greater than or equal to one, or less than or equal to zero.
#'
#' \emph{This function only evaluates for \eqn{p}{power} greater than or equal to one.}
#'
#' Special cases include the normal (\eqn{p = 0}{power = 0}), Poisson (\eqn{p = 1}{power = 1} with \eqn{\phi = 1}{phi = 1}), gamma (\eqn{p = 2}{power = 2}) and inverse Gaussian (\eqn{p = 3}{power = 3}) distributions. For other values of \code{power}, the distributions are still defined but cannot be written in closed form, and hence evaluation is very difficult.
#'
#' When \eqn{1 < p < 2}{1 < power < 2}, the distribution are continuous for \eqn{Y}{Y} greater than zero, with a positive mass at \eqn{Y = 0}{Y = 0}. For \eqn{p > 2}{power > 2}, the distributions are continuous for \eqn{Y}{Y} greater than zero.
#'
#' This function displays the integrand that is evaluated for computing the Fourier inversion, for the PDF or CDF.
#'
#' @return A list containing the real and imaginary parts of \eqn{k(t)}{k(t)}, \code{Real} and \code{Imag} respectively, plus the values of the integrand as \code{IG}.
#' The main purpose of the function is the side-effect of producing a \eqn{2\times2}{2x2} grid of plots.
#' The first is the imaginary parts of \eqn{k(t)}{k(t)}.
#' The second is \eqn{\cos\Im k(t)}{sinIm[k(t)]} (for \code{type=="PDF"}) or \eqn{\sin\Im k(t)}{sinIm[k(t)]} (for \code{type=="CDF"}).
#' The third is the real part of \eqn{\Re k(t)}{Re[k(t)]}
#' The fourth is the integrand, with the envelope shown as a dashed line.
#'
#' @author Peter Dunn (\email{pdunn2@usc.edu.au})
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (2008). Evaluation of Tweedie exponential dispersion model densities by Fourier inversion. \emph{Statistics and Computing}, \bold{18}, 73--86. \doi{10.1007/s11222-007-9039-6}
#'
#' Dunn, Peter K and Smyth, Gordon K (2005). Series evaluation of Tweedie exponential dispersion model densities \emph{Statistics and Computing}, \bold{15}(4). 267--280. \doi{10.1007/s11222-005-4070-y}
#'
#' Dunn, Peter K and Smyth, Gordon K (2001). Tweedie family densities: methods of evaluation. \emph{Proceedings of the 16th International Workshop on Statistical Modelling}, Odense, Denmark, 2--6 July
#'
#' Jorgensen, B. (1987). Exponential dispersion models. \emph{Journal of the Royal Statistical Society}, B, \bold{49}, 127--162.
#'
#' Jorgensen, B. (1997). \emph{Theory of Dispersion Models}. Chapman and Hall, London.
#'
#' Tweedie, M. C. K. (1984). An index which distinguishes between some important exponential families. \emph{Statistics: Applications and New Directions. Proceedings of the Indian Statistical Institute Golden Jubilee International Conference} (Eds. J. K. Ghosh and J. Roy), pp. 579-604. Calcutta: Indian Statistical Institute.
#'
#' @examples
#' tweedie_integrand(2, power = 3, mu = 1, phi = 1,
#'                   plot_args = list(lwd = 2))
#'
#' @seealso \code{\link{dtweedie}}
#' 
#' @importFrom graphics par mtext abline axis lines
#' @keywords models
#' 
#' @export
tweedie_integrand <- function(y, power, mu, phi, 
                              t = seq(0, 5, length = 200), 
                              type = "PDF", whichPlots = 1:4, 
                              plot_args  = list()) {

  # BEGIN: Define function to be used
  k <- function(p, mu, phi, y, t){
    front <- ( mu^(2 - p) ) / (phi * (2 - p) )
    alpha <- (2 - p)/(1 - p) 
    omega <- atan( ( (1 - p) * t * phi) / mu^(1 - p) )

    Real <- front * ( cos(omega * alpha) / ( cos(omega)^alpha) - 1) 
    Imag <- front *   sin(omega * alpha) / ( cos(omega)^alpha ) - (t * y)

    return( list(Real = Real,
                 Imag = Imag) )
  }
  # END: Define function to be used
  
  
  
  # BEGIN: Check for some given parameters supplied via args
  plot_defaults <- list(
    las = 1,
    xlab = expression(Values~of~italic(t)),
    col = "black", 
    lwd = 2)
  
  plot_args  <- utils::modifyList(plot_defaults, 
                                  plot_args)
  # END: Check for some given parameters supplied via args

  
  kvals <- k(p = power, 
             mu = mu, 
             phi = phi,
             y = y, 
             t = t)
  k_Imag <- kvals$Imag
  k_Real <- kvals$Real
  
  if(type == "PDF" ){
    igrand <- exp(k_Real) * cos(k_Imag)
    envelope <- exp(k_Real)
  } else {
    igrand <- exp(k_Real) * sin(k_Imag)/t
    envelope <- exp(k_Real) / t
  }
  
  

  whichPlots <- unique(whichPlots)
  if ( length(whichPlots) > 2 ) {
    par( mfrow = c(2, 2))
  } else {
    ifelse( length(whichPlots) == 2,
            par( mfrow = c(1, 2)),
            par( mfrow = c(1, 1))
    )
  }

  ###
  ### PLOT 1: Im k(t) vs t
  ###
  if (1 %in% whichPlots) {
    do.call(plot,
            c(list(
              x = t,
              y = k_Imag,
              type = "l",
              main = expression( bold(Imaginary~part~of~italic(k)*(italic(t)))),
              ylab = expression(Im*"("*italic(k)*")")),
              plot_args
            ))
    
    # Determine  m  values to display
    mMax <- max(0, max(kvals$Imag)) / pi + ifelse(type == "CDF", 0, pi/2)
    mMin <- min(0, min(kvals$Imag)) / pi - ifelse(type == "CDF", 0, pi/2)
    mValues <- ( floor(mMin) : floor(mMax) )
  
    # Adornments
    mtext(text = expression(atop(italic(m), " ") ),
          side = 4,
          las = 1,
          cex = 0.9,
          line = 2,
          adj = 1,
          at = max(k_Imag) )
    abline(h = 0, 
           col="grey")
    abline(h = mValues * pi + ifelse(type=="PDF", pi/2, 0),
           lty = 2,
           col="grey")
    axis(side = 4, 
         at = mValues * pi + ifelse(type=="PDF", pi/2, 0),
         las = 1,
         cex = 0.8,
         labels = mValues )
  }
  
  ###
  ### PLOT 2: cos( Im k(t) ) (PDF) or sin( Im k(t) ) (CDF) vs t
  ###
  if (2 %in% whichPlots){
    if (type == "PDF"){
      do.call(plot,
              c(list(
                x = t,
                y = cos(k_Imag),
                type = "l",
                main = expression(cos(Im*"("*italic(k)*")")),
                ylab = expression(cos(Im*"("*italic(k)*")"))),
                plot_args
              ))
    } else {
      do.call(plot,
              c(list(
                x = t,
                y = sin(k_Imag),
                type = "l",
                main = expression(sin(Im*"("*italic(k)*")")),
                ylab = expression(cos(Im*"("*italic(k)*")"))),
                plot_args 
              ))
    }
    abline(h = 0, 
           col="grey")
  }
  
  
  ###
  ### PLOT 3: sin( Re k(t) ) vs t
  ###
  if (3 %in% whichPlots) {
    do.call(plot,
            c(list(
              x = t,
              y = k_Real,
              type = "l",
              main = "Real part of k(t)",
              ylab = expression(Re*"("*italic(k)*")")),
              plot_args 
            ))
    abline(h = 0, 
           col="grey")
    }
  
  
  ###
  ### PLOT 4: Integrand
  ###
  if (4 %in% whichPlots) {  
    do.call(plot,
            c(list(
              x = t,
              y = igrand,
              type = "l",
              main = "Integrand",
              ylab = "Integrand"),
              plot_args 
            ))

    abline(h = 0, 
           col="grey")
    lines(x = t,
          y = envelope,
          lty = 2,
          col = "grey")
    lines(x = t,
          y = -envelope,
          lty = 2,
          col = "grey")
  }
  
  return(invisible( list(Real = k_Real,
                         Imag = k_Imag,
                         IG = igrand)) )
  
}