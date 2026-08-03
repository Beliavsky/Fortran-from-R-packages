#' An S4 class to represent the lambda distribution
#' 
#' The \code{ecld} class serves as an object-oriented interface for the lambda distribution. 
#' The \code{ecld} prefix is also used as the namespace for many analytic formulai
#' derived in lambda distribution, especially when lambda = 1,2,3.
#' Because of the extensive use of analytic formulai and enhanced precision through
#' the unit distribution, MPFR is not needed in most cases. This makes option pricing
#' calculation in \code{ecld} much faster than its counterpart built on the more
#' general-purpose \code{ecd} library.
#'
#' @slot call the match.call slot
#' @slot lambda numeric
#' @slot sigma numeric
#' @slot beta  numeric
#' @slot mu  numeric
#' @slot use.mpfr logical, whether to use mpfr for ecld object. If any of the above parameters
#'                is mpfr, then this flag is set to \code{TRUE}.
#' @slot is.sged logical, if \code{TRUE}, interpret parameters as SGED.
#' @slot ecd     the companion object of ecd class (optional)
#' @slot mu_D    the risk-neutral drift, optional, but preferred to have value
#'               if the object is to engage with OGF calculation.
#' @slot epsilon the residual risk, optional as a storage for lambda transformation
#' @slot rho     the momentum shift, optional as a storage for lambda transformation
#' @slot ecd_RN  the risk-neutral companion object of ecd class (optional)
#' @slot status  numeric, bitmap recording the state of the calculation layers.
#'               1: bare bone; 2: ecd; 4: mu_D; 8: ecd_RN
#'
#' @include ecd-class.R
#' @keywords ecld
#' 
#' @author Stephen H. Lihn
#' 
#' @section Details:
#'   The lambda distribution is defined by a depressed polynomial of \eqn{\lambda}-th order, \deqn{
#'     {\left|y(z)\right|}^\lambda + ... - \beta z y(z) = z^2
#'   }
#'   where \eqn{y(z)} must approach minus infinity as \eqn{z} approaches plus or minus infinity.
#'   The density function is defined as \deqn{
#'     P\left(x; \lambda, \sigma, \beta, \mu\right)
#'     \equiv\, \frac{1}{C\,\sigma} 
#'     e^{y\left(\left|\frac{x-\mu}{\sigma}\right|\right)},
#'   }
#'   and \eqn{C} is the normalization constant,\deqn{
#'     C = \int_{-\infty}^{\infty}e^{y(z)}\,dz,
#'   }
#'   where 
#'   \eqn{\lambda} is the shape parameter, 
#'   \eqn{\sigma} is the scale parameter, 
#'   \eqn{\beta} is the asymmetric parameter, 
#'   \eqn{\mu} is the location parameter.
#'   \cr
#'   The distribution is symmetric when \eqn{\beta=0}, which becomes \deqn{
#'     P\left(x; \lambda, \sigma, \mu\right)
#'     \equiv\, \frac{1}{\lambda \Gamma\left(\frac{2}{\lambda}\right) \sigma} 
#'     e^{-{\left|\frac{x-\mu}{\sigma}\right|}^{\frac{2}{\lambda}}}.
#'   }
#'   This functional form is not unfamiliar and has appeared under several names, such as
#'   generalized normal distribution and power exponential distribution, where
#'   \eqn{\lambda < 2}.
#'   \cr
#'   However, we are most interested in \eqn{\lambda >= 2}, which is called the "local regime".
#'   In this regime, the MGF diverges which requires regularization aka truncation of the right tail.
#'   The \eqn{\lambda} option model pays special attention to \eqn{\lambda=2,3,4}
#'   where many closed form solutions can be obtained.
#'   In particular, SPX options fit best at \eqn{\lambda=4}, which is called "quartic lambda".
#'   \cr
#'   Since option model often has to deal with very small numbers which are closed to the machine error 
#'   of double precision calculation, the method supports MPFR. As soon as one of the \code{ecld} parameters
#'   becomes MPFR (by simply multiplying \code{ecd.mp1}), the subsequent calculations will use MPFR. 
#'   
#' @references
#'   For lambda distribution and option pricing model, see 
#'   Stephen Lihn (2015). 
#'   \emph{The Special Elliptic Option Pricing Model and Volatility Smile}.
#'   SSRN: \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2707810}.
#'   \cr
#'   Closed form solutions are derived in 
#'   Stephen Lihn (2016). \emph{Closed Form Solution and Term Structure for SPX Options}.
#'   SSRN: \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2805769}
#'   and \cr 
#'   Stephen Lihn (2017). \emph{From Volatility Smile to Risk Neutral Probability and 
#'   Closed Form Solution of Local Volatility Function}.
#'   SSRN: \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2906522}
#'
#' @exportClass ecld
setClass("ecld",
         representation(call = "call",
                        lambda = "numericMpfr",
                        sigma = "numericMpfr",
                        beta = "numericMpfr",
                        mu = "numericMpfr",
                        use.mpfr = "logical",
                        is.sged = "logical",
                        ecd = "ecd",
                        mu_D = "numericMpfr",
                        epsilon = "numericMpfr",
                        rho = "numericMpfr",
                        ecd_RN = "ecd",
                        status = "numeric"),
          prototype(call = call("ecld"),
                    lambda = NaN,
                    sigma = NaN,
                    beta = NaN,
                    mu = NaN,
                    use.mpfr = FALSE,
                    is.sged = FALSE,
                    ecd = NULL,
                    mu_D = NaN,
                    epsilon = NaN,
                    rho = NaN,
                    ecd_RN = NULL,
                    status = 0)
)
### <---------------------------------------------------------------------->
