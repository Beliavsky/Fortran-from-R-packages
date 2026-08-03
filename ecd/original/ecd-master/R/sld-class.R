#' An S4 class to represent the stable lambda distribution
#'
#' The \code{sld} class serves as an object-oriented interface for the stable lambda distribution.
#' The \code{sld} prefix is also used as the namespace for the related analytic formulae
#' derived in stable lambda distribution.
#'
#' @slot call the match.call slot
#' @slot t numeric
#' @slot nu0 numeric
#' @slot theta numeric
#' @slot convo  numeric
#' @slot beta.a  numeric
#' @slot mu  numeric
#' @slot lambda numeric, this is default to 4.
#'
#' @include ecd-class.R
#' @keywords ecld
#'
#' @author Stephen H. Lihn
#'
#' @section Details:
#'   See \code{\link{dsl}} for definition of stable lambda distribution.
#'
#' @seealso
#'   \code{\link{dstablecnt}} for \eqn{N_\alpha(.)},
#'   and \code{\link{dstdlap}} for \eqn{f_L^{(m)}(.)}.
#'
#' @references
#'   For more detail, see Section 8.2 of
#'   Stephen Lihn (2017).
#'   \emph{A Theory of Asset Return and Volatility
#'   under Stable Law and Stable Lambda Distribution}.
#'   SSRN: 3046732, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3046732}.
#'
#'
#' @exportClass sld
setClass("sld",
         representation(call = "call",
                        t = "numeric",
                        nu0 = "numeric",
                        theta = "numeric",
                        convo = "numeric",
                        beta.a = "numeric",
                        mu = "numeric",
                        lambda = "numeric"),
          prototype(call = call("sld"),
                    t = NaN,
                    nu0 = NaN,
                    theta = NaN,
                    convo = NaN,
                    beta.a = NaN,
                    mu = NaN,
                    lambda = 4)
)
### <---------------------------------------------------------------------->
