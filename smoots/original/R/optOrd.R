#' Optimal Order Selection
#'
#' From a matrix with values of an information criterion for different orders
#' \eqn{p} and \eqn{q} of an autoregressive-moving-average (ARMA) model, the
#' optimal orders are selected.
#'
#'@param mat a numeric matrix, whose rows represent the AR orders
#'\eqn{p = 0, 1, ..., p_{max}}{p = 0, 1, ..., p_max} and whose columns
#'represent the MA orders \eqn{q = 0, 1, q_{max}}{q = 0, 1, q_max}; the
#'elements of the matrix are then the values of an information criterion
#'calculated for ARMA models with the different order combinations; a matrix
#'returned by the function \code{\link{critMatrix}} of the \code{smoots}
#'package shares these characteristics.
#'@param restr a single expression (not a character object) that defines
#'further restrictions; the standard logical operators, e.g. \code{>=},
#'\code{&} or \code{==}, can be used; refer to the rows with \code{p} and to
#'the columns with \code{q}; is set to \code{NULL} by default, i.e. no
#'restrictions are imposed.
#'@param sFUN the selection function; is set to \code{min}, i.e. the minimal
#'value that meets the restrictions \code{restr} is selected and the
#'corresponding orders \eqn{p} and \eqn{q} are returned.
#'
#'@details
#'Given a matrix \code{mat} filled with the values of an information criterion
#'for different estimated ARMA(\eqn{p,q}) models, where the rows represent
#'different orders \eqn{p = 0, 1, ..., p_{max}}{p = 0, 1, ..., p_max} and where
#'the columns represent the orders \eqn{q = 0, 1, ..., q_{max}}{q = 0, 1, ...,
#'q_max}, the function returns a vector with the optimal orders
#'\eqn{p} and \eqn{q}. Further selection restrictions can be passed to the
#'argument \code{restr} as an expression. To implement a restriction, the rows
#'and columns are addressed via \code{p} and \code{q}, respectively. Moreover,
#'standard boolean operators such as \code{==}, \code{>=} or \code{&} can be
#'used. See the Section \emph{Examples} for examples of different restrictions.
#'In many cases, the minimum value of a criterion is considered to indicate
#'the best model. However, in some other cases a different selection approach
#'might be appropriate. Therefore, a selection function can be considered by
#'means of the argument \code{sFUN}. The default is \code{sFUN = min}, i.e. the
#'function \code{\link[base:Extremes]{min}} is applied to select the optimal
#'orders.
#'
#'@md
#'
#'@export
#'
#'@return
#'The function returns a vector with two elements. The first element is the
#'optimal order \eqn{p}, whereas the second element is the selected optimal
#'order \eqn{q}.
#'
#'@author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}
#'
#'@examples
#'\dontrun{
#'set.seed(21)
#'Xt <- arima.sim(model = list(ar = c(1.2, -0.5), ma = 0.7), n = 1000) + 7
#'mat <- smoots::critMatrix(Xt)
#'optOrd(mat)  # without restrictions
#'optOrd(mat, p <= q)  # with one restriction
#'optOrd(mat, p >= 1 & q >= 4)  # with two restrictions
#'}
#'

optOrd <- function(mat, restr = NULL, sFUN = min){
  restr = substitute(restr)
  if(!is.null(restr)){
    p <- row(mat)-1
    q <- col(mat)-1
    ord.opt <- c(which(mat == sFUN(mat[eval(restr)]), arr.ind = TRUE) - 1)
  } else {
    ord.opt <- c(which(mat == sFUN(mat), arr.ind = TRUE) - 1)
  }
  message("The optimal orders are p = ", ord.opt[[1]], " and q = ",
    ord.opt[[2]], ".")
  names(ord.opt) <- c("p", "q")
  ord.opt
}
