#' @importFrom stats runif rgeom
#'
#' @name cfbivgeo
#' @aliases cfbivgeo
#'
#' @title Cross-factorial Moment for the Basu-Dhar Bivariate Geometric Distribution
#'
#' @description This function computes the cross-factorial moment for the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values.
#'
#' @author Ricardo P. Oliveira \email{rpuziol.oliveira@gmail.com}
#' @author Jorge Alberto Achcar \email{achcar@fmrp.usp.br}
#'
#' @references
#'
#' Basu, A. P., & Dhar, S. K. (1995). Bivariate geometric distribution. \emph{Journal of Applied Statistical Science}, \bold{2}, 1, 33-44.
#'
#' Li, J., & Dhar, S. K. (2013). Modeling with bivariate geometric distributions. \emph{Communications in Statistics-Theory and Methods}, 42, \bold{2}, 252-266.
#'
#' Achcar, J. A., Davarzani, N., & Souza, R. M. (2016). Basu–Dhar bivariate geometric distribution in the presence of covariates and censored data: a Bayesian approach. \emph{Journal of Applied Statistics}, \bold{43}, 9, 1636-1648.
#'
#' de Oliveira, R. P., & Achcar, J. A. (2018). Basu-Dhar's bivariate geometric distribution in presence of censored data and covariates: some computational aspects. \emph{Electronic Journal of Applied Statistical Analysis}, \bold{11}, 1, 108-136.
#'
#' de Oliveira, R. P., Achcar, J. A., Peralta, D., & Mazucheli, J. (2018). Discrete and continuous bivariate lifetime models in presence of cure rate: a comparative study under Bayesian approach. \emph{Journal of Applied Statistics}, 1-19.
#'
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{3}} of the Basu-Dhar bivariate Geometric distribution. For real data applications, use the maximum likelihood estimates or Bayesian estimates to get the cross-factorial moment.
#'
#' @return \code{\link[BivGeo]{cfbivgeo}} computes the cross-factorial moment for the Basu-Dhar bivariate geometric distribution for arbitrary parameter values.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' cfbivgeo(theta)
#'
#' @details
#'
#' The cross-factorial moment between X and Y, assuming the Basu-Dhar bivariate geometric distribution, is given by,
#'
#' \deqn{E[XY] = \frac{1 - \theta_1 \theta_2 \theta_{3}^2}{(1 - \theta_1\theta_3)(1 - \theta_2\theta_3)(1 - \theta_1 \theta_2 \theta_{3})}}
#'
#' Note that the cross-factorial moment is always positive.
#'
#' @examples
#'
#' cfbivgeo(theta = c(0.5, 0.5, 0.7))
#' # [1] 2.517483
#' cfbivgeo(theta = c(0.2, 0.5, 0.7))
#' # [1] 1.829303
#' cfbivgeo(theta = c(0.8, 0.9, 0.1))
#' # [1] 1.277864
#' cfbivgeo(theta = c(0.9, 0.9, 0.9))
#' # [1] 35.15246
#'
#' @source
#'
#' \code{\link[BivGeo]{cfbivgeo}} is calculated directly from the definition.
#'
#' @rdname cfbivgeo
#' @export

cfbivgeo <- function(theta)
{
	if(theta[1] <= 0 | theta[1] >= 1) return('theta1 out of bounds')
	if(theta[2] <= 0 | theta[2] >= 1) return('theta2 out of bounds')
	if(theta[3] <= 0 | theta[3] > 1)  return('theta12 out of bounds')

	p1 	<- 1 -  theta[1] * theta[2] * theta[3]^2
	p2 	<- (1 - theta[1] * theta[3]) * (1 - theta[2] * theta[3]) * (1 - theta[1] * theta[2] * theta[3])

	cf <- p1/p2

	return(cf)
}
