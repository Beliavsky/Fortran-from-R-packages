#' @importFrom stats runif
#'
#' @name pbivgeo
#' @aliases pbivgeo
#'
#' @title Joint Cumulative Function for the Basu-Dhar Bivariate Geometric Distribution
#'
#' @description This function computes the joint cumulative function of the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values.
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
#' @param x matrix or vector containing the data. If x is a matrix then it is considered as x the first column and y the second column (y argument need be setted to NULL). Additional columns and y are ignored.
#'
#' @param y vector containing the data of y. It is used only if x is also a vector. Vectors x and y should be of equal length.
#'
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_3} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_3 \le 1}.
#'
#' @param lower.tail logical; If TRUE (default), probabilities are \eqn{P(X \le x, Y \le y)} otherwise \eqn{P(X > x, Y > y)}.
#'
#' @return \code{\link[BivGeo]{pbivgeo}} gives the values of the cumulative function.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' pbivgeo(x, y, theta, lower.tail = TRUE)
#'
#' @details
#'
#' The joint cumulative function for a random vector (\eqn{X}, \eqn{Y}) following a Basu-Dhar bivariate geometric distribution could be written as:
#' \deqn{P(X \le x, Y \le y) = 1 - (\theta_{1}\theta_3)^{x} - (\theta_{2}\theta_3)^{y} + \theta_{1}^{x}\theta_{2}^{y} \theta_{3}^{\max(x,y)}}
#' and the joint survival function is given by:
#' \deqn{P(X > x, Y > y) = \theta_{1}^{x}\theta_{2}^{y} \theta_{3}^{\max(x,y)}}
#'
#' @examples
#'
#' # If x and y are integer numbers:
#'
#' pbivgeo(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), lower.tail = TRUE)
#' # [1] 0.79728
#'
#' # If x is a matrix:
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbivgeo(x = matr, y = NULL, theta = c(0.2,0.4,0.7), lower.tail = TRUE)
#' # [1] 0.8424384 0.9787478
#'
#' # If lower.tail = FALSE:
#'
#' pbivgeo(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), lower.tail = FALSE)
#' # [1] 0.01568
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#' pbivgeo(x = matr, y = NULL, theta = c(0.9,0.4,0.7), lower.tail = FALSE)
#' # [1] 0.01975680 0.00139404
#'
#' @source
#'
#' \code{\link[BivGeo]{pbivgeo}} is calculated directly from the definition.
#'
#' @seealso
#'
#' \code{\link[stats]{Geometric}} for the univariate geometric distribution.
#'
#' @rdname pbivgeo
#' @export

pbivgeo <- function(x, y = NULL, theta = c(), lower.tail = TRUE)
{
	if(is.matrix(x))
	{
		x0	<-	x[,1]
		y0	<-	x[,2]
	}
	else if(is.vector(x) & is.vector(y))
	{
		if(length(x)==length(y))
		{
			x0	<-	x
			y0	<-	y
		}
		else
		{
			stop('lengths of x and y are not equal')
		}
	}
	else
	{
		stop('x is not a matrix or x and y are not vectors')
	}

	if(theta[1] <= 0 | theta[1] >= 1) return('theta1 out of bounds')
	if(theta[2] <= 0 | theta[2] >= 1) return('theta2 out of bounds')
	if(theta[3] <= 0 | theta[3] > 1)  return('theta3 out of bounds')

  	# Max values between x and y
  	z 	  <-  pmax(x0, y0)

  	sf    <-  theta[1]^x0 * theta[2]^y0 * theta[3]^z

  	sfx   <-  (theta[1] * theta[3])^x0
  	sfy   <-  (theta[2] * theta[3])^y0

  	cdf   <-  1 - sfx - sfy + sf

  	if(lower.tail) return(cdf) else return(sf)
}
