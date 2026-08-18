#' @importFrom stats runif
#'
#' @name dbivgeo
#' @aliases dbivgeo1 dbivgeo2
#'
#' @title Joint Probability Mass Function for the Basu-Dhar Bivariate Geometric Distribution
#'
#' @description This function computes the joint probability mass function of the Basu-Dhar bivariate geometric distribution for arbitrary parameter values.
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
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{3}} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_{3} \le 1}.
#'
#' @param log logical argument for calculating the log probability or the probability function. The default value is FALSE.
#'
#' @return \code{\link[BivGeo]{dbivgeo1}} gives the values of the probability mass function using the first form of the joint pmf.
#'
#' @return \code{\link[BivGeo]{dbivgeo2}} gives the values of the probability mass function using the second form of the joint pmf.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' dbivgeo1(x, y = NULL, theta, log = FALSE)
#' dbivgeo2(x, y = NULL, theta, log = FALSE)
#'
#' @details
#'
#' The joint probability mass function for a random vector (\eqn{X}, \eqn{Y}) following a Basu-Dhar bivariate geometric distribution could be written in two forms. The first form is described by:
#' \deqn{P(X = x, Y = y) = \theta_{1}^{x - 1} \theta_{2}^{y - 1} \theta_{3}^{z_1} - \theta_{1}^{x} \theta_{2}^{y - 1} \theta_{3}^{z_2} - \theta_{1}^{x - 1} \theta_{2}^{y} \theta_{2}^{z_3} + \theta_{1}^{x} \theta_{2}^{y} \theta_{3}^{z_4}}
#' where \eqn{x,y > 0} are positive integers and \eqn{z_1 = \max(x - 1, y - 1),z_2 = \max(x, y - 1), z_3 = \max(x - 1, y), z_4 = \max(x, y)}. The second form is given by the conditions:
#'
#' If X < Y, then
#' \deqn{P(X = x, Y = y) = \theta_1^{x - 1} (\theta_2 \theta_{3})^{y - 1}(1 - \theta_{2} \theta_{3}) (1 - \theta_1)}
#' If X = Y, then
#' \deqn{P(X = x, Y = y) = (\theta_1 \theta_2 \theta_{3})^{x - 1}(1 - \theta_1 \theta_{3} - \theta_2 \theta_{3} + \theta_1 \theta_2 \theta_{3}) }
#' If X > Y, then
#' \deqn{P(X = x, Y = y) = \theta_2^{y - 1} (\theta_1 \theta_{3})^{x - 1}(1 - \theta_{1} \theta_{3}) (1 - \theta_2)}
#'
#' @examples
#'
#' # If x and y are integer numbers:
#'
#' dbivgeo1(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), log = FALSE)
#' # [1] 0.16128
#' dbivgeo2(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), log = FALSE)
#' # [1] 0.16128
#'
#' # If x is a matrix:
#'
#' matr 	<- 	 matrix(c(1,2,3,5), ncol = 2)
#'
#' dbivgeo1(x = matr, y = NULL, theta = c(0.2,0.4,0.7), log = FALSE)
#' # [1] 0.0451584000 0.0007080837
#' dbivgeo2(x = matr, y = NULL, theta = c(0.2,0.4,0.7), log = FALSE)
#' # [1] 0.0451584000 0.0007080837
#'
#' # If log = TRUE:
#'
#' dbivgeo1(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), log = TRUE)
#' # [1] -1.824613
#' dbivgeo2(x = 1, y = 2, theta = c(0.2, 0.4, 0.7), log = TRUE)
#' # [1] -1.824613
#'
#' @source
#'
#' \code{\link[BivGeo]{dbivgeo1}} and \code{\link[BivGeo]{dbivgeo2}} are calculated directly from the definitions.
#' @seealso
#'
#' \code{\link[stats]{Geometric}} for the univariate geometric distribution.
#'
#' @rdname dbivgeo
#' @export

dbivgeo1 <- function(x, y = NULL, theta = c(), log = FALSE)
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

	# Max values

	z1		<-	pmax(x0 - 1, y0 - 1)
	z2		<-  pmax(x0, y0 - 1)
	z3 		<-  pmax(x0 - 1, y0)
	z4 		<-  pmax(x0, y0)

  	# Joint pmf parts
  	p1	  	<-  theta[1]^(x0 - 1) * theta[2]^(y0 - 1) * theta[3]^z1
  	p2	  	<-  theta[1]^x0 	  * theta[2]^(y0 - 1) * theta[3]^z2
  	p3	  	<-  theta[1]^(x0 - 1) * theta[2]^y0       * theta[3]^z3
  	p4	  	<-  theta[1]^x0 	  * theta[2]^y0       * theta[3]^z4

	# Joint pmf
  	pmf   <-  p1 - p2 - p3 + p4

  	if(log) return(log(pmf)) else return(pmf)
}

#' @export

dbivgeo2 <- function(x, y = NULL, theta = c(), log = FALSE)
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

	# Parameters relations

	gamma1	<- 	theta[1] * theta[3]
	gamma2	<- 	theta[2] * theta[3]
	gamma3	<- 	theta[1] * theta[2] * theta[3]

	# Joint pmf parts

	p1    	<-  theta[1]^(x0 - 1) * (gamma2)^(y0 - 1)
	p2    	<-  (1 - gamma2) * (1 - theta[1])
	P1    	<-  p1 * p2

	p3    	<-  (gamma3)^(x0 - 1)
	p4    	<-  (1 - gamma1 - gamma2 + gamma3)
	P2    	<-  p3 * p4

	p5    	<-  theta[2]^(y0 - 1) * (gamma1)^(x0 - 1)
	p6    	<-  (1 - gamma1) * (1 - theta[2])
	P3  	<-  p5 * p6

  	# Joint pmf

  	pmf 	<- ifelse(x0 < y0, P1, ifelse(x0 > y0, P3, P2))

  	if(log) return(log(pmf)) else return(pmf)
}
