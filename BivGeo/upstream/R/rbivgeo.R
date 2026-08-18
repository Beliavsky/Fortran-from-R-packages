#' @importFrom stats runif rgeom
#'
#' @name rbivgeo
#' @aliases rbivgeo1 rbivgeo2
#'
#' @title Generates Random Deviates from the Basu-Dhar Bivariate Geometric Distribution
#'
#' @description This function generates random values from the Basu-Dhar bivariate geometric distribution assuming arbitrary parameter values.
#'
#' @author Ricardo P. Oliveira \email{rpuziol.oliveira@gmail.com}
#' @author Jorge Alberto Achcar \email{achcar@fmrp.usp.br}
#'
#' @references
#'
#' Marshall, A. W., & Olkin, I. (1967). A multivariate exponential distribution. \emph{Journal of the American Statistical Association}, \bold{62}, 317, 30-44.
#'
#' Basu, A. P., & Dhar, S. K. (1995). Bivariate geometric distribution. \emph{Journal of Applied Statistical Science}, \bold{2}, 1, 33-44.
#'
#' Li, J., & Dhar, S. K. (2013). Modeling with bivariate geometric distributions. \emph{Communications in Statistics-Theory and Methods}, \bold{42}, 2, 252-266.
#'
#' Achcar, J. A., Davarzani, N., & Souza, R. M. (2016). Basu–Dhar bivariate geometric distribution in the presence of covariates and censored data: a Bayesian approach. \emph{Journal of Applied Statistics}, \bold{43}, 9, 1636-1648.
#'
#' de Oliveira, R. P., & Achcar, J. A. (2018). Basu-Dhar's bivariate geometric distribution in presence of censored data and covariates: some computational aspects. \emph{Electronic Journal of Applied Statistical Analysis}, \bold{11}, 1, 108-136.
#'
#' @param n number of observations. If length(n) \eqn{> 1}, the length is taken to be the number required.
#' @param theta vector (of length 3) containing values of the parameters \eqn{\theta_1, \theta_2} and \eqn{\theta_{3}} of the Basu-Dhar bivariate Geometric distribution. The parameters are restricted to \eqn{0 < \theta_i < 1, i = 1,2} and \eqn{0 < \theta_{3} \le 1}.
#'
#' @return \code{\link[BivGeo]{rbivgeo1}} and \code{\link[BivGeo]{rbivgeo2}} generate random deviates from the Bash-Dhar bivariate geometric distribution. The length of the result is determined by n, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @return Invalid arguments will return an error message.
#'
#' @usage
#'
#' rbivgeo1(n, theta)
#' rbivgeo2(n, theta)
#'
#' @details
#'
#' The conditional distribution of X given Y is given by:
#'
#' If X < Y, then
#' \deqn{P(X = x | Y = y) = \theta_1^{x - 1}(1 - \theta_1)}
#' If X = Y, then
#' \deqn{P(X = x | Y = y) = \frac{\theta_1^{x - 1}(1 - \theta_1 \theta_{3} - \theta_2 \theta_{3} + \theta_1 \theta_2 \theta_{3})}{1 - \theta_2 \theta_{3}} }
#' If X > Y, then
#' \deqn{P(X = x | Y = y) = \frac{\theta_1^{x - 1} \theta_{3}^{x - y}(1 - \theta_{1} \theta_{3}) (1 - \theta_2)}{1 - \theta_2 \theta_{3}}}
#'
#' @examples
#'
#' rbivgeo1(n = 10, theta = c(0.5, 0.5, 0.7))
#' #       [,1] [,2]
#' #  [1,]    2    1
#' #  [2,]    3    1
#' #  [3,]    1    1
#' #  [4,]    1    1
#' #  [5,]    2    2
#' #  [6,]    1    3
#' #  [7,]    2    2
#' #  [8,]    1    1
#' #  [9,]    1    1
#' # [10,]    2    2
#'
#' rbivgeo2(n = 10, theta = c(0.5, 0.5, 0.7))
#' #       [,1] [,2]
#' #  [1,]    1    1
#' #  [2,]    2    1
#' #  [3,]    2    1
#' #  [4,]    4    1
#' #  [5,]    1    1
#' #  [6,]    2    2
#' #  [7,]    3    2
#' #  [8,]    3    1
#' #  [9,]    3    2
#' # [10,]    1    1
#'
#' @source
#'
#' \code{\link[BivGeo]{rbivgeo1}} generates random deviates using the inverse transformation method. Returns a matrix that the first column corresponds to X generated random values and the second column corresponds to Y generated random values.
#'
#' \code{\link[BivGeo]{rbivgeo2}} generates random deviates using the shock model. Returns a matrix that the first column corresponds to X generated random values and the second column corresponds to Y generated random values. See Marshall and Olkin (1967) for more details.
#'
#' @seealso
#'
#' \code{\link[stats]{Geometric}} for the univariate geometric distribution.
#'
#' @rdname rbivgeo
#' @export

rbivgeo1 <- function(n, theta)
{

	theta1 	<- theta[1]
	theta2 	<- theta[2]
	theta3 	<- theta[3]

	if(theta[1] <= 0 | theta[1] >= 1) return('theta1 out of bounds')
	if(theta[2] <= 0 | theta[2] >= 1) return('theta2 out of bounds')
	if(theta[3] <= 0 | theta[3] > 1)  return('theta3 out of bounds')

  # Generating X

  x 		<- 	 rgeom(n, 1 - theta1 * theta3) + 1
  x 		<-   array(x, dim = c(n, 1))

  # Generating Y

  u 		<- 	runif(n)
  u 		<- 	array(u, dim = c(n, 1))
  y 		<- 	array(c(0), dim = c(n, 1))

  for (i in 1:n)
  {
    if(x[i] == 1 && u[i] < 1 - theta2 * theta3 * (1 - theta1)/(1 - theta1 * theta3)) {y[i] = 1}
    else
      if(x[i] == 1)
      {
        for(s in 2:n) if(u[i] >= 1 - (1 - theta1) * (theta2 * theta3)^(s - 1)/(1 - theta1 * theta3) && u[i] < 1 - (1 - theta1) * (theta2 * theta3)^s/(1 - theta1 * theta3)) {y[i] = s}
      }
    else
      for (j in 1:x[i] - 2) if(u[i] >= 1 - theta2^j && u[i] < 1 - theta2^(j + 1)) {y[i] = j + 1}
    else
      if(u[i] >= 1 - theta2^(x[i] - 1) && u[i] < 1 - (theta2^x[i]) * theta3 * (1 - theta1)/(1 - theta1 * theta3)) {y[i] = x[i]}
    else
      if(u[i] >= 1 - (theta2^x[i]) * theta3 * (1 - theta1)/(1 - theta1 * theta3) && u[i] < 1 - (1 - theta1) * (theta2^(x[i] + 1)) * theta3^2/(1 - theta1 * theta3)) {y[i] = x[i] + 1}
    else
      if(u[i] >= 1 - (1 - theta1) * (theta2^(x[i] + 1)) * theta3^2/(1 - theta1 * theta3))
      {
        for(k in 1 + x[i]:n) if(u[i] >= 1 - (1 - theta1) * theta3^(-x[i] + 1) * (theta2 * theta3)^k/(1 - theta1 * theta3) && u[i] < 1 - (1 - theta1) * theta3^(-x[i] + 1) * (theta2 * theta3)^(k+1)/(1 - theta1 * theta3)) {y[i] = k + 1}
      }
  }

  r.val         <- cbind(x,y)
  return(r.val)
}
#' @export

rbivgeo2 <- function(n, theta)
{

  theta1 	<- theta[1]
  theta2 	<- theta[2]
  theta3 	<- theta[3]

  if(theta[1] <= 0 | theta[1] >= 1) return('theta1 out of bounds')
  if(theta[2] <= 0 | theta[2] >= 1) return('theta2 out of bounds')
  if(theta[3] <= 0 | theta[3] > 1)  return('theta3 out of bounds')

  r1    <-   rgeom(n, 1 - theta1) + 1
  r2    <-   rgeom(n, 1 - theta2) + 1
  r3    <-   rgeom(n, 1 - theta3) + 1

  x 		<- 	 pmin(r1,r3)
  x 		<-   array(x, dim = c(n, 1))

  y 		<- 	 pmin(r2,r3)
  y 		<-   array(y, dim = c(n, 1))

  r.val         <- cbind(x,y)
  return(r.val)
}
