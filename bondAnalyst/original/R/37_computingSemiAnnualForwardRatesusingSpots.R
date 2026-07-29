#'Calculates Semi-Annual Forward Rates using the given Spot Rates.
#'@details
#'According to information provided by Frank J. Fabozzi (2008), the method \code{saForwards()} is developed to calculate Semi-Annual Forward Rates using the given Spot Rates. Here, \code{spots} is vector of given 1-period (6-month) spot rates, \code{BGN} is six monthly period in which Forward Rate begins (for example  a value of 1 would mean 1-period from now (that is six-months from now) and value of 2 would means two six-monthly periods from now (or 1 years into future) and so on), and \code{APLY} means number of six-monthly periods for which the Forward Rate Applies (for example, a value for 1 means for six-months and APLY value of 4 means for two years so, BGN=2 and APLY=1 means computing implied Forward Rate of 1 year from now for six-months; this is also called a six-months forward rate from one year into future), \code{t} is a vector of number of six-month periods ranging from 1 to any specified number of six-month periods for which the Spot Rates are available, and \code{n} is number of six-month periods under consideration.
#'@param spots A vector.
#'@param BGN A number.
#'@param APLY A number.
#'@param times A vector.
#'@param n A number.
#'@return Input values to five arguments  \code{spots} ,\code{BGN},\code{APLY}, \code{times} and \code{n}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Fabozzi, F. J. (2008). Handbook of Finance: Financial Markets and Instruments. John Wiley & Sons.
#'@examples
#'saForwards(spots=c(0.05,0.054,0.058,0.064,0.070,0.072,0.074,0.078),BGN=1,APLY=1,t=c(1:8),n=8)
#'saForwards(spots=c(0.05,0.054,0.058,0.064,0.070,0.072,0.074,0.078),BGN=2,APLY=1,t=c(1:8),n=8)
#'saForwards(spots=c(0.05,0.054,0.058,0.064,0.070,0.072,0.074,0.078),BGN=6,APLY=1,t=c(1:8),n=8)
#'saForwards(spots=c(0.05,0.054,0.058,0.064,0.070,0.072,0.074,0.078),BGN=2,APLY=4,t=c(1:8),n=8)
#'saForwards(spots=c(0.05,0.054,0.058,0.064,0.070,0.072,0.074,0.078),BGN=4,APLY=2,t=c(1:8),n=8)
#'@export
saForwards<-function (spots,BGN,APLY,times, n)
{
  all = list(spots, times,n )
  if (!is.vector(times) | !is.numeric(times))
    stop("times represents number of payments must be 1 to number of years as a numeric vector.")
  if (any(times <= 0))
    stop("occurance of spotRates cannot have negative values in times.")
  if (length(spots) == 0 | length(times) == 0)
    stop("Insufficient information.")
  if (BGN == 0 | APLY == 0)
    stop("Enter valid Information.")
  if (length(spots) != length(times))
    stop("Number of SpotRates not equal to times.")
  spots <- (spots/2)
  computed_m_a_FR_t_b=((((1+spots[BGN+APLY])^(BGN+APLY))/((1+spots[BGN])^BGN))^(1/APLY))-1
  computed_m_a_FR_t_b <- (computed_m_a_FR_t_b *2)
  (computed_m_a_FR_t_b = round(computed_m_a_FR_t_b, digits=6))
}



