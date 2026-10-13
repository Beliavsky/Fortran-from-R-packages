

#' Call Premium
#'
#' Calculates the premium of a European-style call option using the Black-Scholes option pricing model
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the value of the call option
#' @export
#'
#'@importFrom stats pnorm
#'
#' @examples callpremium(100, 100, 0.20, (45/365), 0.02, 0.02)
callpremium <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  s*exp(-d*t)*pnorm(d1)-x*exp(-r*t)*pnorm(d2)
}


#' Put Premium
#'
#' Calculates the premium of a European-style put option using the Black-Scholes option pricing model
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the value of the put option
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples putpremium(100, 100, 0.20, (45/365), 0.02, 0.02)
putpremium <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  x*exp(-r*t)*pnorm(-d2) - s*exp(-d*t) * pnorm(-d1)
}

#' Call Delta
#'
#' Calculates the delta of the European- style call option
#'
#' The delta of an option can be defined as the rate of change of the option value given a $1 change in the underlying asset price.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the call delta
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples calldelta(100, 100, 0.20, (45/365), 0.02, 0.02)
calldelta <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  exp(-d*t)*pnorm(d1)
}


#' Put Delta
#'
#' Calculates the delta of the European- style put option
#'
#' The delta of an option can be defined as the rate of change of the option value given a $1 change in the underlying asset price.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the put delta
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples putdelta(100, 0.20, (45/365), 0.02, 0.02)
putdelta <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  exp(-d*t)*(pnorm(d1)-1)
}


#' Call Theta
#'
#' Calculates the theta of the European- style call option
#'
#' Theta is the "time-decay" of the option value measured as a daily value
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the call theta
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples calltheta(100, 100, 0.20, (45/365), 0.02, 0.02)
calltheta <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - (sigma*sqrt(t))
  (1/365)*(1/t)*(-1*(((s*sigma*exp(-d*t))/2*sqrt(t))* (1/sqrt(2*pi)) * exp((-1*d1^2)/2))
                 - (r*x*exp(-r*t)*pnorm(d2)) + (d*s*exp(-d*t)*pnorm(d1)))
}

#' Put Theta
#'
#' Calculates the theta of the European- style put option
#'
#' Theta is the "time-decay" of the option value measured as a daily value.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the put theta
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples puttheta(100, 100, 0.20, (45/365), 0.02, 0.02)
puttheta <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - (sigma*sqrt(t))
  (1/365)*(1/t)*(-1*(((s*sigma*exp(-d*t))/2*sqrt(t))* (1/sqrt(2*pi)) * exp((-1*d1^2)/2))
                 + (r*x*exp(-r*t)*pnorm(-d2)) - (d*s*exp(-d*t)*pnorm(-d1)))
}

#' Call Rho
#'
#' Calculates the rho of the European- style call option
#'
#' Rho measures the change in the option's value given a 1% change in the interest rate.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the call rho
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples callrho(100, 100, 0.20, (45/365), 0.02, 0.02)
callrho <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  (x*(t)*exp(-r*(t))*pnorm(d2))/100
}

#' Put Rho
#'
#' Calculates the rho of the European- style put option
#'
#' Rho measures the change in the option's value given a 1% change in the interest rate.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the put rho
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples putrho(100, 100, 0.20, (45/365), 0.02, 0.02)
putrho <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  (-1/100)*(x*(t)*exp(-r*(t))*pnorm(-d2))
}

#' Option Gamma
#'
#' Calculates the gamma of a European- style call and put option
#'
#' Gamma is the rate of change of the option's delta given a $1 change in the underlying asset.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the option gamma
#' @export
#'
#' @examples optiongamma(100, 100, 0.20, (45/365), 0.02, 0.02)
optiongamma <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  (exp(-d*t)/(s*sigma*sqrt(t)))*(1/sqrt(2*pi))*(exp((-d1^2)/2))
}

#' Option Vega
#'
#' Calculates the vega of a European- style call and put option
#'
#' Vega measures the change in the option's value given a 1% change in the implied volatility.
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use r.cont
#' @param d Annual continuously-compounded dividend yield, use r.cont
#'
#' @return Returns the option vega
#' @export
#'
#' @examples optionvega(100, 100, 0.20, (45/365), 0.02, 0.02)
optionvega <- function(s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  d1 <- ((log(s/x)+((r-d+(sigma^2/2))*t)))/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  (1/100)*(s*exp(-d*t)*sqrt(t)*(1/sqrt(2*pi))*exp((-d1^2)/2))
}

#' Call Option Evaluation
#'
#'Creates a data.frame containing call option greeks; delta, gamma, vega, theta, rho and the call premium
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns a data.frame containing the option premium and greeks:
#' \itemize{
#' \item Premium
#' \item Delta
#' \item Gamma
#' \item Vega
#' \item Theta
#' \item Rho
#' }
#' @export
#' @author John T. Buynak
#'
#' @examples calleval(100, 100, 0.20, (45/365), 0.02, 0.02)
calleval <- function(s, x, sigma, t, r, d = 0){

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  data.frame(Premium = callpremium(s, x, sigma, t, r, d),
             Delta = calldelta(s, x, sigma, t, r, d),
             Gamma = optiongamma(s, x, sigma, t, r, d),
             Vega = optionvega(s, x, sigma, t, r, d),
             Theta = calltheta(s, x, sigma, t, r, d),
             Rho = callrho(s, x, sigma, t, r, d))

}

#' Put Option Evaluation
#'
#'Creates a data.frame containing put option greeks; delta, gamma, vega, theta, rho and the putpremium
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns a data.frame containing the option premium and greeks:
#' \itemize{
#' \item Premium
#' \item Delta
#' \item Gamma
#' \item Vega
#' \item Theta
#' \item Rho
#' }
#' @export
#' @author John T. Buynak
#'
#' @examples puteval(100, 100, 0.20, (45/365), 0.02, 0.02)
puteval <- function(s, x, sigma, t, r, d = 0){

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  data.frame(Premium = putpremium(s, x, sigma, t, r, d),
             Delta = putdelta(s, x, sigma, t, r, d),
             Gamma = optiongamma(s, x, sigma, t, r, d),
             Vega = optionvega(s, x, sigma, t, r, d),
             Theta = puttheta(s, x, sigma, t, r, d),
             Rho = putrho(s, x, sigma, t, r, d))

}

#' Dual Option Evaluation
#'
#'Creates a data.frame containing both call and put option greeks; delta, gamma, vega, theta, rho and the option premium
#'
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns a data.frame containing the call and put option premium and greeks:
#' \itemize{
#' \item Premium
#' \item Delta
#' \item Gamma
#' \item Vega
#' \item Theta
#' \item Rho
#' }
#' @export
#'
#' @examples opteval(100, 100, 0.20, (45/365), 0.02, 0.02)
opteval <- function(s, x, sigma, t, r, d = 0){
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm


  data.frame(Option = c("Call", "Put"), Premium = c(callpremium(s, x, sigma, t, r, d),
                                                    putpremium(s, x, sigma, t, r, d)), Delta = c(calldelta(s, x, sigma, t, r, d),
                                                                                                 putdelta(s, x, sigma, t, r, d)),
             Gamma = c(optiongamma(s, x, sigma, t, r, d), optiongamma(s, x, sigma, t, r, d)),
             Vega = c(optionvega(s, x, sigma, t, r, d), optionvega(s, x, sigma, t, r, d)),
             Theta = c(calltheta(s, x, sigma, t, r, d), puttheta(s, x, sigma, t, r, d)),
             Rho = c(callrho(s, x, sigma, t, r, d ), putrho(s, x, sigma, t, r, d)))

}

#' Call Option Greek
#'
#' Computes the selected option greek, including premium
#'
#' @param greek String value, desired option greek to return
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the desired option greek, including premium
#' @export
#'
#' @examples callgreek("delta", 100, 100, 0.20, (45/365), 0.02, 0.02)
#' callgreek("gamma", 100, 100, 0.20, (45/365), 0.02, 0.02)
callgreek <- function(greek = c("delta", "gamma", "theta", "vega", "rho", "premium"),
                      s, x, sigma, t, r, d = 0) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  if(!is.character(greek)){
    stop("The greek argument must be a character string!", call. = FALSE)
  }else{

    if(greek == "delta"){
      calldelta(s, x, sigma, t, r, d)
    }else if(greek == "gamma") {
      optiongamma(s, x, sigma, t, r, d)
    } else if(greek == "vega"){
      optionvega(s, x, sigma, t, r, d)
    } else if(greek == "theta") {
      calltheta(s, x, sigma, t, r, d)
    }else if(greek == "rho"){
      callrho(s, x, sigma, t, r, d)
    }else if(greek=="premium"){
      callpremium(s, x, sigma, t, r, d)
    }
  }
}

#' Put Option Greek
#'
#' Computes the selected option greek, including premium
#'
#' @param greek String value, desired option greek to return
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Returns the dired option greek, including premium
#' @export
#'
#' @examples putgreek("vega", 100, 100, 0.20, (45/365), 0.02, 0.02)
putgreek <- function(greek = c("delta", "gamma", "theta", "vega", "rho", "premium"),
                     s, x, sigma, t, r, d = 0) {

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  if(!is.character(greek)){
    stop("The greek argument must be a character string!", call. = FALSE)
  }else{

    if(greek == "delta"){
      putdelta(s, x, sigma, t, r, d)
    }else if(greek == "gamma") {
      optiongamma(s, x, sigma, t, r, d)
    } else if(greek == "vega"){
      optionvega(s, x, sigma, t, r, d)
    } else if(greek == "theta") {
      puttheta(s, x, sigma, t, r, d)
    }else if(greek == "rho"){
      putrho(s, x, sigma, t, r, d)
    }else if(greek=="premium"){
      putpremium(s, x, sigma, t, r, d)
    }
  }
}

#' Probability Between
#'
#'Calculates the probability of the underlying asset value falling between two prices in a designated time frame, given the daily standard devaiation of the underlying returns.
#'
#' This function has two separate possible operations:
#' 1. Calculates the probability of the underlying asset value falling between two prices in a designated time frame, given the daily standard devaiation of the underlying returns.
#' 2. Calculates the probable price range, given a set probability
#'
#' @param spot Current price of the underlying asset
#' @param lower Lower price of the range
#' @param upper Upper price of the range
#' @param mean The average daily price movement, default = 0
#' @param asd Annualized standard deviation of the underlying returns
#' @param dsd Daily standard deviation of the underlying returns (Annual vol/sqrt(256)), used as an alternative to the asd parameter in conjuction with the dte parameter
#' @param dte Days until expiration, designated time frame
#' @param p Designated probability
#' @param quantile Logical. If True, calculates the probable price range
#' @param tradedays Number of trade days in a year, default = 262
#'
#' @return Returns a probability (if quantile = FALSE), Returns a data.frame (if quantile = TRUE)
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#'
#' @examples prob.btwn(spot = 100, lower = 90, upper = 110, mean = 0, dsd = 0.01, dte = 45)
#' @examples prob.btwn(spot = 100, mean = 0, dsd = 0.01, dte = 45, p = 0.75, quantile = TRUE)
prob.btwn <- function(spot, lower, upper, asd = 0, dsd = 0, dte = 0, mean = 0, p, quantile = FALSE, tradedays = 262) {


  if(dsd == 0 & asd > 0 ){
    dsd <- asd/sqrt(tradedays)
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd > 0 & dte > 0){
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd == 0 & dte == 0){

    stop("Please evaluate inputs")

  }


  if(quantile == TRUE){

    #xsd <- dsd * sqrt(dte)
    p2 <- (1 - p)/2
    #1 - qnorm(p, mean, asd)
    tprob <- c(qnorm(p2, mean, asd), -1 * qnorm(p2, mean, asd))
    data.frame("probability" = p, "percent change" = tprob, "price" = spot * (1+ tprob))



  }else{

    #xsd <- dsd * sqrt(dte)

    lower.pc <- (lower - spot)/spot
    upper.pc <- (upper - spot)/spot

    plower <- pnorm(lower.pc, mean, asd)
    pupper <- pnorm(upper.pc, mean, asd)

    pupper - plower
  }
}


#' Probability Below
#'
#'Calculates the probability of the underlying asset value remaining below a price level in a designated time frame, given the daily standard devaiation of the underlying returns.
#'
#' This function has two separate possible operations:
#' 1. Calculates the probability of the underlying asset value remaining below a price level in a designated time frame, given the daily standard devaiation of the underlying returns.
#' 2. Calculates the price the asset will remain below, given the designated probability
#'
#' @param spot Current price of the underlying asset
#' @param upper Upper price of the range
#' @param mean The average daily price movement, default = 0
#' @param asd Annualized standard deviation of the underlying returns
#' @param dsd Daily standard deviation of the underlying returns (Annual vol/sqrt(256)), used as an alternative to the asd parameter in conjuction with the dte parameter
#' @param dte Days until expiration, designated time frame
#' @param p Designated probability
#' @param quantile Logical. If True, calculates the price the asset will remain below, given the designated probability
#' @param tradedays Number of trade days in a year, default = 262
#'
#'
#' @return Returns a probability (if quantile = FALSE), Returns a data.frame (if quantile = TRUE)
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#'
#' @examples prob.below(spot = 100, upper = 110, mean = 0, dsd = 0.01, dte = 45)
#' @examples prob.below(spot = 100, mean = 0, dsd = 0.01, dte = 45, p = 0.75, quantile = TRUE)
prob.below <- function(spot, upper, mean = 0, asd = 0, dsd = 0, dte = 0, p, quantile = FALSE, tradedays = 262) {

  if(dsd == 0 & asd > 0 ){
    dsd <- asd/sqrt(tradedays)
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd > 0 & dte > 0){
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd == 0 & dte == 0){

    stop("Please evaluate inputs")

  }

  if(quantile == TRUE){




    tprob <- c(qnorm(p, mean, asd))
    data.frame("probability" = p, "percent change" = tprob, "price" = spot * (1+ tprob))

  }else{


    upper.pc <- (upper - spot)/spot

    pupper <- pnorm(upper.pc, mean, asd)

    as.numeric(pupper)
  }
}

#' Probability Above
#'
#'Calculates the probability of the underlying asset value remaining above a price level in a designated time frame, given the daily standard devaiation of the underlying returns.
#'
#' This function has two separate possible operations:
#' 1. Calculates the probability of the underlying asset value remaining above a price level in a designated time frame, given the daily standard devaiation of the underlying returns.
#' 2. Calculates the price the asset will remain above, given the designated probability
#'
#' @param spot Current price of the underlying asset
#' @param lower Lower price of the range
#' @param mean The average daily price movement, default = 0
#' @param asd Annualized standard deviation of the underlying returns
#' @param dsd Daily standard deviation of the underlying returns (Annual vol/sqrt(256)), used as an alternative to the asd parameter in conjuction with the dte parameter
#' @param dte Days until expiration, designated time frame
#' @param p Designated probability
#' @param quantile Logical. If True, calculates the price the asset will remain above, given the designated probability
#' @param tradedays Number of trade days in a year, default = 262
#'
#' @return Returns a probability (if quantile = FALSE), Returns a data.frame (if quantile = TRUE)
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#'
#' @examples prob.above(spot = 100, lower = 110, mean = 0, dsd = 0.01, dte = 45)
#' @examples prob.above(spot = 100, mean = 0, dsd = 0.01, dte = 45, p = 0.75, quantile = TRUE)
prob.above <- function(spot, lower, mean = 0, asd = 0, dsd = 0, dte = 0, p, quantile = FALSE, tradedays = 262) {

  if(dsd == 0 & asd > 0 ){
    dsd <- asd/sqrt(tradedays)
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd > 0 & dte > 0){
    asd <- dsd * sqrt(dte)
  }else if(asd == 0 & dsd == 0 & dte == 0){

    stop("Please evaluate inputs")

  }


  if(quantile == TRUE){

    #xsd <- dsd * sqrt(dte)


    tprob <- c(-1 * qnorm(p, mean, asd))
    data.frame("probability" = p, "percent change" = tprob, "price" = spot * (1+ tprob))
  }else {

    #xsd <- dsd * sqrt(dte)

    lower.pc <- (lower - spot)/spot

    plower <- pnorm(lower.pc, mean, asd)

    1 - as.numeric(plower)
  }
}


#' Time Difference
#'
#'Computes the difference in time between two dates
#'
#' @param date1 Earlier date
#' @param date2 Later date
#' @param period String value, either "days", or "years"
#'
#' @return Returns a numeric value
#' @export
#'
#' @examples tdiff("2018-01-01", "2018-06-30", "days")
tdiff <- function(date1, date2, period = c("days, years")) {
  if(period == "days"){
    as.numeric(difftime(date2, date1, units = "days"))
  } else{
    as.numeric(difftime(date2, date1, units = "days"))/365
  }
}


#' Continuously Compounded Rate
#'
#' Convert a given nominal rate to a continuously compounded rate
#'
#' @param r nominal rate
#' @param n number of times compounded each year
#'
#' @return Returns a continuously compounded rate
#' @export
#'
#' @examples r.cont(0.12, 2)
r.cont <- function(r, n) {

  log((1+r/n)^n)
}


#' Plot Bull Call Spread
#'
#' Plot a bull call spread (debit spread)
#'
#' @param s Spot price of the underlying asset
#' @param x1 Lower-strike option price (long option)
#' @param x2 Higher-strike option price (short option)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the lower-strike option
#' @param sigma2 Annualized implied volatility of the higher-strike option
#' @param d Annual continuously compounded risk-free rate
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a vertical call spread (debit spread).
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotbullcall(s= 100, x1 = 95, x2 = 105, t = (45/365), r = 0.02,
#' sigma = 0.20, sigma2 = 0.20, d = 0, ll = 0.75, ul = 1.25)
plotbullcall <- function(s, x1, x2, t, r, sigma, sigma2 = sigma, d = 0, ll = 0.75, ul=1.25,
                         xlab = "spot", ylab = "profit/loss", main = "Bull Call Spread", ...) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  bullcall <- data.frame(spot = c((s*ll):(s*ul)))
  bullcall$bce <- ((callpremium(bullcall$spot, x1, sigma, t = 1e-5, r) -
                      callpremium(bullcall$spot, x2, sigma2, t = 1e-5, r)) -
                     callpremium(s, x1, sigma, t, r) +
                     callpremium(s, x2, sigma2, t, r))
  bullcall$bct0.5 <- ((callpremium(bullcall$spot, x1, sigma, t*0.5, r) -
                         callpremium(bullcall$spot, x2, sigma2, t*0.5, r)) -
                        callpremium(s, x1, sigma, t, r) +
                        callpremium(s, x2, sigma2, t, r))
  bullcall$bct <- ((callpremium(bullcall$spot, x1, sigma, t, r) -
                      callpremium(bullcall$spot, x2, sigma2, t, r)) -
                     callpremium(s, x1, sigma, t, r) +
                     callpremium(s, x2, sigma2, t, r))

  par(mfrow = c(1,1))
  plot(bullcall$spot, bullcall$bce, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
  abline("h" = 0)
  abline("v" = s)
  lines(bullcall$spot, bullcall$bct0.5, col = "red")
  lines(bullcall$spot, bullcall$bct, col = " blue")


}

#' Plot Bear Call Spread
#'
#' Plot a bear call spread (credit spread)
#'
#' @param s Spot price of the underlying asset
#' @param x1 Lower-strike option price (short option)
#' @param x2 Higher-strike option price (long option)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the lower-strike option
#' @param sigma2 Annualized implied volatility of the higher-strike option
#' @param d Annual continuously compounded risk-free rate
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a vertical call spread (credit spread).
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotbearcall(s= 100, x1 = 95, x2 = 105, t = (45/365), r = 0.02,
#' sigma = 0.20, sigma2 = 0.20, d = 0, ll = 0.75, ul = 1.25)
plotbearcall <- function(s, x1, x2, t, r, sigma, sigma2 = sigma,  d = 0, ll = 0.75, ul=1.25,
                         xlab = "spot", ylab = "Profit/Loss", main = "Bear Call Spread", ...) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  bearcall <- data.frame(spot = c((s*ll):(s*ul)))
  bearcall$bce <- ((callpremium(bearcall$spot, x2, sigma2, t = 1e-5, r) -
                      callpremium(bearcall$spot, x1, sigma, t = 1e-5, r)) -
                     callpremium(s, x2, sigma2, t, r) +
                     callpremium(s, x1, sigma, t, r))
  bearcall$bct0.5 <- ((callpremium(bearcall$spot, x2, sigma2, t*0.5, r) -
                         callpremium(bearcall$spot, x1, sigma, t*0.5, r)) -
                        callpremium(s, x2, sigma2, t, r) +
                        callpremium(s, x1, sigma, t, r))
  bearcall$bct <- ((callpremium(bearcall$spot, x2, sigma2, t, r) -
                      callpremium(bearcall$spot, x1, sigma, t, r)) -
                     callpremium(s, x2, sigma2, t, r) +
                     callpremium(s, x1, sigma, t, r))

  par(mfrow = c(1,1))
  plot(bearcall$spot, bearcall$bce, type = "l", xlab = xlab, ylab = ylab, main = main, ...) # first plot
  abline("h" = 0)
  abline("v" = s)
  lines(bearcall$spot, bearcall$bct0.5, col = "red")
  lines(bearcall$spot, bearcall$bct, col = " blue")


}

#' Plot Bull Put Spread
#'
#' Plot a bull put spread (credit spread)
#'
#' @param s Spot price of the underlying asset
#' @param x1 Lower-strike option price (long option)
#' @param x2 Higher-strike option price (short option)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the lower-strike option
#' @param sigma2 Annualized implied volatility of the higher-strike option
#' @param d Annual continuously compounded risk-free rate
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a vertical put spread (credit spread).
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotbullput(s= 100, x1 = 95, x2 = 105, t = (45/365), r = 0.02,
#' sigma = 0.20, sigma2 = 0.20, d = 0, ll = 0.75, ul = 1.25)
plotbullput <- function(s, x1, x2, t, r, d = 0, sigma, sigma2 = sigma, ll = 0.75, ul = 1.25,
                        xlab = "spot", ylab = "Profit/Loss", main = "Bull Put Spread", ...) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  bullput <- data.frame(spot = c((s*ll):(s*ul)))
  bullput$bce <- ((putpremium(bullput$spot, x1, sigma, t = 1e-5, r) -
                     putpremium(bullput$spot, x2, sigma2, t = 1e-5, r)) -
                    putpremium(s, x1, sigma, t, r) +
                    putpremium(s, x2, sigma2, t, r))
  bullput$bct0.5 <- ((putpremium(bullput$spot, x1, sigma, t*0.5, r) -
                        putpremium(bullput$spot, x2, sigma2, t*0.5, r)) -
                       putpremium(s, x1, sigma, t, r) +
                       putpremium(s, x2, sigma2, t, r))
  bullput$bct <- ((putpremium(bullput$spot, x1, sigma, t, r) -
                     putpremium(bullput$spot, x2, sigma2, t, r)) -
                    putpremium(s, x1, sigma, t, r) +
                    putpremium(s, x2, sigma2, t, r))

  par(mfrow = c(1,1))
  plot(bullput$spot, bullput$bce, type = "l", xlab = xlab, ylab = ylab, main = main, ...) # first plot
  abline("h" = 0)
  abline("v" = s)
  lines(bullput$spot, bullput$bct0.5, col = "red")
  lines(bullput$spot, bullput$bct, col = " blue")


}

#' Plot Bear Put Spread
#'
#' Plot a bear put spread (debit spread)
#'
#' @param s Spot price of the underlying asset
#' @param x1 Lower-strike option price (short option)
#' @param x2 Higher-strike option price (long option)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the lower-strike option
#' @param sigma2 Annualized implied volatility of the higher-strike option
#' @param d Annual continuously compounded risk-free rate
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a vertical put spread (debit spread).
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotbearput(s= 100, x1 = 95, x2 = 105, t = (45/365), r = 0.02,
#' sigma = 0.20, sigma2 = 0.20, d = 0, ll = 0.75, ul = 1.25)
plotbearput <- function(s, x1, x2, t, r, sigma, sigma2 = sigma, d = 0, ll = 0.75, ul = 1.25,
                        xlab = "spot", ylab = "Profit/Loss", main = "Bear Put Spread", ...) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  bearput <- data.frame(spot = c((s*ll):(s*ul)))
  bearput$bce <- ((putpremium(bearput$spot, x2, sigma2, t = 1e-5, r) -
                     putpremium(bearput$spot, x1, sigma, t = 1e-5, r)) -
                    putpremium(s, x2, sigma2, t, r) +
                    putpremium(s, x1, sigma, t, r))
  bearput$bct0.5 <- ((putpremium(bearput$spot, x2, sigma2, t*0.5, r) -
                        putpremium(bearput$spot, x1, sigma, t*0.5, r)) -
                       putpremium(s, x2, sigma2, t, r) +
                       putpremium(s, x1, sigma, t, r))
  bearput$bct <- ((putpremium(bearput$spot, x2, sigma2, t, r) -
                     putpremium(bearput$spot, x1, sigma, t, r)) -
                    putpremium(s, x2, sigma2, t, r) +
                    putpremium(s, x1, sigma, t, r))

  par(mfrow = c(1,1))
  plot(bearput$spot, bearput$bce, type = "l", xlab = xlab, ylab = ylab, main = main, ...) # first plot
  abline("h" = 0)
  abline("v" = s)
  lines(bearput$spot, bearput$bct0.5, col = "red")
  lines(bearput$spot, bearput$bct, col = " blue")


}

#' Plot Double Vertical Spread
#'
#' Plot a double vertical spread (credit spread)
#'
#' The double vertical spread consists of a credit put spread and a credit debit spread.
#'
#' @param s Spot price of the underlying asset
#' @param x1 Lower-strike put option price (long option)
#' @param x2 Higher-strike put option price (short option)
#' @param x3 Lower-strike call option price (short option)
#' @param x4 Higher-strike call option price (long option)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the lower-strike put option
#' @param sigma2 Annualized implied volatility of the higher-strike put option
#' @param sigma3 Annualized implied volatility of the lower-strike call option
#' @param sigma4 Annualized implied volatility of the higher-strike call option
#' @param d Annual continuously compounded risk-free rate
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a double vertical spread (credit spread).
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotdv(s= 100, x1 = 90, x2 = 95, x3 = 105, x4 = 110, t = (45/365), r = 0.02, sigma = 0.20)
plotdv <- function(s, x1, x2, x3, x4, t, r, sigma, sigma2 = sigma, sigma3 = sigma, sigma4 =sigma,
                   d = 0, ll = 0.75, ul = 1.25,
                   xlab = "spot", ylab = "Profit/Loss", main = "Double Vertical Spread", ...) {

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  dv <- data.frame(spot = c((s*ll):(s*ul)))
  dv$sume <- ((putpremium(dv$spot, x1, sigma, t = 1e-5, r) -
                 putpremium(dv$spot, x2, sigma2, t = 1e-5, r)) -
                putpremium(s, x1, sigma, t, r) +
                putpremium(s, x2, sigma2, t, r)) +
    ((callpremium(dv$spot, x4, sigma4, t = 1e-5, r) -
        callpremium(dv$spot, x3, sigma3, t = 1e-5, r)) -
       callpremium(s, x4, sigma4, t, r) +
       callpremium(s, x3, sigma3, t, r))

  dv$sumt0.5 <- ((putpremium(dv$spot, x1, sigma, t*0.5, r) -
                    putpremium(dv$spot, x2, sigma2, t*0.5, r)) -
                   putpremium(s, x1, sigma, t, r) +
                   putpremium(s, x2, sigma2, t, r)) +
    ((callpremium(dv$spot, x4, sigma4, t*0.5, r) -
        callpremium(dv$spot, x3, sigma3, t*0.5, r)) -
       callpremium(s, x4, sigma4, t, r) +
       callpremium(s, x3, sigma3, t, r))

  dv$sumt <- ((putpremium(dv$spot, x1, sigma, t, r) -
                 putpremium(dv$spot, x2, sigma2, t, r)) -
                putpremium(s, x1, sigma, t, r) +
                putpremium(s, x2, sigma2, t, r)) +
    ((callpremium(dv$spot, x4, sigma4, t, r) -
        callpremium(dv$spot, x3, sigma3, t, r)) -
       callpremium(s, x4, sigma4, t, r) +
       callpremium(s, x3, sigma3, t, r))


  par(mfrow = c(1,1))
  plot(dv$spot, dv$sume, type = "l", xlab = xlab, ylab = ylab, main = main, ...) # first plot
  abline("h" = 0)
  abline("v" = s)
  lines(dv$spot, dv$sumt0.5, col = "red")
  lines(dv$spot, dv$sumt, col = " blue")


}


#' Plot Custom Vertical Spread
#'
#' @param options String argument, either "call" or "put"
#' @param s Spot price of the underlying asset
#' @param x1 Short strike (either higher or lower)
#' @param x2 Long strike (either higher or lower)
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Annualized implied volatility of the short option
#' @param sigma2 Annualized implied volatility of the long option
#' @param d Annual continuously compounded dividend yield
#' @param ll Lower-limit of the plot, set as (desired price/spot)
#' @param ul Upper-limit of the plot, set as (desired price/spot)
#' @param xlab X-Axis Label
#' @param ylab Y-Axis Label
#' @param main Title of the plot
#' @param ... Additional plot parameters
#'
#' @return Returns a plot of a custom vertical spread.
#' Black line: The profit(loss) at expiration.
#' Red line: The profit(loss) at (1/2) time "t" ~ half-way to expiration.
#' Blue line: The profit(loss) at inception.
#' @author John T. Buynak
#' @export
#'
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#'
#' @examples plotvertical("call", 100, 90, 110, (45/365), 0.02, 0.20)
plotvertical <- function(options = c("call", "put"), s, x1, x2,
                         t, r, sigma, sigma2 = sigma, d = 0, ll = 0.75, ul = 1.25,
                         xlab = "spot", ylab = "profit/loss", main = "Vertical Spread", ...) {
  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  if(!is.character(options)){
    stop("The options argument must be a character string!", call. = FALSE)
  }else{

    if(options == "put") {
      vert <- data.frame(spot = c((s*ll):(s*ul)))

      vert$exp <- ((putpremium(vert$spot, x2, sigma2, t = 1e-5, r) -
                      putpremium(vert$spot, x1, sigma, t = 1e-5, r)) -
                     putpremium(s, x2, sigma2, t, r) +
                     putpremium(s, x1, sigma, t, r))

      vert$half <- ((putpremium(vert$spot, x2, sigma2, t*0.5, r) -
                       putpremium(vert$spot, x1, sigma, t*0.5, r)) -
                      putpremium(s, x2, sigma2, t, r) +
                      putpremium(s, x1, sigma, t, r))

      vert$begin <- ((putpremium(vert$spot, x2, sigma2, t, r) -
                        putpremium(vert$spot, x1, sigma, t, r)) -
                       putpremium(s, x2, sigma2, t, r) +
                       putpremium(s, x1, sigma, t, r))


    }

    else {
      vert <- data.frame(spot = c((s*ll):(s*ul)))

      vert$exp <- ((callpremium(vert$spot, x2, sigma2, t = 1e-5, r) -
                      callpremium(vert$spot, x1, sigma, t = 1e-5, r)) -
                     callpremium(s, x2, sigma2, t, r) +
                     callpremium(s, x1, sigma, t, r))

      vert$half <- ((callpremium(vert$spot, x2, sigma2, t*0.5, r) -
                       callpremium(vert$spot, x1, sigma, t*0.5, r)) -
                      callpremium(s, x2, sigma2, t, r) +
                      callpremium(s, x1, sigma, t, r))

      vert$begin <- ((callpremium(vert$spot, x2, sigma2, t, r) -
                        callpremium(vert$spot, x1, sigma, t, r)) -
                       callpremium(s, x2, sigma2, t, r) +
                       callpremium(s, x1, sigma, t, r))

    }
  }
  par(mfrow = c(1,1))
  plot(vert$spot, vert$exp, type = "l", xlab = xlab, ylab = ylab, main = main, ...) # first plot
  abline("h" = 0)
  abline("v" = s)
  lines(vert$spot, vert$half, col = "red")
  lines(vert$spot, vert$begin, col = " blue")
}

#' Vertical Spread Analytics
#'
#' Calculates the key analytics of a vertical spread
#'
#' @param options Character string. Either "call", or "put"
#' @param s Spot price of the underlying asset
#' @param x1 Strike price of the short option
#' @param x2 Strike price of the long option
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Implied volatility of the short option (annualized)
#' @param sigma2 Implied volatility of the long option (annualized)
#' @param vol Manual over-ride for the volatility of the underlying asset (annualized)
#' @param d Annual continuously compounded dividend yield
#'
#' @return Returns a data.frame
#' @export
#'
#' @examples vertical("call", s = 100, x1 = 90, x2 = 110, t = (45/365), r =  0.025, sigma = 0.20, vol = 0.25)
vertical <- function(options = c("call", "put"), s, x1, x2, t, r, sigma, sigma2 = sigma, vol = sigma, d = 0){

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  if(!is.character(options)){
    stop("The options argument must be a character string!", call. = FALSE)
  }

  if(options == "call"){
    c1 <- callpremium(s, x1, sigma, t, r, d)
    c2 <- callpremium(s, x2, sigma2, t, r, d)



    if(x1 > x2){
      be <- x2 + abs(c1 - c2)
      dsd <- (vol/sqrt(256))

      vert <- data.frame(Spot = s, Short.Strike = x1, Long.Strike = x2)
      vert$Max.Profit <- ((x1 - x2)+(c1 - c2))
      vert$Max.Loss <- (c1 - c2)
      vert$Breakeven <- x2 + abs(c1 - c2)
      vert$Prob.BE <- prob.above(s, lower = be, dsd = dsd, dte = t*365)
      vert$Prob.Max.Profit <- prob.above(s, lower = x1, dsd = dsd, dte = t*365)
      vert$Prob.Max.Loss <- prob.below(s, upper = x2, dsd = dsd, dte = t*365)
      vert$Initial.DC <- c1 - c2

      as.data.frame(t(round(vert,2)))
    }

    else{
      be <- x1 + (c1 - c2)
      dsd <- (vol/sqrt(256))

      vert <- data.frame(Spot = s, Short.Strike = x1, Long.Strike = x2)

      vert$Max.Profit <- (c1 - c2)
      vert$Max.Loss <- ((x1 - x2)+(c1 - c2))
      vert$Breakeven <- x1 + (c1 - c2)
      vert$Prob.BE <- prob.below(s, upper = be, dsd=dsd, dte= t*365)
      vert$Prob.Max.Profit <- prob.below(s, upper = x1, dsd=dsd, dte= t*365)
      vert$Prob.Max.Loss <- prob.above(s, lower = x2, dsd=dsd, dte= t*365)
      vert$Initial.DC <- c1 - c2

      as.data.frame(t(round(vert,2)))
    }

  }else{
    c1 <- putpremium(s, x1, sigma, t, r, d)
    c2 <- putpremium(s, x2, sigma2, t, r, d)


    if(x1 < x2){
      be <- x2 - abs(c1 - c2)
      dsd <- (vol/sqrt(256))

      vert <- data.frame(Spot = s, Short.Strike = x1, Long.Strike = x2)
      vert$Max.Profit <- ((x2 - x1)+(c1 - c2))
      vert$Max.Loss <- (c1 - c2)
      vert$Breakeven <- x2 - abs(c1 - c2)
      vert$Prob.BE <- prob.below(s, upper = be, dsd= dsd, dte= t*365)
      vert$Prob.Max.Profit <- prob.below(s, upper = x1, dsd= dsd, dte= t*365)
      vert$Prob.Max.Loss <- prob.above(s, lower = x2, dsd= dsd, dte= t*365)
      vert$Initial.DC <- c1 - c2

      as.data.frame(t(round(vert,2)))
    }

    else{
      be <- x1 - (c1 - c2)
      dsd <- (vol/sqrt(256))

      vert <- data.frame(Spot = s, Short.Strike = x1, Long.Strike = x2)

      vert$Max.Profit <- (c1 - c2)
      vert$Max.Loss <- ((x2 - x1)+(c1 - c2))
      vert$Breakeven <- x1 - (c1 - c2)
      vert$Prob.BE <- prob.above(s, lower = be, dsd = dsd, dte = t*365)
      vert$Prob.Max.Profit <- prob.above(s, lower = x1, dsd = dsd, dte = t*365)
      vert$Prob.Max.Loss <- prob.below(s, upper = x2, dsd = dsd, dte = t*365)
      vert$Initial.DC <- c1 - c2

      as.data.frame(t(round(vert,2)))
    }
  }

}

#' Double Vertical Spread Analytics
#'
#' Calculates the key analytics of a Double Vertical Credit Spread
#'
#' @param s Spot price of the underlying asset
#' @param x1 Strike price of the lower strike (long) put option
#' @param x2 Strike price of the higher strike (short) put option
#' @param x3 Strike price of the lower strike (short) call option
#' @param x4 Strike price of the higher strike (long) call option
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param sigma Implied volatility of the lower strike (long) put option (annualized)
#' @param sigma2 Implied volatility of the higher strike (short) put option (annualized)
#' @param sigma3 Implied volatility of the lower strike (short) call option (annualized)
#' @param sigma4 Implied volatility of the higher strike (long) call option (annualized)
#' @param vol Manual over-ride for the volatility of the underlying asset (annualized)
#' @param d Annual continuously compounded dividend yield
#'
#' @return Returns a data.frame
#' @export
#'
#' @examples dv(s = 100, x1 = 90, x2 = 95, x3 = 105, x4 = 110, t = 0.08, r = 0.02, sigma = 0.2, vol = 0.3)
dv <- function(s, x1, x2, x3, x4, t, r, sigma, sigma2 = sigma, sigma3 = sigma, sigma4 = sigma, vol = sigma, d = 0){

  ttm <- t
  if(t == 0){
    ttm <- 1e-200
  }else{
    ttm <- t
  }

  t <- ttm

  o1 <- putpremium(s, x1, sigma, t, r, d)
  o2 <- putpremium(s, x2, sigma2, t, r, d)
  o3 <- callpremium(s, x3, sigma3, t, r, d)
  o4 <- callpremium(s, x4, sigma4, t, r, d)

  be1 <- x2 - ((o2 + o3) - (o1 + o4))
  be2 <- x3 + ((o2 + o3) - (o1 + o4))

  dsd <- vol/sqrt(256)
  dte <- t * 365

  dvert <- data.frame(Spot = s, Long.Put = x1, Short.Put = x2, Short.Call = x3, Long.Call = x4)
  dvert$Lower.BE <- be1
  dvert$Higher.BE <- be2
  dvert$Prob.BE <- prob.btwn(s, be1, be2, dsd = dsd, dte = dte)
  dvert$Prob.Max.Profit <- prob.btwn(s, x2, x3, dsd = dsd, dte = dte)
  dvert$Prob.Max.Loss <- (prob.below(s, x1, dsd = dsd, dte = dte)) + (prob.above(s, x4, dsd = dsd, dte = dte))

  as.data.frame(t(round(dvert,2)))

}



#' Implied Volatility Calculation
#'
#' Computes the implied volatility of an option, either a call or put, given the option premium and key parameters
#'
#' @param type String argument, either "call" or "put"
#' @param price Current price of the option
#' @param s Spot price of the underlying asset
#' @param x Strike Price of the underlying asset
#' @param t Time to expiration in years
#' @param r Annual continuously compounded risk-free rate
#' @param d Annual continuously compounded dividend yield
#'
#' @return Returns a single option's implied volatility
#' @export
#'
#' @examples iv.calc(type = "call", price = 2.93, s = 100, x = 100, t = (45/365), r = 0.02, d = 0)
iv.calc <- function(type, price, s, x, t, r, d=0) {

  dvol <- function(type, price, volatility=0) {
    if(type == "call"){
      price -  callpremium(s = s, x = x, sigma = volatility, t = t, r = r, d = d)
    } else if(type == "put"){
      price -  putpremium(s = s, x = x, sigma = volatility, t = t, r = r, d = d)
    } else{
      stop("type must be either 'call' or 'put' ")
    }
  }



  volchart <- data.frame(vol = seq(0, 5, 0.001))
  volchart$distance <- dvol(type, price, volchart$vol)


  volchart$vol[match(min(abs(volchart$distance)), abs(volchart$distance))]

}


#' Lambda
#'
#' Calculates the Lambda of the call or put option
#'
#' Lambda, or elasticity is the percentage change in the option valueper percentage change in the underlying price. It is a measure of leverage.
#'
#' @param type Character string, either "call" or "put"
#' @param s Spot price of the underlying asset
#' @param x Strike price of the option
#' @param sigma Implied volatility of the underlying asset price, defined as the annualized standard deviation of the asset returns
#' @param t Time to maturity in years
#' @param r Annual continuously-compounded risk-free rate, use the function r.cont
#' @param d Annual continuously-compounded dividend yield, use the function r.cont
#'
#' @return Calculates the Lambda of the option contract
#' @export
#'
#' @examples lambda(type = "put", s = 100, x = 100, sigma = 0.15, t = 45/365, r = 0.02)
lambda <- function(type = "call", s, x, sigma, t, r, d = 0){

  if(type == "call"){
    calldelta(s = s, x = x, t = t, sigma = sigma, r =r, d = d) *
      ( s/callpremium(s = s, x = x, t = t, sigma = sigma, r =r, d = d))
  }else{

    putdelta(s = s, x = x, t = t, sigma = sigma, r =r, d = d) *
      ( s/putpremium(s = s, x = x, t = t, sigma = sigma, r =r, d = d))
  }

}
