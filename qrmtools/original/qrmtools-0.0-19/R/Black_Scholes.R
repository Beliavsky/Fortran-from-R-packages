### Compute the Black--Scholes formula and the Greeks ##########################

##' @title Compute the Black--Scholes Formula
##' @param t initial/current time (in years)
##' @param S stock price at time t
##' @param r risk-free annual interest rate
##' @param sigma annual volatility (standard deviation)
##' @param K strike
##' @param T maturity (in years)
##' @param type character string indicating whether the price
##'        of a call or put option is to be computed
##' @return option price
##' @author Marius Hofert
Black_Scholes <- function(t, S, r, sigma, K, T, type = c("call", "put"))
{
    d1 <- (log(S/K) + (r + sigma^2/2) * (T-t)) / (sigma * sqrt(T-t))
    d2 <- d1 - sigma * sqrt(T-t)
    type <- match.arg(type)
    switch(type,
           "call" = {
               S * pnorm(d1) - K * exp(-r*(T-t)) * pnorm(d2)
           },
           "put" =  {
               -S * pnorm(-d1) + K * exp(-r*(T-t)) * pnorm(-d2)
           },
           stop("Wrong type"))
}

##' @title Compute the Greeks
##' @param t initial/current time (in years)
##' @param S stock price at time t
##' @param r risk-free annual interest rate
##' @param sigma annual volatility (standard deviation)
##' @param K strike
##' @param T maturity (in years)
##' @param type character string indicating whether the price
##'        of a call or put option is to be computed
##' @return Greeks
##' @author Marius Hofert
##' @note See https://en.wikipedia.org/wiki/Greeks_(finance)#Gamma for q = 0, tau = T-t
Black_Scholes_Greeks <- function(t, S, r, sigma, K, T, type = c("call", "put"))
{
    ## Basics
    d1 <- (log(S/K) + (r + sigma^2/2) * (T-t)) / (sigma * sqrt(T-t))
    d2 <- d1 - sigma * sqrt(T-t)

    ## 'type'-independent Greeks
    ## First-order
    vega <- S * dnorm(d1) * sqrt(T-t)
    ## Second-order (only the most important ones)
    gamma <- dnorm(d1) / (S * sigma * sqrt(T-t))
    vanna <- -dnorm(d1) * d2 / sigma
    vomma <- vega * d1 * d2 / sigma

    ## 'type'-dependent Greeks
    type <- match.arg(type)
    switch(type,
           "call" = {
               ## First-order derivatives
               delta <- pnorm(d1)
               theta <- -(S * dnorm(d1) * sigma) / (2 * sqrt(T-t)) -
                   r * K * exp(-r * (T-t)) * pnorm(d2)
               rho <- K * (T-t) * exp(-r * (T-t)) * pnorm(d2)
           },
           "put" =  {
               ## First-order derivatives
               delta <- -pnorm(-d1)
               theta <- -(S * dnorm(d1) * sigma) / (2 * sqrt(T-t)) +
                   r * K * exp(-r * (T-t)) * pnorm(-d2)
               rho <- -K * (T-t) * exp(-r * (T-t)) * pnorm(-d2)
           },
           stop("Wrong type"))

    ## Return
    res <- cbind(
        ## First-order
        delta = delta,
        theta = theta,
        rho = rho,
        vega = vega,
        ## Second-order
        gamma = gamma,
        vanna = vanna,
        vomma = vomma
    )
    if(nrow(res) == 1) drop(res) else res
}
