## confint method by Lucca Scrucca
confint.maxLik <- function(object, parm, level = 0.95, ...)
{
  cf <- coef(object)
  if(missing(parm)) 
    parm <- seq_along(cf)
  pnames <- names(cf)
  if(is.null(pnames)) 
    pnames <- parm
  else if(is.numeric(parm)) 
         parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  q <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), 
              dimnames = list(parm, pct))
  se <- sqrt(diag(vcov(object)))[parm]
  ci[] <- cf[parm] + se %o% q
  return(ci)
}

format.perc <- function(probs, digits) 
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
