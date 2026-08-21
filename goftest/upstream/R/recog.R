##  recog.R
##
## $Revision: 1.4 $ $Date: 2014/06/24 02:13:35 $
##

recogniseCdf <- function(s="punif") {
  if(!is.character(s) || length(s) != 1 || nchar(s) == 0) return(NULL)
  #' strip off the leading 'p' if present
  root <- if(substr(s,1,1) == "p") substr(s, 2, nchar(s)) else s
  a <- switch(root,
              beta     = "beta",
              binom    = "binomial",
              birthday = "birthday coincidence",
              cauchy   = "Cauchy",
              chisq    = "chi-squared",
              exp      = "exponential",
              f        = "F",
              gamma    = "Gamma",
              geom     = "geometric",
              hyper    = "hypergeometric",
              lnorm    = "log-normal",
              logis    = "logistic",
              nbinom   = "negative binomial",
              norm     = "Normal",
              pois     = "Poisson",
              t        = "Student's t",
              tukey    = "Tukey (Studentized range)",
              unif     = "uniform",
              weibull  = "Weibull",
              NULL)
  if(!is.null(a))
    return(paste(a, "distribution"))
  b <- switch(root,
              AD     = "Anderson-Darling",
              CvM    = "Cramer-von Mises",
              wilcox = "Wilcoxon Rank Sum",
              NULL)
  if(!is.null(b))
    return(paste("null distribution of", b, "Test Statistic"))
  return(NULL)
}

#' not exported

getfunky <- function(fname) {
  a <- mget(fname, mode="function", ifnotfound=list(NULL), inherits=TRUE)[[1]]
  return(a)
}

getCdf <- function(s="punif", fatal=TRUE) {
  sname <- deparse(substitute(s), nlines=1L)
  if(is.function(s)) return(s)
  if(is.character(s) && length(s) == 1 && nchar(s) > 0) {
    #' first try adding a leading 'p' (to catch the case s="t")
    if(substr(s,1,1) != "p") {
      f <- getfunky(paste0("p", s))
      if(is.function(f))
        return(f)
    }
    f <- getfunky(s)
    if(is.function(f))
      return(f)
  }
  if(fatal)
    stop(paste("Argument", sQuote(sname),
               "should be a function, or the name of a function"),
         call.=FALSE)
  return(NULL)
}
