
require_tibble_package <- function () {
  if (! requireNamespace("tibble", quietly = TRUE)) {
    stop("The `tibble` package must be installed to use tidy() or glance() methods")
  }
}


tidy.maxLik <- function (x,  ...) {
  require_tibble_package()
  
  s <- summary(x)
  ret <- tibble::as_tibble(s$estimate, rownames = "term")
  colnames(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")
  
  ret
}


glance.maxLik <- function (x, ...) {
  require_tibble_package()

  ll <- logLik(x)
  nobs <- tryCatch(nObs(x), error = function(e) NA)
                           # nobs = NA in case of error
  ret <- tibble::tibble(
           df     = attr(ll, "df"),
           logLik = as.numeric(ll),
           AIC    = AIC(x),
           nobs   = nobs
         )
  
  ret
}
