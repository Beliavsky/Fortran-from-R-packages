#'
#'   braun.R
#'
#'  Braun (1980) method for composite null hypothesis
#'

braun <- function(U, simpletest, m) {
  n <- length(U)
  if(n < 2 * m) stop("Unsufficient data for Braun's method")
  #' split data into m groups
  group <- factor(sample(seq_len(n) %% m))
  #' apply the simple-null test to each subset
  zz <- by(data=U, INDICES=group, FUN=simpletest, simplify=FALSE)
  statistics <- sapply(zz, getElement, "statistic")
  pvalues    <- sapply(zz, getElement, "pvalue")
  statname   <- zz[[1]]$statname
  #' combine
  statistic <- max(statistics)
  pvalue <- 1 - (1 - min(pvalues))^m
  statname <- paste0(statname, "max")
  return(list(statistic=statistic, pvalue=pvalue, statname=statname))
}

