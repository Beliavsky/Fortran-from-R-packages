.wecdf_aux <- function(x, weights) {
  ox <- order(x)
  x <- x[ox]
  ow <- weights[ox]
  n <- length(x)
  rval <- approxfun(x, cumsum(ow), method = "constant", yleft = 0, yright = 1, 
    f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
.wecdf <- function(x, weights) {
  sw <- sum(weights)
  if (sw == 0) {
    ans <- rep(0, length(x))
  }
  else {
    ans <- .wecdf_aux(x, weights/sw)(x)
  }
  return(ans)
}
.rwhich.max <- function(x) {
  xmax <- max(x, na.rm = TRUE)
  tmp <- which(x == xmax)
  if (length(tmp) == 1) {
    return(tmp)
  }
  else {
    return(sample(tmp, size = 1, replace = FALSE))
  }
}
.rwhich.min <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  tmp <- which(x == xmin)
  if (length(tmp) == 1) {
    return(tmp)
  }
  else {
    return(sample(tmp, size = 1, replace = FALSE))
  }
}
.CountUniqueRows <- function(x) {
  if (is.vector(x)) {
    if (length(x) == 0) {
      return(0)
    }
    else {
      return(length(unique(x)))
    }
  }
  else {
    n <- nrow(x)
    if (n == 0) {
      return(0)
    }
    else {
      if (ncol(x) == 1) {
        return(length(unique(x[, 1])))
      }
      else {
        sx <- x[do.call(order, as.list(as.data.frame(x))), ]
        flags <- rep(0L, times = n)
        for (i in 2:n) {
          flags[i] <- prod(sx[i, ] == sx[i - 1, ])
        }
        return({
          n - sum(flags)
        })
      }
    }
  }
}
.Cluster2Assign <- function(cluster) {
  cl <- sort(unique(cluster[cluster != 0]))
  G <- length(cl)
  N <- length(cluster)
  A <- matrix(0, nrow = N, ncol = 1 + G)
  if (sum(cluster == 0) > 0) {
    A[cluster == 0, 1] <- 1
  }
  for (j in 1:G) {
    A[cluster == cl[j], j + 1] <- 1
  }
  return(A)
}
