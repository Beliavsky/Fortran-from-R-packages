rimle <- function(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 20, 
  iter.max = 500, tol = 1e-06) {
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.matrix(data)) {
    data <- data.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("\n'data'must be numeric.\n\n")
  }
  if (!is.null(initial) & {
    !is.numeric(initial) | any(initial < 0) | {
      length(initial) != nrow(data)
    }
  }) {
    stop("\"initial\" must be an integer vector of length equal to nrow(data).")
  }
  G <- as.integer(G)
  if (G < 1) {
    stop("G must be an integer greater or equal to one.")
  }
  if (!is.null(initial)) {
    if (G != max(initial)) {
      stop("G is not equal to max(initial).")
    }
  }
  if (!is.null(logicd) & {
    !is.numeric(logicd) | length(logicd) != 1
  }) {
    stop("\"logicd\" must be numeric of length=1.")
  }
  if (npr.max == 0) {
    logicd <- -Inf
  }
  else if (npr.max < 0 | npr.max >= 1) {
    stop("\"npr.max\" must belog to the interval [0,1). Tipically npr.max<=0.5.")
  }
  if (erc < 1) {
    stop("\"ecr\" must be larger or equal to 1.")
  }
  if (iter.max < 1) {
    stop("\"iter.max\" must be an integer larger or equal to 1.")
  }
  if (tol < 0) {
    stop("\"tol\" must be positive.")
  }
  if (is.null(initial)) {
    initial <- .Cluster2Assign(InitClust(data = data, G = G, knnd.trim = npr.max))
  }
  else {
    initial <- .Cluster2Assign(initial)
    G <- ncol(initial) - 1
    nj <- rep(0, times = G)
    for (j in 1:G) {
      nj[j] <- .CountUniqueRows(data[initial[, 1 + j] == 1, ])
    }
    if (any(nj < 2)) {
      StopMsg <- paste("Some of the initial clusters contain less then 2 distinct points. Initial parameters cannot be computed. See details in help(InitClust) and try a different initial partition.\n\n", 
        sep = "")
      cat("\n")
      stop(StopMsg)
    }
  }
  if (is.null(logicd)) {
    N <- nrow(data)
    P <- ncol(data)
    G <- ncol(initial) - 1
    eps <- sqrt(.Machine$double.eps)
    gausscost <- P * log(2 * pi)
    qc <- qchisq({
      1 - 0.05
    }^{
      1/N
    }, df = P)
    logicds <- rep(0, G)
    for (j in 1:G) {
      idx <- initial[, 1 + j] == 1
      u <- eigen(cov(as.matrix(data[idx, ])), symmetric = TRUE, only.values = TRUE)$values
      if ({
        prod(abs(u)) < .Machine$double.xmin
      } | any(u <= 0) | any(!is.finite(u))) {
        Lmax <- max(u)
        Lmin <- min(u)
        u[u < Lmax/eps] <- Lmax/eps
      }
      logicds[j] <- -{
        gausscost + sum(log(u)) + qc
      }/2
      if (!is.finite(logicds[j])) {
        logicds[j] <- log(.Machine$double.xmin)
      }
    }
    logicd <- min(logicds)
  }
  nn <- .CountUniqueRows(data)
  if (logicd == -Inf) {
    if (nn <= G) {
      cat("\n")
      StopMsg <- paste0("There are m = ", nn, " distinct data points in the data set. With logicd = -Inf, RIMLE existence requires m > G. Try to decrease 'G'.\n\n")
      stop(StopMsg)
    }
  }
  if (logicd > -Inf) {
    if (nn <= {
      G + ceiling(nrow(data) * npr.max)
    }) {
      cat("\n")
      StopMsg <- paste0("There are m = ", nn, " distinct data points in the data set. With  logicd > -Inf, RIMLE existence requires m > G + ceiling(nrow(data)*npr.max). Try to vary 'npr.max', 'logicd', and/or 'G'. \n\n")
      stop(StopMsg)
    }
  }
  ans <- .ECM(data = data, initial = initial, logicd = logicd, npr.max = npr.max, 
    erc = erc, beta = 0, iter.max = iter.max, tol = tol, grid.operation = FALSE)
  class(ans) <- "rimle"
  return(ans)
}
