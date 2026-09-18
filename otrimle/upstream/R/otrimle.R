otrimle <- function(data, G, initial = NULL, logicd = NULL, npr.max = 0.5, erc = 20, 
  beta = 0, iter.max = 500, tol = 1e-06, ncores = NULL, monitor = TRUE) {
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }
  if (!is.matrix(data)) {
    data <- data.matrix(data)
  }
  if (!is.numeric(data)) {
    stop("\"data\" must be numeric.")
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
  if (npr.max == 0) {
    logicd <- -Inf
  }
  else if (npr.max < 0 | npr.max >= 1) {
    stop("\"npr.max\" must belog to the interval [0,1). Tipically npr.max<=0.5.")
  }
  if (erc < 1) {
    stop("\"ecr\" must be larger or equal to 1.")
  }
  if (length(beta) > 1 | any(beta < 0)) {
    stop("\"beta\" must be a non-negative number.")
  }
  if (iter.max < 1) {
    stop("\"iter.max\" must be an integer larger or equal to 1.")
  }
  if (tol < 0) {
    stop("\"tol\" must be positive.")
  }
  if (is.null(ncores)) {
    DetectedCores <- detectCores()
    if (DetectedCores <= 2) {
      ncores <- 1
    }
    else {
      ncores <- {
        DetectedCores - 1
      }
    }
  }
  else {
    ncores <- as.integer(ncores)
    if (ncores <= 0) {
      stop("\"ncores\" must be an integer larger or equal to 1.")
    }
  }
  if (!is.null(logicd) & !is.numeric(logicd)) {
    stop("\"logicd\" must be  numeric vector, typically length(logicd) > 1.")
  }
  if (!is.null(logicd)) {
    logicd <- sort(logicd)
  }
  else {
    logicd <- c(-Inf, seq(-700, -100, by = 50), seq(-95, -55, by = 5), seq(-50, 
      -10, by = 2.5), seq(-9, 0, by = 1))
  }
  nGrid <- length(logicd)
  nn <- .CountUniqueRows(data)
  if (all(logicd == -Inf)) {
    if (nn <= G) {
      cat("\n")
      StopMsg <- paste0("There are m=", nn, " distinct data points in the data set. With logicd = -Inf, RIMLE existence requires m > G. Try to decrease 'G'.\n\n")
      stop(StopMsg)
    }
  }
  if (any(logicd > -Inf)) {
    if (nn <= {
      G + ceiling(nrow(data) * npr.max)
    }) {
      cat("\n")
      StopMsg <- paste0("There are m = ", nn, " distinct data points in the data set. With  logicd > -Inf, RIMLE existence requires m > G + ceiling(nrow(data)*npr.max). Try to vary 'npr.max', 'logicd', and/or 'G'.\n\n")
      stop(StopMsg)
    }
  }
  if (monitor) {
    message("\nOTRIMLE:")
  }
  if (is.null(initial)) {
    if (monitor) {
      message("........computing the initial partition")
    }
    initial <- .Cluster2Assign(InitClust(data = data, G = G, knnd.trim = npr.max))
  }
  else {
    if (monitor) {
      message("........checking the initial partition")
    }
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
  A <- data.frame(logicd = logicd, criterion = NA, iloglik = NA, enpr = NA, code = NA, 
    flag = NA, row.names = NULL, stringsAsFactors = FALSE)
  rownames(A) <- 1:nGrid
  if (monitor) {
    message("........searching for the optimal log(icd)")
  }
  if (ncores == 1) {
    for (k in 1:nGrid) {
      A[k, ] <- .ECM(data = data, initial = initial, logicd = A[k, 1], npr.max = npr.max, 
        erc = erc, beta = beta, iter.max = iter.max, tol = tol, grid.operation = TRUE)
    }
  }
  if (ncores > 1) {
    registerDoParallel(ncores)
    B <- foreach(k = 1:nGrid) %dopar% {
      .ECM(data = data, initial = initial, logicd = A[k, 1], npr.max = npr.max, 
        erc = erc, beta = beta, iter.max = iter.max, tol = tol, grid.operation = TRUE)
    }
    stopImplicitCluster()
    for (k in 1:nGrid) {
      A[k, ] <- B[[k]]
    }
  }
  if (all(is.na(A[, 2]))) {
    ans <- list()
    ans$code <- 0
    ans$optimization <- A
  }
  else {
    x <- A[which.min(A[, 2]), 1]
    if (monitor) {
      message("........computing the optimal solution")
    }
    ans <- .ECM(data = data, initial = initial, logicd = x, npr.max = npr.max, 
      erc = erc, beta = beta, iter.max = iter.max, tol = tol, grid.operation = FALSE)
    ans$optimization <- A
    if (monitor) {
      message("........done\n")
    }
  }
  class(ans) <- "otrimle"
  return(ans)
}
