InitClust <- function(data, G, k = 3, knnd.trim = 0.5, modelName = "VVV") {
  n <- nrow(data)
  knnd.trim <- min(c(max(c(0, floor(n - 2 * G)/n)), knnd.trim))
  distances <- dist(data)
  kthNND <- rep(0, n)
  Labels <- 1:(n - 1)
  kthNND[1] <- sort(distances[Labels])[k]
  Labels[(2):(n - 1)] <- Labels[(2):(n - 1)] + (n - 1 - 1)
  for (i in 2:n) {
    kthNND[i] <- sort(distances[Labels])[k]
    Labels[1:(i - 1)] <- Labels[1:(i - 1)] + 1
    Labels[(i + 1):(n - 1)] <- Labels[(i + 1):(n - 1)] + (n - i - 1)
  }
  kthNND <- kthNND[1:n]
  idx.reg <- {
    kthNND <= sort(kthNND)[n - floor(n * knnd.trim)]
  }
  redo <- FALSE
  ic <- NULL
  try({
    H <- hc(data = data[idx.reg, ], minclus = G, modelName = modelName)
    Hcluster <- as.vector(hclass(H, G))
    ic <- rep(0, n)
    ic[idx.reg] <- Hcluster
  }, silent = TRUE)
  if (is.null(ic)) {
    redo <- TRUE
  }
  else {
    iam <- .Cluster2Assign(ic)
    nj <- rep(0, times = G)
    for (j in 1:G) {
      nj[j] <- .CountUniqueRows(data[iam[, 1 + j] == 1, ])
    }
    if (any(nj < 2)) {
      redo <- TRUE
    }
    else {
      ans <- ic
      redo <- FALSE
    }
  }
  if (redo) {
    ic <- NULL
    try({
      dmat <- dist(data[idx.reg, ], method = "manhattan")
      H <- hclust(dmat, method = "average")
      Hcluster <- cutree(H, k = G)
      ic <- rep(0, n)
      ic[idx.reg] <- Hcluster
    }, silent = TRUE)
    if (is.null(ic)) {
      redo <- TRUE
    }
    else {
      iam <- .Cluster2Assign(ic)
      nj <- rep(0, times = G)
      for (j in 1:G) {
        nj[j] <- .CountUniqueRows(data[iam[, 1 + j] == 1, ])
      }
      if (any(nj < 2)) {
        redo <- TRUE
      }
      else {
        redo <- FALSE
        ans <- ic
      }
    }
  }
  if (redo) {
    ic <- NULL
    iter <- 1L
    while (iter <= 10 & redo) {
      kf <- cut(1:n, breaks = G, labels = FALSE)
      ric <- sample(kf, size = n, replace = FALSE)
      ic <- rep(0, n)
      ic[idx.reg] <- ric[idx.reg]
      iam <- .Cluster2Assign(ic)
      nj <- rep(0, times = G)
      for (j in 1:G) {
        nj[j] <- .CountUniqueRows(data[iam[, 1 + j] == 1, ])
      }
      if (any(nj < 2)) {
        iter <- iter + 1L
      }
      else {
        ans <- ic
        redo <- FALSE
      }
    }
  }
  if (redo) {
    StopMsg <- "\nInitClust() could not find a valid initial partition. A an initial partition may be provided based on a different strategy [see detail in help(InitClust)]\n\n"
    cat("\n")
    stop(StopMsg)
  }
  else {
    return(ans)
  }
}
