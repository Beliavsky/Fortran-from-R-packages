otrimleg <- function(dataset, G = 1:6, multicore = TRUE, ncores = detectCores(logical = FALSE) - 
  1, erc = 20, beta0 = 0, fixlogicd = NULL, monitor = 1, dmaxq = qnorm(0.9995)) {
  rung <- function(k) {
    if (monitor >= 1) 
      cat("Number of clusters: ", k, "\n")
    if (is.null(fixlogicd)) {
      solution <- otrimle(data = X, G = k, , ncores = 1, erc = erc, beta = beta0, 
        , monitor = monitor)
      cat("fixlogicd is null; code= ", solution$code, "iloglik=", solution$iloglik, 
        "\n")
    }
    else {
      solution <- rimle(data = X, G = k, logicd = fixlogicd[k])
      cat("logicd is fixed I; code= ", solution$code, "iloglik=", solution$iloglik, 
        "\n")
      if (solution$code == 0) {
        solution <- otrimle(data = X, G = k, ncores = 1, erc = erc, beta = beta0, 
          , monitor = monitor)
        cat("otrimle; code= ", solution$code, "iloglik=", solution$iloglik, 
          "\n")
      }
    }
    solution
  }
  X <- as.matrix(dataset)
  p <- ncol(X)
  n <- nrow(X)
  if (length(fixlogicd) == 1) 
    fixlogicd <- rep(fixlogicd, max(G))
  out <- init <- list()
  out$G <- G
  out$npar <- rep(NA, max(G))
  for (k in G) out$npar[k] <- k + k * p + k * (p + 1) * p/2
  out$ibic <- out$criterion <- out$iloglik <- out$noiseprob <- out$logicd <- out$denscrit <- rep(NA, 
    max(G))
  out$solution <- out$ddpm <- list()
  if (multicore) 
    slist <- mclapply(G, rung, mc.cores = ncores)
  else slist <- lapply(G, rung)
  i <- 1
  if (monitor >= 1) 
    cat("Computing criteria...\n")
  for (k in G) {
    out$solution[[k]] <- slist[[i]]
    i <- i + 1
    if (class(out$solution[[k]]) == "try-error" | out$solution[[k]]$code == 0) 
      out$criterion[k] <- out$noiseprob[k] <- out$denscrit[k] <- out$ddpm[[k]] <- out$iloglik[k] <- out$ibic[k] <- NA
    else {
      out$iloglik[k] <- out$solution[[k]]$iloglik
      out$ibic[k] <- -2 * out$iloglik[k] + out$npar[k] * log(n)
      out$criterion[k] <- out$solution[[k]]$criterion
      out$logicd[k] <- out$solution[[k]]$logicd
      out$noiseprob[k] <- out$solution[[k]]$exproportion[1]
      ddc <- kerndenscluster(X, fit = out$solution[[k]], maxq = dmaxq)
      out$denscrit[k] <- ddc$measure
      out$ddpm[[k]] <- ddc$ddpm
    }
  }
  out
}
otrimlesimg <- function(dataset, G = 1:6, multicore = TRUE, ncores = detectCores(logical = FALSE) - 
  1, erc = 20, beta0 = 0, simruns = 20, sim.est.logicd = FALSE, monitor = 1) {
  onewdata <- function(ext.input, modelnumber) {
    if (monitor >= 1) 
      cat("Number of clusters=", modelnumber, " simulation run ", ext.input, 
        "\n")
    cat(modelnumber, "class ", class(out$result$solution[[modelnumber]]), "\n")
    ndata <- generator.otrimle(X, out$result$solution[[modelnumber]])
    Grange <- modelnumber
    if (sim.est.logicd) 
      fixlogicd <- NULL
    else {
      fixlogicd <- rep(NA, max(Grange))
      for (g in G) if (!class(out$result$solution[[g]]) == "try-error") 
        fixlogicd[g] <- out$result$solution[[g]]$logicd
    }
    ond <- otrimleg(dataset = ndata$data, G = Grange, multicore = FALSE, erc = erc, 
      beta0 = beta0, fixlogicd = fixlogicd, monitor = monitor)
    ond
  }
  X <- as.matrix(dataset)
  p <- ncol(X)
  n <- nrow(X)
  out <- list()
  out$result <- otrimleg(dataset = X, G = G, multicore = multicore, ncores = ncores, 
    erc = erc, beta0 = beta0, monitor = monitor)
  out$simresult <- list()
  for (g in G) if (!class(out$result$solution[[g]]) == "try-error") {
    if (multicore) 
      out$simresult[[g]] <- mclapply(seq_len(simruns), onewdata, mc.cores = ncores, 
        modelnumber = g)
    else out$simresult[[g]] <- lapply(seq_len(simruns), onewdata, modelnumber = g)
  }
  out$simruns <- simruns
  class(out) <- "otrimlesimgdens"
  out
}
summary.otrimlesimgdens <- function(object, noisepenalty = 0.05, sdcutoff = 2, ...) {
  out <- list()
  out$G <- object$result$G
  out$simeval <- list()
  out$ssimruns <- object$simruns
  out$npr <- out$nprdiff <- out$logicd <- out$denscrit <- rep(NA, max(object$result$G))
  out$simeval$denscritmatrix <- matrix(NA, nrow = max(object$result$G), ncol = object$simruns)
  out$simeval$meandens <- out$simeval$sddens <- out$simeval$standens <- numeric(0)
  out$simeval$errors <- rep(NA, max(object$result$G))
  out$simeval$errors[object$result$G] <- 0
  out$denscrit <- object$result$denscrit
  for (i in object$result$G) {
    if (!class(object$result$solution[[i]]) == "try-error") {
      out$npr[i] <- object$result$solution[[i]]$exproportion[1]
      out$nprdiff[i] <- min(object$result$solution[[i]]$exproportion[2:(i + 
        1)]) - out$npr[i]
      out$logicd[i] <- object$result$solution[[i]]$logicd
      for (j in 1:object$simruns) {
        if (class(object$simresult[[i]][[j]]) == "try-error") {
          out$simeval$errors[i] <- out$simeval$errors[i] + 1
        }
        else {
          out$simeval$denscritmatrix[i, j] <- object$simresult[[i]][[j]]$denscrit[i]
        }
      }
      robustdens <- scaleTau2(out$simeval$denscritmatrix[i, !is.na(out$simeval$denscritmatrix[i, 
        ])], mu.too = TRUE)
      out$simeval$meandens[i] <- robustdens[1]
      out$simeval$sddens[i] <- robustdens[2]
      out$simeval$standens[i] <- (object$result$denscrit[i] - out$simeval$meandens[i])/out$simeval$sddens[i]
    }
  }
  out$peng <- object$result$G + (out$npr/noisepenalty)[object$result$G]
  lgv <- length(object$result$G)
  out$penorder <- order(out$peng)
  fcriterion <- out$simeval$standens[object$result$G]
  fcriterion[is.na(fcriterion)] <- 0
  g <- 1
  out$bestG <- object$result$G[which.min(fcriterion)]
  repeat {
    if (fcriterion[out$penorder[g]] <= sdcutoff) {
      out$bestG <- object$result$G[out$penorder[g]]
      break
    }
    g <- g + 1
    if (g > lgv) 
      break
  }
  out$sdcutoff <- sdcutoff
  out$bestresult <- object$result$solution[[out$bestG]]
  out$cluster <- object$result$solution[[out$bestG]]$cluster
  class(out) <- "summary.otrimlesimgdens"
  out
}
print.summary.otrimlesimgdens <- function(x, ...) {
  cat("Summary object of otrimlesimg (density version, useG=one)\n")
  cat("Numbers of clusters: ", x$bestresult$G, "\n")
  cat("Noise proportions: ", x$npr, "\n")
  cat("Penalised numbers of clusters: ", x$peng, "\n")
  cat("Order of penalised numbers of clusters: ", x$penorder, "\n\n")
  cat("Density criterion (small is good): ", x$denscrit, "\n")
  cat("Standardised density criterion (small is good; cutoff=", x$sdcutoff, "): ")
  for (i in x$G) {
    if (is.na(x$simeval$standens[i])) 
      cat("0/0 ")
    else cat(x$simeval$standens[i], " ")
  }
  cat("\n")
  if (x$bestG <= max(x$G)) 
    cat("Smallest fitting (best) number of clusters: ", x$bestG, "\n\n")
  else cat("No number of clusters <= ", max(x$bestresult$G), " was found to fit well.\n\n")
  cat("OTRIMLE criterion value: ", x$bestresult$criterion, "\n")
  cat("Difference between min cluster and noise proportion: ", x$nprdiff, "\n")
  cat("log ICD: ", x$logicd, "\n")
  cat("Simulation runs: ", x$ssimruns, "\n")
  cat("Errors in runs (by number of clusters): ", x$simeval$errors, "\n")
}
plot.summary.otrimlesimgdens <- function(x, plot = "criterion", penx = NULL, peny = NULL, 
  pencex = 1, cutoff = TRUE, ylim = NULL, ...) {
  opar <- par(no.readonly = TRUE)
  if (plot == "criterion") {
    if (is.null(ylim)) {
      rangevec <- c(x$denscrit, as.vector(x$simeval$denscritmatrix))
      ylim <- c(0, max(rangevec))
    }
    plot(x$G, x$denscrit, type = "l", col = 1, xlab = "Number of clusters", ylab = "Density criterion", 
      ylim = ylim, ...)
    for (i in x$G) {
      points(rep(i, x$ssimruns), x$simeval$denscritmatrix[i, ])
      if (cutoff) 
        points(i, x$simeval$meandens[i] + x$sdcutoff * x$simeval$sddens[i], 
          pch = "X", col = 2, cex = 2)
    }
  }
  else {
    par(mar = c(5, 4, 4, 4) + 0.3)
    if (is.null(ylim)) 
      ylim <- c(0, max(x$npr))
    plot(x$G, x$npr, type = "l", xlab = "Number of clusters", ylab = "Noise proportion", 
      ylim = ylim, ...)
    par(new = TRUE)
    plot(x$G, x$peng, type = "l", col = 2, axes = FALSE, xlab = "", ylab = "")
    par(col.axis = 2)
    axis(side = 4, at = pretty(range(x$peng)))
    mtext("Simplicity", side = 4, line = 3, col = 2)
    ylim <- range(x$peng)
  }
  if (is.null(penx)) 
    penx <- max(x$G)/4
  if (is.null(peny)) 
    peny <- 3 * max(ylim)/4
  pentext <- "Penalised order of n.c.:"
  for (q in x$penorder) pentext <- paste(pentext, " ", q)
  if (!identical(penx, FALSE)) 
    legend(penx, peny, pentext, cex = pencex)
  par(opar)
}
kmeanfun <- function(n) {
  out <- exp(-2.4562342) * n^{
    -0.3524332
  }
  out[n < 20] <- exp(-1.83084) * n[n < 20]^{
    -0.55323
  }
  out
}
ksdfun <- function(n) {
  out <- exp(-2.6575388) * n^{
    -0.4547216
  }
  out[n < 12.5] <- exp(-2.1482546) * n[n < 12.5]^{
    -0.4639833
  }
  out
}
kerndensmeasure <- function(x, weights = rep(1, nrow(as.matrix(x))), maxq = qnorm(0.9995), 
  kernn = 100) {
  sw <- sum(weights)
  wstan <- weights/sw
  x <- x - weighted.mean(x, wstan)
  kdn <- density(x, weights = wstan, from = -maxq, to = maxq, n = kernn)
  sk <- seq(1, (kernn - 1), by = 2)
  yor <- sort(kdn$y, decreasing = TRUE)
  yorm <- matrix(yor, nrow = 2)
  yorh <- colMeans(yorm)
  classcrit <- sqrt((sum((kdn$y[length(sk):1] - yorh)^2) + sum((kdn$y[(length(sk) + 
    1):kernn] - yorh)^2))/kernn)
  out <- list()
  out$cp <- yorh
  out$cpx <- kdn$y
  out$measure <- classcrit
  out
}
kerndensp <- function(x, weights = rep(1, nrow(as.matrix(x))), siglist, maxq = qnorm(0.9995), 
  kernn = 100) {
  p <- ncol(as.matrix(x))
  n <- nrow(as.matrix(x))
  nw <- sum(weights)
  psig <- princomp(x, covmat = siglist)
  out <- list()
  out$cml <- list()
  out$cm <- out$stanmeasure <- numeric(0)
  out$pca <- psig
  for (i in 1:p) {
    out$cml[[i]] <- kerndensmeasure(psig$scores[, i]/psig$sdev[i], weights, maxq = maxq, 
      kernn = kernn)
    out$cm[i] <- out$cml[[i]]$measure
    out$stanmeasure[i] <- (out$cm[i] - kmeanfun(nw))/ksdfun(nw)
  }
  out$measure <- sum(out$stanmeasure[out$stanmeasure > 0]^2)/p
  out
}
kerndenscluster <- function(x, fit, maxq = qnorm(0.9995), kernn = 100) {
  nc <- length(fit$exproportion) - 1
  ddpi <- out <- list()
  ddpm <- numeric(0)
  for (i in 1:nc) {
    siglist <- list(cov = fit$cov[, , i], center = fit$mean[, i], n.obs = nrow(x))
    ddpi[[i]] <- kerndensp(x, weights = fit$tau[, i + 1], siglist = siglist, 
      maxq = maxq, kernn = kernn)
    ddpm[i] <- ddpi[[i]]$measure
  }
  out$ddpi <- ddpi
  out$ddpm <- ddpm
  out$measure <- sqrt(sum(fit$exproportion[2:(nc + 1)] * ddpm^2))
  out
}
generator.otrimle <- function(data, fit) {
  U <- fit
  N <- nrow(data)
  P <- nrow(U$mean)
  G <- ncol(U$mean)
  pr <- U$exproportion
  data <- as.matrix(data)
  X <- matrix(0, nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
  cs <- rmultinom(n = 1, size = N, prob = pr)
  if (cs[1] >= 1) {
    X[1:cs[1], ] <- data[sample(1:N, size = cs[1], replace = TRUE, prob = U$tau[, 
      1]/sum(U$tau[, 1])), ]
  }
  for (j in 1:G) {
    if (cs[j + 1] > 0) {
      dataidx <- {
        sum(cs[1:j]) + 1
      }:{
        sum(cs[1:{
          j + 1
        }])
      }
      X[dataidx, ] <- rmvnorm(cs[j + 1], U$mean[, j], U$cov[, , j])
    }
  }
  ans <- list()
  ans$data <- X
  ans$clustering <- integer(0)
  for (i in 1:(G + 1)) if (cs[i] > 0) 
    ans$clustering <- c(ans$clustering, rep(i, cs[i]))
  return(ans)
}
