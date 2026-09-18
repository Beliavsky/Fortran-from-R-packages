.OtrimleCols <- function(...) {
  Cols <- c("#000000", "#E0162B", "#0052A5", "#31A853", "#F48024", "#8E44AD", "#437DCC", 
    "#F5BD2F", "#566473", "#DC143C", "#915C83", "#CDBFAC", "#87A96B", "#6E7F80", 
    "#00B1C4", "#FFE135", "#BCD4E6", "#DE5D83", "#CC0000", "#006B3C", "#92A1CF", 
    "#2A52BE", "#8C92AC", "#FFBCD9", "#E52B50")
  return(Cols)
}
.PlotOtrimleClustering <- function(data, cluster, margins = {
  1:ncol(data)
}, sym.noise = "+", pictitle = "") {
  palette(.OtrimleCols())
  if (min(cluster) == 0) {
    pch <- c(0:9, letters, LETTERS)[cluster + 1]
  }
  else {
    pch <- c(1:9, letters, LETTERS)[cluster]
  }
  pch[cluster == "0"] <- "+"
  clustcol <- cluster + 1
  if (ncol(data) == 1) {
    plot(cbind(data, 2), col = clustcol, pch = pch, yaxt = "n", ylab = "", xlab = "data", 
      ylim = c(1, 3), main = pictitle)
    axis(side = 1)
  }
  else if (length(margins) == 2) {
    if (is.null(colnames(data))) {
      xlab <- paste("data[ ,", margins[1], "]", sep = "")
      ylab <- paste("data[ ,", margins[2], "]", sep = "")
    }
    else {
      xlab <- colnames(data)[margins[1]]
      ylab <- colnames(data)[margins[2]]
    }
    plot(data[, margins], pch = pch, col = clustcol, main = pictitle, xlab = xlab, 
      ylab = ylab, cex.main = 1)
  }
  else if (length(margins) > 2) {
    pairs(data[, margins], pch = pch, col = clustcol, main = pictitle, cex.main = 1)
  }
}
.PlotOtrimleFit <- function(x, cluster = NULL) {
  palette(.OtrimleCols())
  G <- dim(x$cov)[3]
  P <- dim(x$cov)[2]
  if (!{
    is.null(cluster) || cluster %in% {
      1:G
    }
  }) {
    stop("\"cluster\" is equal to \"NULL\" or an integer from {1,2,...,number of clusters} ")
  }
  if (is.null(cluster) & G > 1) {
    if (G == 2) {
      nr <- 1
      nc <- 2
    }
    else if (G > 2) {
      nc <- nr <- ceiling(sqrt(G))
    }
    par(mfrow = c(nr, nc))
    for (j in 1:G) {
      plot(x = .wecdf(x = x$smd[, j], weights = x$tau[, 1 + j]), y = pchisq(x$smd[, 
        j], df = P), col = j + 1, pch = 20, main = paste("Cluster ", j, sep = ""), 
        xlab = "Weighted ECDF of SMD", ylab = expression(paste("Prob(SMD) under ", 
          {
          chi[p]
          }^2)), cex.main = 1)
      abline(coef = c(0, 1))
    }
    par(mfrow = c(1, 1))
  }
  else {
    plot(x = .wecdf(x = x$smd[, cluster], weights = x$tau[, 1 + cluster]), y = pchisq(x$smd[, 
      cluster], df = P), col = cluster + 1, pch = 20, main = paste("Cluster ", 
      cluster, sep = ""), xlab = "Weighted ECDF of SMD", ylab = expression(paste("Prob(SMD) under ", 
      {
        chi[p]
      }^2)), cex.main = 1)
    abline(coef = c(0, 1))
  }
}
.PlotOtrimleOptimization <- function(x) {
  palette(.OtrimleCols())
  ok <- !is.na(x$optimization[, 2])
  X <- x$optimization[ok, 1]
  Y <- x$optimization[ok, 2]
  if (length(X) == 1) {
    cat("\n")
    cat("OTRIMLE object contains a single solution at logicd=", x$logicd, ", see deails:\n\n", 
      sep = "")
    print(x$optimization)
    cat("\n")
  }
  else {
    if (X[1] == -Inf) {
      if (length(X[-1]) == 1) {
        infval <- X[-1] - 0.3 * abs(X[-1])
      }
      else {
        infval <- min(X[-1]) - 0.3 * max(diff(X[-1]))
      }
      X[1] <- infval
      x.lim <- c(infval, max(X))
      y.lim <- range(Y)
      plot(X[-1], Y[-1], t = "b", xlim = x.lim, ylim = y.lim, pch = 20, lwd = 2, 
        col = 3, main = "OTRIMLE: criterion profiling", axes = FALSE, frame = TRUE, 
        xlab = "log(icd)", ylab = "Criterion", cex.main = 1)
      x.ticks <- axTicks(1)
      y.ticks <- axTicks(2)
      axis(1, las = 1, at = c(X[1], x.ticks[-1]), labels = c(expression(paste("-", 
        infinity)), x.ticks[-c(1)]))
      axis(2, las = 0, at = y.ticks, labels = y.ticks)
      points(X[1], Y[1], pch = 20, col = 2, cex = 1)
      segments(x0 = X[1], y0 = 0, x1 = X[1], y1 = Y[1], col = 2, lty = 2)
    }
    else {
      plot(X, Y, t = "b", pch = 20, lwd = 2, col = 3, cex = 1, main = "OTRIMLE: criterion profiling", 
        xlab = "log(icd)", ylab = "Criterion", cex.main = 1)
    }
  }
}
.PlotOtrimleIloglik <- function(x) {
  palette(.OtrimleCols())
  ok <- !is.na(x$optimization[, 2])
  X <- x$optimization[ok, 1]
  Y <- x$optimization[ok, 3]
  if (length(X) == 1) {
    cat("\n")
    cat("OTRIMLE object contains a single solution at logicd=", x$logicd, ", see deails:\n\n", 
      sep = "")
    print(x$optimization)
    cat("\n")
  }
  else {
    if (X[1] == -Inf) {
      if (length(X[-1]) == 1) {
        infval <- X[-1] - 0.3 * abs(X[-1])
      }
      else {
        infval <- min(X[-1]) - 0.3 * max(diff(X[-1]))
      }
      X[1] <- infval
      x.lim <- c(infval, max(X))
      y.lim <- range(Y)
      plot(X[-1], Y[-1], t = "b", xlim = x.lim, ylim = y.lim, pch = 20, lwd = 2, 
        col = 3, main = "OTRIMLE: improper log-likelihood profiling", axes = FALSE, 
        frame = TRUE, xlab = "log(icd)", ylab = "Improper log-likelihood", 
        cex.main = 1)
      x.ticks <- axTicks(1)
      y.ticks <- axTicks(2)
      axis(1, las = 1, at = c(X[1], x.ticks[-1]), labels = c(expression(paste("-", 
        infinity)), x.ticks[-c(1)]))
      axis(2, las = 0, at = y.ticks, labels = y.ticks)
      points(X[1], Y[1], pch = 20, col = 2, cex = 1)
      segments(x0 = X[1], y0 = 0, x1 = X[1], y1 = Y[1], col = 2, lty = 2)
    }
    else {
      plot(X, Y, t = "b", pch = 20, lwd = 2, col = 3, cex = 1, main = "OTRIMLE: improper log-likelihood profiling", 
        xlab = "log(icd)", ylab = "Improper log-likelihood", cex.main = 1)
    }
  }
}
plot.otrimle <- function(x, what = c("criterion", "iloglik", "fit", "clustering"), 
  data = NULL, margins = NULL, cluster = NULL, ...) {
  if (!inherits(x, "otrimle")) {
    stop("x is not of class \"otrimle\"")
  }
  if (x$code == 0) {
    stop("OTRIMLE object does not contain a solution (x$code=0)")
  }
  if (length(what) > 1) {
    what <- what[1]
  }
  if (what == "clustering") {
    if (is.null(data)) {
      stop("argument \"data\" (vector, matrix or data.frame) is required for  what=\"clustering\"")
    }
    if (is.vector(data)) {
      data <- matrix(data, ncol = 1)
    }
    if (!is.matrix(data)) {
      data <- data.matrix(data)
    }
    if (!is.numeric(data)) {
      stop("\"data\" must be numeric.")
    }
    if (nrow(data) != length(x$cluster)) {
      stop("nrow(data) is not equal to length(cluster)")
    }
    if (is.null(margins)) {
      margins <- 1:ncol(data)
    }
    .PlotOtrimleClustering(data = data, cluster = x$cluster, margins = margins, 
      pictitle = "OTRIMLE: clustering")
  }
  else if (what == "criterion") {
    .PlotOtrimleOptimization(x)
  }
  else if (what == "iloglik") {
    .PlotOtrimleIloglik(x)
  }
  else if (what == "fit") {
    .PlotOtrimleFit(x, cluster = cluster)
  }
}
plot.rimle <- function(x, what = c("fit", "clustering"), data = NULL, margins = NULL, 
  cluster = NULL, ...) {
  if (!inherits(x, "rimle")) {
    stop("x is not of class \"rimle\"")
  }
  if (x$code == 0) {
    stop("RIMLE object does not contain a solution (object$code=0)")
  }
  if (length(what) > 1) {
    what <- what[1]
  }
  if (what == "clustering") {
    if (is.null(data)) {
      stop("argument \"data\" (vector, matrix or data.frame) is required for  what=\"clustering\"")
    }
    if (is.vector(data)) {
      data <- matrix(data, ncol = 1)
    }
    if (!is.matrix(data)) {
      data <- data.matrix(data)
    }
    if (!is.numeric(data)) {
      stop("\"data\" must be numeric.")
    }
    if (nrow(data) != length(x$cluster)) {
      stop("nrow(data) is not equal to length(cluster)")
    }
    if (is.null(margins)) {
      margins <- 1:ncol(data)
    }
    .PlotOtrimleClustering(data = data, cluster = x$cluster, margins = margins, 
      pictitle = "RIMLE: clustering")
  }
  else if (what == "fit") {
    .PlotOtrimleFit(x, cluster = cluster)
  }
}
