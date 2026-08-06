mtvgarchSim <- function (n, m = 2, order.g = c(1,1), order.h = c(1,1,0, 1,1,0), 
                         order.x = NULL, intercept.g = c(1.2,1), size = c(3,5), 
                         speed = c(10,25), location = c(0.5,0.8), 
                         intercept.h = c(0.2,0.3), arch = c(0.10,0.05),
                         garch = c(0.80,0.90), asym = NULL, xtv = NULL, 
                         xreg = NULL, par.xreg = NULL, R = c(1,0.6,0.6,1), 
                         dcc = FALSE, par.dcc = NULL, opt = 0, as.zoo = TRUE, 
                         verbose = FALSE, innovations = NULL)
{
  if (m == 1) stop("For univariate processes use tvgarchSim().")
  names.ID <- paste("y", 1:m, sep = "")
  if (dcc == TRUE && is.null(par.dcc)){
    stop("dcc = TRUE but par.dcc = NULL.")
  }
  if (!is.null(order.g)) {
    if (!any(order.g != 0)) order.g <- NULL
    else {
      if (any(is.null(c(intercept.g, size, speed, location)))) {
      stop ("order.g neither zero nor NULL but intercept.g, size, speed and/or 
      location are/is missing.")
    }
      order.g <- as.matrix(order.g)
      rownames(order.g) <- names.ID
      max.s <- ncol(order.g)
      if (sum(order.g) > length(location[location != 0])) {
        stop("Not enough location parameters.")
      }
      colnames(order.g) <- paste("G", 1:max.s, sep = "")
      names.g <- c("intercept", paste("size", 1:max.s, sep = ""), 
                   paste("speed", 1:max.s, sep = ""))
      for (j in 1:max.s) {
        names.g <- c(names.g, paste("location", j, 1:max(order.g[,j]), sep = ""))
      }
      par.g <- matrix(NA, m, 1+2*max.s+max(rowSums(order.g)))
      colnames(par.g) <- names.g
      rownames(par.g) <- names.ID
      intercept.g <- matrix(intercept.g, m, 1)
      size <- matrix(size, m, max.s, byrow = TRUE)
      speed <- matrix(speed, m, max.s, byrow = TRUE)
      location <- matrix(location, m, max(rowSums(order.g)), byrow = TRUE)
    }
  }
  if (is.null(order.g)) par.g <- NULL
  order.h <- matrix(order.h, m, 3, byrow = TRUE)
  if (any(order.h[,2] == 0)) {
    order.h[which(order.h[,2] == 0),2] = 1
    warning("Models with no ARCH effects (order.h[,2] = 0) not allowed. 
             order.h[,2] replaced by 1. ")
  }
  colnames(order.h) <- c("GARCH", "ARCH", "ASYM")
  rownames(order.h) <- names.ID
  max.p <- max(order.h[,1])
  max.q <- max(order.h[,2])
  max.r <- max(order.h[,3])
  if (sum(order.h[,1]) != length(garch[garch != 0])) {
    stop("Argument garch or order.h[,1] needs revision.")
  }
  if (sum(order.h[,2]) != length(arch[arch != 0])) {
    stop("Argument arch or order.h[,2] needs revision.")
  }
  if (sum(order.h[,3]) != length(asym[asym != 0])) {
    stop("Argument asym or order.h[,3] needs revision.")
  }
  if (!is.null(order.x) && is.null(xreg)) {
    spillovers <- TRUE
    if (any(c(max.q, max.q, max.q) > 1)) {
      stop("Maximum order.h = (1,1,1) with volatility spillovers when simulating 
           variances.")
    }
  }
  else spillovers <- FALSE
  intercept.h <- matrix(intercept.h, m, 1)
  par.h <- intercept.h
  names.h <- c("intercept", paste("arch", 1:max.q, sep = ""))
  arch <- matrix(arch, m, max.q, byrow = TRUE)
  par.h <- cbind(par.h, arch)
  if (max.p != 0) {
    names.h <- c(names.h, paste("garch", 1:max.p, sep = ""))
    garch <- matrix(garch, m, max.p, byrow = TRUE)
    garch[order.h[,1] == 0] <- 0
    par.h <- cbind(par.h, garch)
  }
  if (max.r != 0) {
    names.h <- c(names.h, paste("asym", 1:max.r, sep = ""))  
    asym <- matrix(asym, m, max.r, byrow = TRUE)
    par.h <- cbind(par.h, asym)
  }
  if (!is.null(order.x)) {
    if (is.null(par.xreg)) {
      stop ("order.x not NULL but par.xreg is missing.")
    }
    if (sum(order.x) != length(par.xreg[par.xreg != 0])) {
      stop("Number of par.xreg unequal to sum(order.x).")
    }
    order.x2 <- order.x
    order.x2[order.x2 == 0] <- 1
    if (is.null(xreg)) {
      if (sum(order.x2) != m^2) {
        stop ("order.x not binary.")
      }
      order.x <- matrix(order.x, m, m, byrow = TRUE)
      par.xreg2 <- par.xreg
      par.xreg <- matrix(0, m, m)  
      par.xreg[which(order.x == 1)] <- par.xreg2  
      diag(order.x) <- 0
      names.h <- c(names.h, paste("lag(", names.ID, "^2)", sep = ""))  
    } 
    if (!is.null(xreg)) {
      if (length(order.x) != m*NCOL(xreg)) {
        stop("Length of order.x unequal to m*NCOL(xreg).")
      }
      if (sum(order.x2) != m*NCOL(xreg)) {
        stop ("order.x not binary.")
      }
      order.x <- matrix(order.x, m, NCOL(xreg), byrow = TRUE)
      par.xreg2 <- par.xreg
      par.xreg <- matrix(0, m, NCOL(xreg)) 
      par.xreg[which(order.x == 1)] <- par.xreg2  
      names.h <- c(names.h, paste("xreg", 1:NCOL(xreg), sep = ""))  
    }
    par.h <- cbind(par.h, par.xreg)
  }
  colnames(par.h) <- names.h
  rownames(par.h) <- names.ID
  if (is.null(innovations)) {
    innovations <- matrix(rnorm(n*m), n, m, byrow = TRUE)
  }
  '
    Constructing correlated innovations
  '
  R <- matrix(R, m, m, byrow = TRUE)
  if (sum(diag(R)) != m || any(R > 1) || any(R < (-1)) || 
      min(eigen(R)$values) < 0) {
    stop("R is not a valid correlation matrix.")
  }
  ID <- 1:m
  namesR <- numeric(0)
  for (i in 1:m) {          
    for (j in 1:m) {    
      namesR <- c(namesR, paste("r", paste(paste(ID[j]), ID[i], sep = ""), 
                                sep = ""))
    }
  }
  namesR <- matrix(namesR, m, m)
  namesR <- namesR[lower.tri(namesR)]
  '
    Constructing constant conditional correlations (ccc)
  '
  if (dcc == FALSE) { 
    chol_P <- chol(R)
    innovations <- innovations %*% chol_P  
    ccc <- matrix(R[lower.tri(R)], n, m*(m-1)/2)
    colnames(ccc) <- namesR
  }
  '
    Constructing dynamic conditional correlations (dcc)
  '
  if (dcc == TRUE) {   
    dcc <- matrix(0, n, m*(m-1)/2)
    Qbar <- (1 - par.dcc[1] - par.dcc[2]) * R
    Qt <- R  
    zt2 <- R
    for (t in 1:n) {
      Qt <- Qbar + par.dcc[1] * zt2 + par.dcc[2] * Qt
      invQt <- diag( 1 / sqrt( diag( Qt ) ) )
      Rt  <- invQt %*% Qt %*% invQt
      dcc[t,] <- Rt[lower.tri(Rt)]
      cholRt <- t(chol(Rt))
      zt <- drop(cholRt %*% innovations[t,])
      zt2 <- zt %o% zt
      innovations[t,] <- zt
    }
    colnames(dcc) <- namesR
  }
  if (is.null(colnames(innovations))) {
    colnames(innovations) <- paste("innovations", 1:m, sep = "")
  }
  '
    Constructing g component 
  '
  if (!is.null(order.g)) { 
    if (is.null(xtv)) xtv <- (1:n) / n
    g <- matrix(1, n, m)
    for (i in 1:m) {
      if (order.g[i,1] != 0) {
        s.i <- length(which(order.g[i,] > 0))
        if (!is.null(intercept.g)) intercept.g.i <- intercept.g[i]
        if (!is.null(size)) size.i <- size[i,1:s.i]
        if (!is.null(speed)) speed.i <- speed[i,1:s.i]
        if (!is.null(location)) location.i <- location[i,1:sum(order.g[i,])]
        G <- matrix(1, n, s.i)
        for (j in 1:s.i) {
          G[,j] <- tv(speed = speed.i[j], 
                      location = location.i[(sum(order.g[i,1:j]) -
                                          order.g[i,j]+1):sum(order.g[i,1:j])], 
                      xtv = xtv, opt = opt, order.g = order.g[i,j])
        }
        g[,i] <- intercept.g.i[1] + G %*% size.i[1:s.i]
        par.g[i,1] <- intercept.g.i
        par.g[i,(1+1):(1+length(order.g[i,]))] <- size.i
        par.g[i,(1+length(order.g[i,])+1):(1+2*length(order.g[i,]))] <- speed.i
        par.g[i,(1+2*length(order.g[i,])+1):(1+2*length(order.g[i,]) +
                                               sum(order.g[i,]))] <- location.i
      }
    }
    colnames(g) <- paste("g", 1:m, sep = "")
  }
  if (is.null(order.g)) g <- matrix(1, n, m)
  '
    Constructing h component
  '
  h <- matrix(1, n, m)
  if (is.null(order.x)) xreg.i <- NULL
  if (spillovers == FALSE) {
    for (i in 1:m) {
      intercept.h.i <- intercept.h[i]
      arch.i <- arch[i,1:order.h[i,2]]
      if (order.h[i,1] != 0) garch.i <- garch[i,1:order.h[i,1]]
      if (order.h[i,1] == 0) garch.i <- NULL
      if (order.h[i,3] != 0) asym.i <- asym[i,1:order.h[i,3]]
      if (order.h[i,3] == 0) asym.i <- NULL
      if (!is.null(order.x)) {
        if (any(order.x[i,] != 0)) xreg.i <- 
            as.matrix(xreg[,which(order.x[i,] == 1)]) %*% 
            par.xreg[i,which(order.x[i,] == 1)]
        if (sum(order.x[i,] == 0)) xreg.i <- NULL
      }
      h[,i] <- garchxSim(n = n, intercept = intercept.h.i, arch = arch.i, 
                         garch = garch.i, asym = asym.i, xreg = xreg.i, 
                         innovations = innovations[,i], as.zoo = FALSE, 
                         verbose = TRUE)[,"sigma2"]
      par.h[i,1] <- intercept.h.i
      par.h[i,(1+1):(1+order.h[i,2])] <- arch.i
      if (any(order.h[,1] != 0)) { 
        if (order.h[i,1] != 0) {
          par.h[i,(1+max.q+1):(1+max.q+order.h[i,1])] <- garch.i
        }
        if (order.h[i,1] == 0) {
          par.h[i,(1+max.q+1):(1+max.q+max.p)] <- 0
        }
      }
      if (any(order.h[,3] != 0)) { 
        if (order.h[i,3] != 0) {
          par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+order.h[i,3])] <- asym.i
        }
        if (order.h[i,3] == 0) {
          par.h[i,(1+max.q+max.p+1):(1+max.q+max.p+max.r)] <- 0
        }
      }
      if (!is.null(order.x)) {
        par.h[i,1+max.q+max.p+max.r+which(order.x[i,] == 1)] <- 
          par.xreg[i,which(order.x[i,] == 1)]
      }
    }
  }
  if (spillovers == TRUE) {
    max.x <- ncol(order.x)
    ARCH <- matrix(0, m, m)
    diag(ARCH) <- arch[1:m,1]
    ARCH[which(order.x == 1)] <- par.xreg[which(order.x == 1)]
    GARCH <- matrix(0, m, m)
    ASYM <- matrix(0, m, m)
    if (max.p == 1) {
      diag(GARCH) <- garch[1:m,1]
    }
    if (max.r == 1) {
      diag(ASYM) <- asym[1:m,1]
    }
    ht <- intercept.h/(1-diag(ARCH)-diag(GARCH))
    if (max(Mod(eigen(ARCH+GARCH)$values)) > 1) {
      stop("The (initial) stationarity condition is not satisfied.")
    }
    h[1,] <- ht
    y2 <- colMeans(innovations^2)
    Ineg <- innovations < 0          
    Ineg.y2 <- rep(0, m)   
    for (t in 2:n) {
      ht <- intercept.h + ARCH%*%y2 + GARCH%*%ht + ASYM%*%Ineg.y2
      h[t,] <- ht
      y2 <- h[t,]*(innovations[t,]^2)
      Ineg.y2 <- Ineg[t,]*y2
    }
  }
  colnames(h) <- paste("h", 1:m, sep = "")
  '
    Output
  '
  sigma2 <- h*g
  colnames(sigma2) <- paste("sigma2", 1:m, sep = "")
  y <- sqrt(sigma2)*innovations
  colnames(y) <- names.ID
  if (spillovers == TRUE) {
    names.ID2 <- paste("phi", 1:m, sep = "")
    m.y2 <- colMeans(y^2)
    xreg <- rbind(m.y2, y[-n,]^2)
    colnames(xreg) <- paste("lag(", names.ID2, ",^2, 1)", sep = "")
  }
  if (as.zoo == TRUE) {
    y <- zoo(y, order.by = index(y))
    sigma2 <- zoo(sigma2, order.by = index(y))
    innovations <- zoo(innovations, order.by = index(y))
    g <- zoo(g, order.by = index(y))
    h <- zoo(h, order.by = index(y))
    xreg <- zoo(xreg, order.by = index(y))
  }
  if (verbose == TRUE) {
    if (any(colSums(par.h) == 0)) {
      par.h = par.h[,-which(colSums(par.h) == 0)]
    }
    result <- list(y = y, innovations = innovations, sigma2 = sigma2, 
                   spillovers = spillovers, g = g, order.g = order.g, 
                   par.g = par.g, h = h, order.h = order.h, 
                   par.h = par.h, xreg = xreg)
    if (is.null(par.dcc)) {
      result$ccc <- ccc
      if (as.zoo == TRUE) result$ccc <- zoo(result$ccc, order.by = index(y))
    }
    if (!is.null(par.dcc)) {
      result$dcc <- dcc
      if (as.zoo == TRUE) result$dcc <- zoo(result$dcc, order.by = index(y))
    }
  }
  if (verbose == FALSE) 
  {
    result <- y
  }
  return(result)
}
