tvgarch <- function (y, order.g = 1, order.h = c(1,1,0), xtv = NULL, 
                     xreg = NULL, initial.values = list(), opt = 2, 
                     upper.speed = NULL, tvgarch = FALSE, turbo = FALSE, 
                     trace = FALSE)
{
  n <- length(y)
  y.index <- index(y)
  if (!is.null(order.g) && order.g == 0) order.g <- NULL 
  npar.h <- 1 + sum(order.h)
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    npar.h <- npar.h + NCOL(xreg)
  } 
  robse.par.h <- numeric(npar.h)
  maxpqr <- max(order.h)
  if (!is.null(order.g)) {
    if (!is.null(xtv)){
      xtv <- matrix(xtv, n, 1)
      if (is.null(colnames(xtv))) colnames(xtv) <- "xtv"
    }
    if (is.null(xtv)) {
      xtv <- matrix((1:n)/n, n, 1)
      colnames(xtv) <- "time"
      xtv <- as.zoo(xtv, order.by = y.index)
    }
    s <- length(order.g)
    npar.g <- 1 + 2 * s + sum(order.g)
  }
  '
    Initial values
  '
  if (!is.null(order.g)) {
    if (is.null(initial.values$intercept.g)) intercept.g <- 1
    else intercept.g <- initial.values$intercept.g
    if (is.null(initial.values$size)) size <- rep(0.1, s)
    else size <- initial.values$size
    if (is.null(initial.values$speed)) {
      speed <- rep(10, s)
      if (opt == 1) speed <- speed/sd(xtv)
      if (opt == 2) speed <- log(speed)
    }
    else speed <- initial.values$speed
    if (!is.null(upper.speed)) {
      if (length(upper.speed) != 1) {
        warning("Using only first upper bound for speed.")
        upper.speed <- upper.speed[1]
      }
      upper.speed <- rep(upper.speed, s)
    }
    if (is.null(upper.speed)) upper.speed <- rep(250, s)
    if (opt == 1) upper.speed <- upper.speed/sd(xtv)
    if (opt == 2) upper.speed <- log(upper.speed)
    if (is.null(initial.values$location)) {
      location <- NULL
      for (i in 1:s) {
        if (order.g[i] == 1) location <- c(location, mean(xtv))
        if (order.g[i] > 1) {
          location <- c(location, seq(from = min(xtv)+0.5*sd(xtv), 
                                      to = max(xtv)-0.5*sd(xtv), 
                                      by = (max(xtv)-min(xtv) - 
                                              sd(xtv))/(order.g[i]-1)))
        }
      }
    }
    else location <- initial.values$location
    if (s != length(size)) {
    stop("Mismatch between the number of transition 
         functions, s, and the number of initial size values.")
    }
    if (s != length(speed)) {
      stop("Mismatch between the number of transition functions, s, and the 
           number of initial speed values.")
    }
    if (sum(order.g) != length(location)) {
      stop("Mismatch between order.gand the initial location values.")
    }
  }
  if (is.null(initial.values$intercept.h)) intercept.h <- 0.1
  if (!is.null(initial.values$intercept.h)) {
    intercept.h <- initial.values$intercept.h
  }
  if (is.null(initial.values$arch)) arch <- rep(0.1/order.h[2], order.h[2])
  if (!is.null(initial.values$arch)) arch <- initial.values$arch
  if (order.h[2] != length(arch)) {
    stop("Mismatch between the number of ARCH-type parameters, q, and the arch 
         initial values.")
  }
  if (order.h[1] != 0) {
    if (is.null(initial.values$garch)) garch <- rep(0.7/order.h[1], order.h[1])
    else garch <- initial.values$garch
    if (order.h[1] != length(garch)) {
      stop("Mismatch between the number of GARCH-type parameters, p, and the 
           garch initial values.")
    }
  }
  if (order.h[1] == 0) garch <- NULL
  if (order.h[3] != 0) {
    if (is.null(initial.values$asym)) asym <- rep(0.02/order.h[3], order.h[3])
    if (!is.null(initial.values$asym)) asym <- initial.values$asym
    if (order.h[3] != length(asym)) {
      stop("Mismatch between the number of GJR-type parameters, r, and the asym 
           initial values.")
    }
  }
  if (order.h[3] == 0)  asym <- NULL
  if (!is.null(xreg)) {
    if (is.null(initial.values$par.xreg)) par.xreg <- rep(0.01, NCOL(xreg))
    if (!is.null(initial.values$par.xreg)) par.xreg <- initial.values$par.xreg
    if (NCOL(xreg) != length(par.xreg)) {
      stop("Mismatch between the number of covariates, X, and the par.xreg 
           initial values.")
    }
  }
  if (is.null(xreg)) par.xreg <- NULL
  '
    Vector of initial parameters for h component
  '
  names.h <- c("intercept.h", paste("arch", paste(seq(1:order.h[2]), sep = ""), 
                                    sep = ""))
  ini.par.h <- c(intercept.h, arch)
  if (order.h[1] != 0) {
    names.h <- c(names.h, paste("garch", paste(seq(1:order.h[1]), sep = ""), 
                                sep = "")) 
    ini.par.h <- c(ini.par.h, garch)
  }
  if (order.h[3] != 0) {
    names.h <- c(names.h, paste("asym", paste(seq(1:order.h[3]), sep = ""), 
                                sep = "")) 
    ini.par.h <- c(ini.par.h, asym)
  }
  if (!is.null(xreg)) { 
    if (!is.null(colnames(xreg))) names.h <- c(names.h, colnames(xreg))
    else names.h <- c(names.h, paste("xreg", paste(seq(1:NCOL(xreg)), sep = ""), 
                                     sep = ""))
    ini.par.h <- c(ini.par.h, par.xreg)
  }
  if (length(ini.par.h) != npar.h) {
    stop("Check initial parameters in the h component.")
  }
  initial.h <- matrix(ini.par.h, 1, npar.h)
  colnames(initial.h) <- names.h
  rownames(initial.h) <- "Value:"
  if (trace == TRUE) {
    cat("\n Initial parameters for the h component\n\n")
    print(round(initial.h, 4))
  }
  '
   g component and maximization by parts
  '
  g <- y
  g[1:n] <- 1
  if (!is.null(order.g)) {
    '
      Vector of initial parameters for g component
    ' 
    ini.par.g <- c(intercept.g, size, speed, location)
    names.g <- 
      c("intercept.g", paste("size", paste(seq(1:s), sep = ""), sep = ""), 
        paste("speed", paste(seq(1:s), sep = ""), sep = ""))
    for (i in 1:s) {
      if (s == 1) {
        names.g <- c(names.g, paste("location", 
                                    paste(seq(1:max(order.g[i])), sep = ""), 
                                    sep = ""))
      }
      else {
        names.g <- c(names.g, paste("location", paste(i, sep = ""),
                                    paste(seq(1:max(order.g[i])), sep = ""), 
                                    sep = ""))
      }
    }
    if (length(ini.par.g) != npar.g) {
      stop("Check initial parameters in the g component.")
    }
    initial.g <- matrix(ini.par.g, 1, npar.g)
    colnames(initial.g) <- names.g
    rownames(initial.g) <- "Value:"
    if (trace == TRUE) {
      cat("\n Initial parameters for the g component:\n\n")
      print(round(initial.g, 4))
    }
    '
      Parameter constraints for g component and iter = 0
    '
    r <- diag(npar.g)
    r1 <- r[1,]
    ui.g <- r1
    ci.g <- 1e-5
    if (s == 1) {
      r2 <- r[2,]
      r2[1] <- 1
      ui.g <- rbind(ui.g, r2)
      ci.g <- c(ci.g, 1e-5)
    } 
    if (s > 1) {
      r.size <- combos(s)
      r3 <- matrix(0, nrow(r.size), npar.g)
      r3[,1] <- 1
      r3[,2:(s+1)] <- r.size
      ui.g <- rbind(ui.g, r3)
      ci.g <- c(ci.g, rep(1e-5, nrow(r.size))) 
    }
    r4 <- r[(2+s):(2*s+1),]
    r5 <- r[(2*s+2):npar.g,]
    ui.g <- rbind(ui.g, r4, -r4, r5, -r5)
    ci.g <- c(ci.g, rep(1e-5, s), -upper.speed, rep(min(xtv)+1e-3,
                                                            sum(order.g)),
              rep(-max(xtv)+1e-3, sum(order.g))) 
    if (any(order.g > 1)) {
      for (i in which(order.g > 1)) {
        r6 <- -r[(1+2*s+sum(order.g[1:i]) - 
                    order.g[i]+1):(1+2*s+sum(order.g[1:i])-1),] + 
          r[(1+2*s+sum(order.g[1:i])-order.g[i]+2):(1+2*s+sum(order.g[1:i])),]
        ui.g <- rbind(ui.g, r6)
        ci.g <- c(ci.g, rep(0, order.g[i]-1)) 
      }
    }
    '
      Maximization by parts: estimating initial TV parameters
    '
    iter <- 0
    h <- rep(1, n)
    iter0.fit.g <- constrOptim(theta = ini.par.g, f = tvObj, grad = NULL, 
                               ui = ui.g, ci = ci.g, y = y, order.g = order.g, 
                               xtv = xtv, opt = opt, fixed.h = h, iter0 = TRUE, 
                               flag = 1)
    par.hat0.g <- iter0.fit.g$par
    g[1:n] <- tvObj(par.g = iter0.fit.g$par, fixed.par.g = NULL, xtv = xtv, 
                    opt = opt, order.g = order.g, fixed.h = h, y = y, 
                    iter0 = TRUE, flag = 2)
    estimates <- matrix(par.hat0.g, 1, length(names.g))
    colnames(estimates) <- names.g
    rownames(estimates) <- "Value:"
    if (trace == TRUE) {
      cat("\n Initial estimates for the g component:\n\n")
      if (turbo == TRUE) print(round(estimates[1,], 4))
      else print(round(estimates, 4))
    }
    '
      Parameter constraints for g component and iter > 0
    '
    par.hat.g <- par.hat0.g[-1]
    par.hat0.g <- par.hat0.g[1]
    r <- diag(length(par.hat.g))
    ui.g <- rbind(r[(s+1):(2*s),], -r[(s+1):(2*s),], r[(2*s+1):(npar.g-1),], 
                  -r[(2*s+1):(npar.g-1),])
    ci.g <- c(rep(1e-5, s), -upper.speed, rep(min(xtv)+1e-3, sum(order.g)), 
              rep(-max(xtv)+1e-3, sum(order.g))) 
    if (s == 1) {
      r1 <- r[1,]
      ui.g <- rbind(ui.g, r1)
      ci.g <- c(ci.g, -par.hat0.g[1]+1e-5)
    } 
    if (s > 1) {
      r.size <- combos(s)
      r1 <- matrix(0, nrow(r.size), length(par.hat.g))
      r1[,1:s] <- r.size
      ui.g <- rbind(ui.g, r1)
      ci.g <- c(ci.g, rep(-par.hat0.g[1]+1e-5, nrow(r.size))) 
    }
    if (any(order.g > 1)) {
      for(i in which(order.g > 1)){
        r2 <- 
          -r[(2*s+sum(order.g[1:i])-order.g[i]+1):(2*s+sum(order.g[1:i])-1),] + 
          r[(2*s+sum(order.g[1:i])-order.g[i]+2):(2*s+sum(order.g[1:i])),]
        ui.g <- rbind(ui.g, r2)
        ci.g <- c(ci.g, rep(0, order.g[i]-1)) 
      }
    }
    '
      Iterative estimation
    '
    iter <- iter + 1
    maxiter <- 1000 
    phi <- y/sqrt(g)
    par.hat.h <- ini.par.h
    repeat{
      if (trace == TRUE) {
        cat(" \n\n")
        cat("Iteration round:",iter, "\n\n")
      }
      '
        Estimating GJR-GARCH-X parameters:
      '
      iter.fit.h <- garchx(y = phi, order = order.h, xreg = xreg, 
                           initial.values = par.hat.h)
      h <- garchxRecursion(pars = as.numeric(iter.fit.h$par), aux = iter.fit.h)
      if (iter > 1) conv.h <- max(abs(par.hat.h-iter.fit.h$par))
      par.hat.h <- iter.fit.h$par
      if (trace == TRUE) {
        cat("Estimates (h component):", par.hat.h, "\n")
        cat("Likelihood:", logLik.garchx(iter.fit.h), "\n")
      }
      '
        Estimating TV parameters:
      '
      iter.fit.g <- constrOptim(theta = par.hat.g, f = tvObj, grad = NULL, 
                                ui = ui.g, ci = ci.g, fixed.par.g = par.hat0.g, 
                                y = y, order.g = order.g, xtv = xtv, opt = opt, 
                                fixed.h = h, iter0 = FALSE, flag = 1)
      conv.g <- max(abs(par.hat.g - iter.fit.g$par))
      par.hat.g <- iter.fit.g$par
      g[1:n] <- tvObj(par.g = par.hat.g, fixed.par.g = par.hat0.g, xtv = xtv, 
                      opt = opt, order.g = order.g, fixed.h = h, y = y, 
                      iter0 = FALSE, flag = 2)
      phi <- y/sqrt(g)
      if (trace == TRUE) {
        cat("Estimates (g component):", par.hat.g, "\n")
        cat("Likelihood:", - iter.fit.g$value, "\n")
        if (iter > 1) cat("Max. convergence:", max(conv.h, conv.g), "\n")
      }
      '
        Next round
      '
      iter <- iter + 1
      if (iter > 2) {
        if (iter >= maxiter || max(conv.g, conv.h) < 1e-5){
          if (iter == maxiter) {
            warning("Convergence has not been reached after ", 
                    paste(maxiter, sep = "")," iterations.")
          }
          break
        }
      }
    }
    par.hat.g <- c(par.hat0.g, iter.fit.g$par)
    if (turbo == TRUE) {
      vcov.g <- NULL
      robse.par.g <- rep(NA_real_, length(par.hat.g))
    } 
    if (turbo == FALSE) {
      '
          TV variance-covariance matrix:
      '
      vcov.g <- matrix(NA_real_, npar.g-1, npar.g-1)
      rownames(vcov.g) <- names.g[-1]
      colnames(vcov.g) <- names.g[-1]
      jac.g <- jacobian(func = tvObj, x = iter.fit.g$par, 
                        fixed.par.g = par.hat0.g, xtv = xtv, opt = opt, 
                        order.g = order.g, fixed.h = h, y = y, iter0 = FALSE, 
                        flag = 0)
      J.g <- crossprod(jac.g)  
      H.g <-  optimHess(par = iter.fit.g$par, fn = tvObj, 
                        fixed.par.g = par.hat0.g, xtv = xtv, opt = opt, 
                        order.g = order.g, fixed.h = h, y = y, iter0 = FALSE, 
                        flag = 1)
      solHG <- solve(-H.g)
      vcov.iter <- solHG %*% J.g %*% solHG
      vcov.g <- vcov.iter
      robse.par.g <- c(NA_real_, sqrt(diag(vcov.g)))
      rownames(vcov.g) <- names.g[-1]
      colnames(vcov.g) <- names.g[-1]
    }
    '
      TV final parameter estimates:
    '
    estimates.g <- as.matrix(rbind(par.hat.g, robse.par.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- names.g
  }
  '
    Final GJR-GARCH-X estimation
  '
  iter.fit.h <- garchx(y = y/sqrt(g), order = order.h, xreg = xreg, 
                       initial.values = ini.par.h)
  aux.h <- iter.fit.h
  h <- garchxRecursion(pars = as.numeric(iter.fit.h$par), aux = iter.fit.h)
  if (turbo == TRUE) {
    vcov.h <- NULL
    robse.par.h <- rep(NA_real_, length(iter.fit.h$par))
  }
  if (turbo == FALSE) {
    '
      GJR-GARCH-X variance-covariance matrix:
    '
    vcov.h <- vcov.garchx(iter.fit.h, vcov.type = "robust")
    robse.par.h <- sqrt(diag(vcov.h))
  }
  '
    GJR-GARCH-X final parameter estimates:
  '
  estimates.h <- as.matrix(rbind(iter.fit.h$par, robse.par.h))
  rownames(estimates.h) <- c("Estimate:", "Std. Error:")
  colnames(estimates.h) <- names.h
  '
    TV-GJR-GARCH-X final parameter estimates:
  '
  if (!is.null(order.g)) {
    par.ini <- c(estimates.g[1,-1],estimates.h[1,])
    r <- diag(length(par.ini))
    ui.gh <- cbind(ui.g, matrix(0, nrow(ui.g), npar.h))
    ci.gh <- ci.g
    ui.gh <- rbind(ui.gh, r[-(1:ncol(ui.g)),])
    ci.gh <- c(ci.gh, rep(0, npar.h)) 
    iter.fit <- constrOptim(theta = par.ini, f = tvgarchObj, grad = NULL, 
                            ui = ui.gh, ci = ci.gh, 
                            fixed.par.g = par.hat0.g[1], y = y, 
                            order.g = order.g, xtv = xtv, opt = opt, 
                            iter.fit.h = iter.fit.h, flag = 1)
    par <- iter.fit$par
    if (turbo == TRUE) {
      vcov <- NULL
      robse.par <- rep(NA_real_, length(iter.fit$par))
    }
    if (turbo == FALSE) {
      jac <- jacobian(func = tvgarchObj, x = iter.fit$par, 
                      fixed.par.g = estimates.g[1,1], y = y, order.g = order.g, 
                      xtv = xtv, opt = opt, iter.fit.h = iter.fit.h, flag = 0)
      JJ <- crossprod(jac)  
      HH <-  optimHess(par = iter.fit$par, fn = tvgarchObj, 
                       fixed.par.g = estimates.g[1,1], y = y, order.g = order.g, 
                       xtv = xtv, opt = opt, iter.fit.h = iter.fit.h, flag = 1)
      solHG <- solve(-HH)
      vcov <- solHG %*% JJ %*% solHG
      colnames(vcov) <- names(iter.fit$par)
      rownames(vcov) <- names(iter.fit$par)
      robse.par <- sqrt(diag(vcov))
    }
  }
  '
    TV-GJR-GARCH-X output:
  '
  par.h = estimates.h[1,]
  se.h = estimates.h[2,]
  if (!is.null(order.g)) {
    if (tvgarch == TRUE) {
      n <- length(xtv)  
      par.g <- c(estimates.g[1,1], par[1:(2*s+sum(order.g))])
      par.h <- par[-(1:(2*s+sum(order.g)))]
      G <- matrix(1, n, s)
      for (i in 1:s) {
        G[,i] <- 
          tv(speed = par.g[1+s+i], 
             location = par.g[(1+2*s+sum(order.g[1:i])-
                                 order.g[i]+1):(1+2*s+sum(order.g[1:i]))], 
             xtv = xtv, opt = opt, order.g = order.g[i])
      }
      g <- par.g[1] + G %*% par.g[2:(s+1)]
      h <- garchxRecursion(pars = as.numeric(par.h), aux = iter.fit.h)
      se.g <- c(NA_real_, robse.par[1:(2*s+sum(order.g))])
      se.h <- robse.par[-(1:(2*s+sum(order.g)))]
    }
    else {
      par.g = estimates.g[1,]
      se.g = estimates.g[2,]
    }
  }
  sigma2 <- as.matrix(h*g)
  colnames(sigma2) <- "sigma2"
  logLik <- sum(dnorm(x = y, mean = 0, sd = sqrt(sigma2), log = TRUE))
  residuals <- as.matrix(y/sqrt(sigma2))
  colnames(residuals) <- "innovations"
  results <- list(par.h = par.h, se.h = se.h, 
                  names.h = names.h, sigma2 = sigma2, residuals = residuals, 
                  h = h, g = g, logLik = logLik, vcov.h = vcov.h, 
                  message.h = iter.fit.h$message, order.g = order.g, 
                  order.h = order.h, xreg = xreg, y = y, y.index = y.index, 
                  date = date(), turbo = turbo, aux.h = aux.h, 
                  iter.fit.h = iter.fit.h, tvgarch = tvgarch)
  if (!is.null(order.g)) {
    results$vcov = vcov
    results$par.g = par.g
    results$se.g = se.g
    results$names.g = names.g
    results$vcov.g = vcov.g
    results$message.g = iter.fit.g$message
    results$xtv = xtv
    results$opt = opt
    results$iter = iter
  }
  class(results) <- "tvgarch"
  return(results)
} 
 


