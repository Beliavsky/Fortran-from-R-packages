## 0. Auxilary functions ##################################


#' Transform a vector of correlations to symmetric correlation matrix
#' @param rho vector of correlations
#' @return symmetric matrix of correlations
#' @author Erik Hintz with code from Yoshiba (2018)
rhoToOmega <- function(rho){
   dim <- (sqrt(8*length(rho)+1)+1)/2 
   Omega <- diag(1/2,dim)
   Omega[lower.tri(Omega)] <- rho
   Omega <- Omega + t(Omega)
   Omega
}

## 1.1. t copula case ##########################################################

#' Transform internal to original parameters (t case)
#' @param para vector of parameters 
#' @return 2-list with the correlation matrix and dof 'df'  
#' @author Erik Hintz with code from Yoshiba (2018)
int_to_org_pars_t <- function(par){
   lengthpar <- length(par)
   df <- par[lengthpar] # grab dof (last element)
   theta <- par[-lengthpar] 
   dim <- (1 + sqrt(1 + 8*(lengthpar - 1)))/2
   LTR <- diag(dim)
   LTR[-1,1] <- cos(theta[1:(dim-1)])
   cumsin <- sin(theta[1:(dim-1)])
   if(dim > 2){
      for(j in 2:(dim-1)){
         LTR[j,j] <- cumsin[1]
         k <- (j - 1)*(dim - j/2) + 1
         thj <- theta[k:(k+dim-j-1)]
         cumsin <- cumsin[-1]
         LTR[((j+1):dim),j] <- cumsin * cos(thj)
         cumsin <- cumsin * sin(thj) 
      } 
   }
   LTR[dim,dim] <- cumsin[1]
   Omega <- LTR %*% t(LTR)
   ## Return
   list(rho = Omega[lower.tri(Omega)], df = df)
}

#' Transform original parameters to internal parameters (t case)
#' @param rho original 'correlations'
#' @param gamma skewness 
#' @param df degrees of freedom
#' @return vector with transformed parameters (ANGLES -- GAMMA -- DF)
#' @author Erik Hintz with code from Yoshiba (2018)
org_to_int_pars_t <- function(rho, df){
   R <- rhoToOmega(rho)
   LTR <- t(chol(R))
   dim <- nrow(LTR)
   theta <- acos(LTR[2:dim, 1])
   cumsin <- sin(theta)[-1]
   if(dim > 2){
      for(j in 2:(dim-1)){
         thj <- acos(LTR[(j+1):dim,j]/cumsin);
         theta <- c(theta,thj);
         cumsin <- (cumsin*sin(thj))[-1];
      } 
   }
   ## Return
   c(theta, df)
}


## 1.2. skew-t copula case #####################################################

#' Transform internal to original parameters (skew t case)
#' @param par vector of parameters 
#' @return 3-list with the correlation matrix, dof 'df' and 'gamma' 
#' @author Erik Hintz with code from Yoshiba (2018)
int_to_org_pars_st <- function(par){
   ntheta <- length(par) - 2 
   dim <- (1 + sqrt(1 + 8*ntheta))/2
   theta <- par[1:ntheta]
   LTR <- diag(dim)
   LTR[-1,1] <- cos(theta[1:(dim-1)])
   cumsin <- sin(theta[1:(dim-1)])
   if(dim > 2){
      for(j in 2:(dim-1)){
         LTR[j,j] <- cumsin[1]
         k <- (j - 1)*(dim - j/2) + 1
         thj <- theta[k:(k+dim-j-1)]
         cumsin <- cumsin[-1]
         LTR[((j+1):dim),j] <- cumsin*cos(thj)
         cumsin <- cumsin * sin(thj) 
      } 
   }
   LTR[dim,dim] <- cumsin[1]
   Omega <- LTR %*% t(LTR)
   gamma <- par[ntheta+2] 
   df <- exp(par[ntheta+1]) + 2 
   ## Return
   list(rho = Omega[lower.tri(Omega)], df = df, gamma = gamma)
}

#' Transform original parameters to internal parameters (skew-t case)
#' @param rho original 'correlations'
#' @param gamma skewness 
#' @param df degrees of freedom
#' @return vector with transformed parameters (ANGLES -- GAMMA -- DF)
#' @author Erik Hintz with code from Yoshiba (2018)
org_to_int_pars_st <- function(rho, gamma, df){
   R <- rhoToOmega(rho)
   LTR <- t(chol(R))
   dim <- nrow(LTR)
   theta <- acos(LTR[2:dim, 1])
   cumsin <- sin(theta)[-1]
   if(dim > 2){
      for(j in 2:(dim-1)){
         thj <- acos(LTR[(j+1):dim,j]/cumsin);
         theta <- c(theta,thj);
         cumsin <- (cumsin*sin(thj))[-1];
      } 
   }
   ## Return
   c(theta, log(df - 2), gamma)
}


## 2. Fitting functions ##################################################

#' Fit P matrix of skew-t copula using an EM-like algorithm
#' @param u (n, d) data-matrix (must be in (0,1))
#' @param pseudoskewt (n, d) matrix of pseudo-observations 
#' @param df dof parameter
#' @param gamma skewness vector
#' @param P starting guess for P
#' @param P_inv inverse of P 
#' @param ldet log(det(P))
#' @param report.ll logical if log-lik before and after should be reported
#' @param P_maxiter maximum number of iterations for 'P' updates
#' @param P_tol relative convergence tolerance for elements in 'P'
#' @return 5-list with 'P_next', 'P_next_inv', 'ldet_next', 'll' 
#' @author Erik Hintz 
fitscaleskewtEM <- function(u, pseudoskewt = NULL, df, gamma, P, P_inv = NULL,
                            ldet = NULL,
                            report.ll = FALSE, P_maxiter = 100, P_tol, ...){
   
   if(!is.matrix(u)) u <- cbind(u)
   d <- ncol(u)
   n <- nrow(u)
   ## If length(gamma) == 1 => equiskewed; make a d-vector for mahalanobis() 
   if(length(gamma) != d){
      stopifnot(length(gamma) == 1)
      gamma <- rep(gamma, d)
   }
   if(is.null(P_inv)){
      P_inv <- pd.solve(P, log.det = TRUE)
      ldet <- attr(P_inv, "log.det")
   } else stopifnot(is.numeric(ldet))
   ## Compute pseudo-sample
   if(!is.null(pseudoskewt)){
      ## Most basic check 
      if(!is.matrix(pseudoskewt)) pseudoskewt <- cbind(pseudoskewt)
   } else {
      pseudoskewt <- if(all(gamma == gamma[1])){
         ## Equiskewness => one call to qskewt1d()
         matrix(qskewt1d_yos(as.vector(u), df = df, gamma = gamma[1]), 
                ncol = d)
      } else {
         sapply(1:d, function(i) qskewt1d_yos(u[, i], gamma = gamma[i], 
                                              df = df))
      }
   }
   pseudoskewt_mean_colvec <- matrix(colMeans(pseudoskewt), ncol = 1)
   ## Update 'scale'
   P_converged <- FALSE
   P_iter <- 1
   if(any(gamma != 0)){
      while(!P_converged & P_iter < P_maxiter){
         mahaXplusdf <- mahalanobis(pseudoskewt, center = FALSE, cov = P_inv,
                                    inverted = TRUE) + df
         mahagam <- mahalanobis(gamma, center = FALSE, cov = P_inv,
                                inverted = TRUE)
         besselarg <- sqrt(mahaXplusdf * mahagam)
         bessel1 <- besselK(besselarg, nu = (-(d + df)/2 + 1))
         bessel2 <- pmax(besselK(besselarg, nu = (-(d + df)/2)), .Machine$double.xmin)
         bessel3 <- besselK(besselarg, nu = (-(d + df)/2 - 1))
         delta <- 1 / sqrt(mahaXplusdf/mahagam) * bessel3 / bessel2
         eta_bar <- sum(sqrt(mahaXplusdf/mahagam) * bessel1 / bessel2)/n
         P_next <- as.matrix(Matrix::nearPD(
            (crossprod(sqrt(delta)*pseudoskewt)/n + eta_bar * tcrossprod(gamma)
             - tcrossprod(pseudoskewt_mean_colvec, gamma) -
                t(tcrossprod(pseudoskewt_mean_colvec, gamma))))$mat)
         P_next_inv <- pd.solve(P_next, log.det = TRUE)
         ldet <- attr(P_next_inv, "log.det")
         reldiff_P <- max(abs(P_next - P)/P)
         P_converged <- (reldiff_P < P_tol)
         P_iter <- P_iter + 1
         P <- P_next
         P_inv <- P_next_inv
      }
   } else {
      ## Standard multivariate t setting
      while(!P_converged & P_iter < P_maxiter){
         maha2_current <- mahalanobis(pseudoskewt, center = FALSE, 
                                      cov = P_inv, inverted = TRUE)
         weights <- (df + d) / (df + maha2_current)
         P_next <- as.matrix(nearPD(
            crossprod(sqrt(weights)*pseudoskewt)/n)$mat) 
         P_next_inv <- pd.solve(P_next, log.det = TRUE)
         ldet <- attr(P_next_inv, "log.det")
         reldiff_P <- max(abs(P_next - P)/P)
         P_converged <- (reldiff_P < P_tol)
         P_iter <- P_iter + 1
         P <- P_next
         P_inv <- P_next_inv
      }
   }
   P_next <- cov2cor(P_next)
   P_next_inv <- pd.solve(P_next, log.det = TRUE)
   ldet_next <- attr(P_next_inv, "log.det")
   ## Compute log-density
   ll <- if(report.ll){
      sum(dskewtcopula(u, scale = P_next, gamma = gamma,
                       df = df, log = TRUE))
   }
   ## Return
   list(P_next = P_next, P_next_inv = P_next_inv, ldet_next = ldet_next, 
        ll = ll)
}


#' Fitting  skew-t copulas
#' @param u (n, d) matrix of copula observations in (0,1)
#' @param fit.method string indicating the fitting method to be used 
#' @param df.init NULL or vector with initial estimates for 'df'; can contain NAs
#' @param gamma.init NULL or vector with initial estimates for 'gamma'
#' @param df.bounds 2-vector giving bounds on the dof parameter
#' @param gamma.bounds 2-vector giving bounds on the gamma parameter
#' @param control see ?get_set_param()
#' @param verbose logical if warnings shall be returned
#' @return S3 object of class 'fitgStudentcopula'
#' @author Erik Hintz
fitskewtcopula <- function(u, fit.method = c("EM-MLE", "Full-MLE"),
                           equiskewness = TRUE, 
                           df.init = NULL, gamma.init = NULL,
                           df.bounds = c(4, 15), gamma.bounds = c(-2, 2),
                           control = list(), verbose = TRUE,
                           optim.method = c("L-BFGS-B", "Nelder-Mead"),
                           control.optim = list(maxit = 250)){
   
   if(!equiskewness)
      stop("Currently, fitting is only available for equiskewness.")
   call <- match.call() # for return
   fit.method <- match.arg(fit.method) 
   optim.method <- match.arg(optim.method)
   ## Checks
   if(!is.matrix(u)) u <- cbind(u)
   d <- ncol(u)
   n <- nrow(u)
   ## Initial parameters
   df_current <- if(!is.null(df.init)) df.init else mean(df.bounds)
   gamma_current <- if(!is.null(gamma.init)) gamma.init else mean(gamma.bounds)
   P_current <- sin(pcaPP::cor.fk(u) * pi/2)
   P_current <- as.matrix(Matrix::nearPD(P_current, corr = TRUE)$mat)
   if(fit.method == "EM-MLE"){
      ECME.tol = list(ll = 1e-5, par = 1e-4,
                      P = 1e-4)
      ## Specify target-function passed to optim():
      ## @param 'par': 2-vector with 'df' and 'gamma'
      ## @return neg.max.ll numeric (>0) which is 
      ## - argmax_{\df, gamma} log l(df, \hat{P}(df))
      ## where \hat{P}(df) is the MLE for P for fixed 'df' and l is the
      ## copula likelihood function 
      P_myfun <- P_current 
      calls <- 0
      negloglik <- function(par){
         if(par[1] < df.bounds[1] | par[1] > df.bounds[2] | 
            par[2] < gamma.bounds[1] | par[2] > gamma.bounds[2]) return(1e99)
         temp <- fitscaleskewtEM(u, pseudoskewt = NULL, df = par[1], 
                                 gamma = rep(par[2], d), 
                                 P = P_myfun, P_inv = NULL, ldet = NULL, 
                                 report.ll = TRUE, 
                                 P_maxiter = 100, P_tol = ECME.tol$par)
         ## Return -max ll achieved 
         P_myfun <<- temp$P_next
         -temp$ll
      }
      opt.obj <- if(optim.method == "L-BFGS-B") optim(c(df_current, gamma_current), negloglik, method = "L-BFGS-B", 
               lower = c(df.bounds[1], gamma.bounds[1]), upper = c(
                  df.bounds[2], gamma.bounds[2]), control = control.optim) else optim(c(df_current, gamma_current), negloglik, 
                           control = control.optim)                                                                
      df_next <- opt.obj$par[1]
      gamma_next <- opt.obj$par[2]
      ll_next <- -opt.obj$value
      P_next <- P_myfun
      list(df = df_next, scale = P_next, ll = ll_next, gamma = gamma_next)
   } else { # full MLE
      ## Compute *internal* starting values
      initial.int.pars <- org_to_int_pars_st(rho = P_current[lower.tri(P_current)], 
                                             df = df_current, gamma = gamma_current)
      lth <- length(initial.int.pars) - 2 # last elements are 'df' and 'gamma'
      ## Construct bounds
      upbounds <- matrix(pi, ncol = d - 1, nrow = d - 1)
      diag(upbounds) <- 2 * pi
      upbounds <- as.vector(upbounds[lower.tri(upbounds, diag = TRUE)])
      up_ <- c(upbounds, log(df.bounds[2]-2), gamma.bounds[2])
      lo_ <- c(rep(0, lth), log(df.bounds[1]-2), gamma.bounds[1])
      ## Set up log-likelihood as a function of *internal* parameters 
      negloglik_ <- function(intpars){
         if(any(intpars < lo_) | any(intpars > up_)) return(Inf)
         orgpars <- int_to_org_pars_st(intpars) 
         scale <- rhoToOmega(orgpars$rho)
         factor <- tryCatch(chol(scale), error = function(e) NULL)
         if(is.null(factor)) return(+1e18) # scale does not have full rank
         scale.inv <- chol2inv(factor, size = d)
         ldet <- sum(log(diag(factor)))
         -sum(dskewtcopula(u, df = orgpars$df, gamma = rep(orgpars$gamma, d), 
                           scale = scale, scale.inv = scale.inv, ldet = ldet, 
                           log = TRUE))
      }
  
      ## Call the optimizer 
      opt.obj <- if(optim.method == "L-BFGS-B") 
         optim(initial.int.pars, negloglik_, method = "L-BFGS-B", lower = c(
            rep(0, lth), log(df.bounds[1] - 2), gamma.bounds[1]), upper = c(
               upbounds, log(df.bounds[2] - 2), gamma.bounds[2])) else optim(
                  initial.int.pars, negloglik_)
      newpars <- int_to_org_pars_st(opt.obj$par) 
      list(df = newpars$df, scale = rhoToOmega(newpars$rho), 
           gamma = newpars$gamma, ll = -opt.obj$value)
   }
}


