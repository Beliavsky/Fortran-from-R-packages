# Calculating betas:
calculateBetasR <- function(data, xt, event_times, variance=FALSE, var_estimator='dN'){

  ncol <- length(event_times)

  Yt_val <- Yt(data, event_times)
  dNt_val <- dNt(data, event_times)

  sample_size = nrow(xt)
  number_covs = ncol(xt)

  # diag_dNt = matrix(0, nrow=sample_size, ncol=sample_size)
  # beta_var = matrix(0, nrow=number_covs, ncol=number_covs)

  betas_var_list <- vector("list", ncol)

  if(variance){
    betas_list0 <- lapply(1:ncol, function(i) fitOLS2(prepareX(Yt_val[[i]], xt), dNt_val[[i]], Yt_val[[i]]))

    if(var_estimator=='dN'){
      betas_var_list <- lapply(1:ncol, function(i) betas_list0[[i]][[2]] %*% diag(dNt_val[[i]]) %*% t(betas_list0[[i]][[2]]))
    } else if(var_estimator=='XdB'){
      betas_var_list <- lapply(1:ncol, function(i) betas_list0[[i]][[2]] %*% diag(as.vector(prepareX(Yt_val[[i]], xt) %*% betas_list0[[i]][[1]])) %*% t(betas_list0[[i]][[2]]))
    }
    betas_list <- lapply(1:ncol, function(i) betas_list0[[i]][[1]])
  } else{
    betas_list <- lapply(1:ncol, function(i) fitOLS2(prepareX(Yt_val[[i]], xt), dNt_val[[i]], Yt_val[[i]])[[1]]
    )
  }


  # for(i in 1:ncol){
  #   xx2 = relsurv:::prepareX(Yt_val[[i]], xt)
  #   betas = relsurv:::fitOLS2(xx2, dNt_val[[i]], Yt_val[[i]])
  #
  #   betas_list[[i]] = betas[[1]]
  #
  #   if(var_estimator==1){
  #     diag(diag_dNt) <- dNt_val[[i]]
  #   } else if(var_estimator == 2){
  #     betas0_vec <- betas[[1]];
  #     diag(diag_dNt) <- xx2 %*% betas0 # mogoce je treba dati na betas0 t()
  #   }
  #
  #   betas1 = betas[[2]]
  #   # betas_var_list[[i]] <- betas1 %*% diag_dNt %*% t(betas1)
  # }

  out <- list(betas_list, betas_var_list)

  return(out)
}

# An R version of the fitOLSconst C function:
fitOLSconstR <- function(mX, mZ, dNt, Yt) {

  no_cov = ncol(mX)
  no_cov_Z = ncol(mZ)
  sample_size = nrow(mX)
  no_at_risk = sum(Yt)

  Xminus = matrix(0, nrow = sample_size, ncol = no_cov)
  H = matrix(0, ncol=no_cov, nrow=no_cov)
  Identity = matrix(0, ncol=sample_size, nrow=sample_size)
  diag(Identity) <- 1

  prvaKomponenta = matrix(0, ncol=no_cov_Z, nrow=no_cov_Z)
  drugaKomponenta = matrix(0, ncol=no_cov_Z, nrow=1)

  out <- list()
  out[[1]] = prvaKomponenta
  out[[2]] = drugaKomponenta
  out[[3]] = Xminus

  if(no_at_risk >= no_cov){
    mXtX = t(mX)%*%mX

    rcf = rcond(mXtX)

    if(rcf != 0){
      Xminus = solve(mXtX)%*%t(mX)

      H = Identity-mX%*%Xminus
      prvaKomponenta = t(mZ)%*%H%*%mZ
      drugaKomponenta = t(mZ)%*%H%*%dNt

      out[[1]] = prvaKomponenta
      out[[2]] = drugaKomponenta
      out[[3]] = Xminus
    }
  }
  return(out)
}

# Calculating betas and gamma:
calculateBetasGammasR <- function(data, xt, zt, event_times, var_estimator='dN'){

  ncol <- length(event_times)

  Yt_val <- Yt(data, event_times)
  dNt_val <- dNt(data, event_times)

  sample_size = nrow(xt)
  number_covs = ncol(xt)

  diff_t <- diff(c(min(data$start), event_times))

  # diag_dNt = matrix(0, nrow=sample_size, ncol=sample_size)
  # beta_var = matrix(0, nrow=number_covs, ncol=number_covs)

  betas_var_list <- vector("list", ncol)

  mod_list <- lapply(1:ncol, function(i) fitOLSconst(prepareX(Yt_val[[i]], xt),
                                                               prepareX(Yt_val[[i]], zt)[,2:(ncol(zt)+1), drop=FALSE],
                                                               dNt_val[[i]],
                                                               Yt_val[[i]]))

  # Prvi integral:
  int_ztHz <- Reduce('+', lapply(1:ncol, function(i) mod_list[[i]][[1]]*diff_t[i]))
  int_ztHz_inv <- solve(int_ztHz)

  # Drugi integral:
  int_ztHdN <- Reduce('+', lapply(1:ncol, function(i) mod_list[[i]][[2]]))

  gamma_coef <- int_ztHz_inv%*%int_ztHdN

  Xminus <- lapply(1:ncol, function(i) mod_list[[i]][[3]])

  betas_list <- lapply(1:ncol, function(i) Xminus[[i]]%*%(matrix(dNt_val[[i]], ncol=1) - zt %*% gamma_coef * diff_t[i]))

  out <- list(betas_list=betas_list, betas_var_list=betas_var_list, gamma_coef=gamma_coef)

  # Reduce('+', lapply(1:49, function(i) solve(mod_list[[i]][[1]]) %*% mod_list[[i]][[2]]))/max(data$Y)
  # sum(lapply(1:ncol, function(i) solve(mod_list[[i]][[1]]) %*% mod_list[[i]][[2]]))/max(data$Y)

  return(out)
}

#' Fit an additive hazards model
#'
#' Fits the Aalen additive hazard model. The function can be used for multi-state model
#' data (as in the package mstate; class msdata) by supplying the start and stop times in the
#' Surv object and adding a strata(trans) object in the formula (where trans denotes the
#' transition in the multi-state model).
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param variance a logical value indicating whether the variances of the hazards should be computed.
#' Default is FALSE.
#' @param var_estimator Choose variance estimator. The default option 'dN' uses dN(t) in the variance formula, see formula 4.63 in Aalen et al. (2008). Option 'XdB' uses X*dB(t), see formula 4.64 in Aalen et al. (2008).
#' @param add.times specific times at which estimates should be obtained.
#' @return An object of class \code{aalen.model}.
#' @seealso \code{\link{rsaalen}}
#' @author Damjan Manevski
#' @keywords survival
#' @examples
#'
#' # Survival:
#' data(rdata)
#' mod <- survaalen(Surv(time, cens)~sex+age, data=rdata)
#' head(mod$coefficients)
#' tail(mod$coefficients)
#'
#' # Multi-state model:
#' data(ebmt1wide)
#' mod <- survaalen(Surv(Tstart, Tstop, status)~age.1+age.2+age.3+strata(trans), data=ebmt1wide)
#' head(mod$coefficients$trans1)
#' head(mod$coefficients$trans2)
#' head(mod$coefficients$trans3)
survaalen <- function(formula, data, variance=FALSE, var_estimator='dN',
                      add.times){

  # Find covariates and strata:
  covs <- attr(terms(formula), 'term.labels')
  strata_obj <- grep("strata\\(", covs, value=TRUE)

  # Prepare objects in case of const():
  formula_new <- formula
  covs_new <- covs
  constTRUE <- grepl('const', deparse(formula[[3]]))
  if(constTRUE){
    covs_const <- grepl('const\\(', covs)
    rnames_gamma <- covs[covs_const]

    covs_wconst <- covs
    covs_wconst[covs_const] <- gsub('const\\(', '', covs_wconst[covs_const])
    covs_wconst[covs_const] <- gsub('\\)', '', covs_wconst[covs_const])

    covs_new <- covs_wconst

    if(length(covs_wconst)>1){
      covs_wconst <- paste0(covs_wconst, collapse = ' + ')
    }
    formula_new <- as.formula(paste0(deparse(formula[[2]]), '~', covs_wconst))
  }

  # Run a multi-state model:
  if(length(strata_obj)>0){

    # Check:
    if(length(strata_obj)>1) stop('You have supplied multiple strata() objects in the formula. Please supply only one.')

    # Save original data:
    data_orig <- data
    data <- as.data.frame(data)

    # Find strata object:
    strata_obj <- gsub("strata\\(|\\)", "", strata_obj)
    strata_levels <- unique(data[,strata_obj])

    coefficients <- list()
    coefficients.var <- list()
    gamma <- vector("list", length = length(strata_levels))
    # gamma.var <- vector("list", length = length(strata_levels))

    all_times <- c()

    # Find max time:
    max_time <- max(data[,as.character(formula_new[[2]][3])])

    # For every strata:
    for(i in 1:length(strata_levels)){
      # Subset data and covs:
      data_tmp <- data[data[,strata_obj]==strata_levels[i], ]
      covs_tmp <- covs[endsWith(covs, paste0(".", strata_levels[i])) | endsWith(covs, paste0(".", strata_levels[i], ')'))] # grep(paste0('.', strata_levels[i]), covs_new, value=TRUE)

      # Prepare formula:
      covs_tmp2 <- paste0(covs_tmp, collapse=' + ')
      formula_tmp <- as.formula(paste0(deparse(formula_new[[2]]), '~', covs_tmp2))

      # Run model:
      if(!missing(add.times)){
        mod_tmp <- survaalen(formula_tmp, data_tmp, variance, add.times=add.times)
      } else{
        mod_tmp <- survaalen(formula_tmp, data_tmp, variance)
      }

      # Remove time 0:
      # if(mod_tmp$coefficients[1,1]==0){
      #   mod_tmp$coefficients <- mod_tmp$coefficients[2:nrow(mod_tmp$coefficients),]
      # }

      # Save coefficients:
      coefficients[[i]] <- mod_tmp$coefficients
      if(variance) coefficients.var[[i]] <- mod_tmp$coefficients.var

      all_times <- c(all_times, mod_tmp$coefficients[,1])

      if(constTRUE){
        if('gamma' %in% names(mod_tmp)){
          gamma[[i]] <- mod_tmp$gamma
          # gamma.var[[i]] <- mod_tmp$gamma.var
        }
      }
    }
    names(coefficients) <- paste0('trans', strata_levels)
    names(gamma) <- paste0('trans', strata_levels)

    if(variance){
      names(coefficients.var) <- paste0('trans', strata_levels)
      # names(gamma.var) <- paste0('trans', strata_levels)
    }

    # Add all times in the mstate model:
    all_times <- sort(unique(c(all_times, max_time)))
    for(i in 1:length(strata_levels)){
      all_times_tmp <- all_times[!(all_times %in% coefficients[[i]][,1])]
      tmp_df <- matrix(NA, nrow=length(all_times_tmp), ncol= ncol(coefficients[[i]]))

      tmp_df[,1] <- all_times_tmp
      colnames(tmp_df) <- colnames(coefficients[[i]])

      coefficients[[i]] <- rbind(coefficients[[i]], tmp_df)
      coefficients[[i]] <- coefficients[[i]][order(coefficients[[i]][,1]),]

      if(variance){
        if(nrow(coefficients.var[[i]]) != 0){
          coefficients.var[[i]] <- rbind(coefficients.var[[i]], tmp_df)
          coefficients.var[[i]] <- coefficients.var[[i]][order(coefficients.var[[i]][,1]),]
        }
      }

      for(j in 2:ncol(coefficients[[i]])){
        coefficients[[i]][,j] <- mstateNAfix(coefficients[[i]][,j], 0)
        if(variance){
          if(nrow(coefficients.var[[i]]) != 0){
            coefficients.var[[i]][,j] <- mstateNAfix(coefficients.var[[i]][,j], 0)
          }
        }
      }
    }

    if(constTRUE){
      if(variance){
        out <- list(coefficients=coefficients, coefficients.var=coefficients.var,
                    gamma=gamma#, gamma.var=gamma.var
                    )
      } else{
        out <- list(coefficients=coefficients, gamma=gamma)
      }
    } else{
      if(variance){
        out <- list(coefficients=coefficients, coefficients.var=coefficients.var)
      } else{
        out <- list(coefficients=coefficients)
      }
    }
    class(out) <- 'aalen.model'

    return(out)

    # Run a usual survival model:
  } else{
    # Save original data:
    data_orig <- data
    sample_size <- nrow(data_orig)

    # Standard relsurv format:
    rform <- rformulate2(formula_new, data)

    data <- rform$data
    xt <- rform$X

    # Times:
    all_times <- unique(c(data$start, data$Y))
    event_times <- unique(data$Y[data$stat==1])

    all_times <- sort(all_times)
    event_times <- sort(event_times)

    # Take only times until last event time:
    all_times <- all_times[all_times <= event_times[length(event_times)]]

    if(!missing(add.times)){
      add.times <- add.times[add.times<=all_times[length(all_times)]]
      all_times <- sort(unique(c(add.times, all_times)))
      event_times <- sort(unique(c(add.times, event_times)))
    }

    # Find Yt / dNt:
    # dNt <- dNt(data, event_times)
    # Yt <- Yt(data, event_times)

    # xx1 <- prepareX(Yt[[1]], as.matrix(xt))
    # fuu <- fitOLS2(xx1, dNt[[1]], Yt[[1]])

    # browser()

    # Find betas:
    # betas <- calculateBetas(data, as.matrix(xt), event_times, 1)
    if(constTRUE){
      zt <- xt[,covs_const, drop=FALSE]
      xt <- xt[,!covs_const, drop=FALSE]

      all_times_w0 <- all_times[2:length(all_times)]

      betas0 <- calculateBetasGammasR(data, as.matrix(xt), as.matrix(zt), all_times_w0, 1)
    } else{
      betas0 <- calculateBetasR(data, as.matrix(xt), event_times, variance, var_estimator)
    }

    # Variance:
    betas_var <- betas0[[2]]
    betas_var <- lapply(betas_var, function(x) diag(x))
    betas_var2 <- do.call(rbind, betas_var)

    # Point estimates:
    betas <- betas0[[1]]
    # Take care of format:
    betas2 <- t(do.call(cbind, betas))
    # betas2 <- rbind(matrix(0, nrow=1, ncol=(ncol(xt)+1)),
    #                 betas2)

    # Cumulative:
    for(iu in 1:ncol(betas2)){
      betas2[,iu] <- cumsum(betas2[,iu])
      if(ncol(betas_var2)!=0){
        betas_var2[,iu] <- cumsum(betas_var2[,iu])#*sample_size
      }
    }

    # Add zero:
    if(0 %in% all_times){
      event_times <- c(0, event_times)
      betas2 <- rbind(matrix(0, nrow=1, ncol=(ncol(xt)+1)),
                      betas2)
      if(ncol(betas_var2)!=0){
        betas_var2 <- rbind(matrix(0, nrow=1, ncol=(ncol(xt)+1)),
                            betas_var2)
      }
    }

    # Add times:
    if(constTRUE){
      betas2 <- cbind(all_times, betas2)
      betas2 <- betas2[all_times %in% c(0, event_times),]
    } else{
      betas2 <- cbind(event_times, betas2)
    }
    if(ncol(betas_var2)!=0){
      betas_var2 <- cbind(event_times, betas_var2)
    }

    # Add naming:
    colnames(betas2) <- c('time', '(Intercept)', colnames(xt))
    if(ncol(betas_var2)!=0){
      colnames(betas_var2) <- c('time', '(Intercept)', colnames(xt))
    }

    # Prepare final object:
    if(variance){
      var.obj <- betas_var2
      # var.obj[,2:ncol(var.obj)] <- abs(var.obj[,2:ncol(var.obj)])/4

      out <- list(coefficients=betas2, coefficients.var=var.obj) # prej je betas_var2
    } else{
      out <- list(coefficients=betas2)
    }

    if(constTRUE){
      rownames(betas0[[3]]) <- rnames_gamma
      colnames(betas0[[3]]) <- 'gamma'

      out[[length(out)+1]] <- betas0[[3]]
      names(out)[length(out)] <- 'gamma'

      # if(variance){
      #   out[[length(out)+1]] <- betas0[[3]] # tmp
      #   colnames(out[[length(out)]]) <- 'gamma.var'
      #   names(out)[length(out)] <- 'gamma.var'
      # }
    }

    out$formula <- formula
    out$terms <- terms(formula)

    class(out) <- 'aalen.model'

    return(out)
  }
}

# TO DO:
# Scheike vprasanje: zakaj premikas case, ce je ob istem casu vec dogodkov?
# Some plot function for mstate output of surv/rsaalen? Za coefficiente?
