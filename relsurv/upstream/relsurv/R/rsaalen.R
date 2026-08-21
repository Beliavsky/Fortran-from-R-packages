# Calculating betas:
calculateBetasRelsurv22R <- function(data, xt, event_times, all_times,
                                     ratetable, atts, data_mat, variance, var_estimator='dN'){

  ncol <- length(all_times)

  Yt_val <- Yt(data, all_times)
  dNt_val <- dNt(data, all_times)
  dLambdaP_val <- dLambdaPR(data_mat, all_times, event_times, ratetable, atts)
  # dLambdaP_val <- c(0, dLambdaP_val)

  sample_size = nrow(xt)
  number_covs = ncol(xt)

  diag_dNt = matrix(0, nrow=sample_size, ncol=sample_size)
  beta_var = matrix(0, nrow=number_covs, ncol=number_covs)

  betas_var_list <- vector("list", ncol)

  if(variance){
    betas_list0 <- lapply(1:ncol, function(i) fitOLS2(prepareX(Yt_val[[i]], xt), dNt_val[[i]] - (dLambdaP_val[,i+1] - dLambdaP_val[,i]), Yt_val[[i]]))

    if(var_estimator=='dN'){
      betas_var_list <- lapply(1:ncol, function(i) betas_list0[[i]][[2]] %*% diag(dNt_val[[i]] - (dLambdaP_val[,i+1] - dLambdaP_val[,i])) %*% t(betas_list0[[i]][[2]]))
    } else if(var_estimator=='XdB'){
      betas_var_list <- lapply(1:ncol, function(i) betas_list0[[i]][[2]] %*% diag(as.vector(prepareX(Yt_val[[i]], xt) %*% betas_list0[[i]][[1]])) %*% t(betas_list0[[i]][[2]]))
    }
    betas_list <- lapply(1:ncol, function(i) betas_list0[[i]][[1]])
  } else{
    betas_list <- lapply(1:ncol, function(i) fitOLS2(prepareX(Yt_val[[i]], xt), dNt_val[[i]] - (dLambdaP_val[,i+1] - dLambdaP_val[,i]), Yt_val[[i]])[[1]])
  }
  # for(i in 1:ncol){
  #   xx2 = prepareX(Yt_val[[i]], xt)
  #   betas = fitOLS2(xx2, dNt_val[[i]], Yt_val[[i]])
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

# Calculating betas and gamma:
calculateBetasGammasRelsurv22R <- function(data, xt, zt, event_times, all_times,
                                     ratetable, atts, data_mat, var_estimator=1){

  ncol <- length(all_times)

  Yt_val <- Yt(data, all_times)
  dNt_val <- dNt(data, all_times)
  dLambdaP_val <- dLambdaPR(data_mat, all_times, event_times, ratetable, atts)
  # dLambdaP_val <- c(0, dLambdaP_val)

  sample_size = nrow(xt)
  number_covs = ncol(xt)

  diag_dNt = matrix(0, nrow=sample_size, ncol=sample_size)
  beta_var = matrix(0, nrow=number_covs, ncol=number_covs)

  diff_t <- diff(c(min(data$start), all_times))

  betas_var_list <- vector("list", ncol)

  mod_list <- lapply(1:ncol, function(i) fitOLSconst(prepareX(Yt_val[[i]], xt),
                                                                 prepareX(Yt_val[[i]], zt)[,2:(ncol(zt)+1), drop=FALSE],
                                                                   dNt_val[[i]] - (dLambdaP_val[,i+1] - dLambdaP_val[,i]),
                                                                 Yt_val[[i]])
  )

  # Prvi integral:
  int_ztHz <- Reduce('+', lapply(1:ncol, function(i) mod_list[[i]][[1]]*diff_t[i]))
  int_ztHz_inv <- solve(int_ztHz)

  # Drugi integral:
  int_ztHdN <- Reduce('+', lapply(1:ncol, function(i) mod_list[[i]][[2]]))

  gamma_coef <- int_ztHz_inv%*%int_ztHdN

  Xminus <- lapply(1:ncol, function(i) mod_list[[i]][[3]])

  betas_list <- lapply(1:ncol, function(i) Xminus[[i]]%*%(matrix(dNt_val[[i]], ncol=1) - zt %*% gamma_coef * diff_t[i]))

  out <- list(betas_list=betas_list, betas_var_list=betas_var_list, gamma_coef=gamma_coef)

  return(out)
}

#' Fit an extended additive hazards model using relative survival
#'
#' Fits the Aalen additive hazard model using relative survival. The function can be used for multi-state model
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
#' @param var_estimator Choose variance estimator, in the same way as in survaalen. The default option 'dN' uses dN(t)-dLambda_P(t) in the variance estimator, equivalent to formula 4.63 in Aalen et al. (2008). Option 'XdB' uses X*dB(t), see formula 4.64 in Aalen et al. (2008).
#' @param ratetable a table of event rates, organized as a \code{ratetable} object, such as \code{slopop}.
#' @param rmap an optional list to be used if the variables are not organized and named in the same way as in the \code{ratetable} object.
#' @param split.transitions only relevant if a multi-state model is fitted. An integer vector containing the numbered transitions that should be split. Use same numbering as in the given transition matrix.
#' @param add.times specific times at which estimates should be obtained.
#' @return An object of class \code{aalen.model}.
#' @seealso \code{\link{survaalen}}
#' @author Damjan Manevski
#' @keywords survival
#' @examples
#'
#' # Survival:
#' data(rdata)
#' mod <- rsaalen(Surv(time, cens)~sex+age, data=rdata, ratetable=slopop,
#'      rmap=list(age=age*365.241))
#' head(mod$coefficients)
#' tail(mod$coefficients)
#'
#' # Multi-state model:
#' data(ebmt1wide)
#' # Generate sex and year data (for illustrative purposes since it is not given in the data):
#' ebmt1wide$sex <- sample(1:2, size = nrow(ebmt1wide), replace = TRUE)
#' ebmt1wide$year <- as.Date('2010-01-01')
#'
#' mod <- rsaalen(Surv(Tstart, Tstop, status)~age.1+age.2+age.3+strata(trans), data=ebmt1wide,
#'                  ratetable = slopop, rmap = list(age=age*365.241), split.transitions = 2:3)
#' head(mod$coefficients$trans1)
#' head(mod$coefficients$trans2)
#' head(mod$coefficients$trans3)
rsaalen <- function(formula, data, variance=FALSE, var_estimator='dN',
                    ratetable=relsurv::slopop, rmap, split.transitions,
                    add.times){

  # Find covariates and strata:
  covs <- attr(terms(formula), 'term.labels')
  strata_obj <- grep("strata\\(", covs, value=TRUE)

  # Prepare objects in case of const():
  formula_new <- formula
  covs_new <- covs
  constTRUE <- any(grepl('const', deparse(formula[[3]])))
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

    relsurv.mods <- vector("list", length(split.transitions))

    # Check:
    if(length(strata_obj)>1) stop('You have supplied multiple strata() objects in the formula. Please supply only one.')
    if(missing('split.transitions')) stop('Please define which transitions should be split in the split.transitions argument. If relative survival is not needed, use the survaalen function.')

    # Take care of rmap:
    if (!missing(rmap)) {
      rmap_tmp <- substitute(rmap)
      if(inherits(rmap_tmp, "call")){
        rmap <- rmap_tmp
      }
    }

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

    stevec <- 1

    # For every strata:
    for(i in 1:length(strata_levels)){
      # Subset data and covs:
      data_tmp <- data[data[,strata_obj]==strata_levels[i], ]
      covs_tmp <- covs[endsWith(covs, paste0(".", strata_levels[i])) | endsWith(covs, paste0(".", strata_levels[i], ')'))] # grep(paste0('.', strata_levels[i]), covs_new, value=TRUE)

      # Prepare formula:
      covs_tmp2 <- paste0(covs_tmp, collapse=' + ')
      formula_tmp <- as.formula(paste0(deparse(formula_new[[2]]), '~', covs_tmp2))

      # Run model:
      if(i %in% split.transitions){
        if(!missing(add.times)){
          mod_tmp <- rsaalen(formula_tmp, data_tmp, variance, var_estimator, ratetable, rmap, add.times=add.times)
        } else{
          mod_tmp <- rsaalen(formula_tmp, data_tmp, variance, var_estimator, ratetable, rmap)
        }

        relsurv.mods[[stevec]] <- mod_tmp
        stevec <- stevec+1
      } else{
        if(!missing(add.times)){
          mod_tmp <- survaalen(formula_tmp, data_tmp, variance, add.times)
        } else{
          mod_tmp <- survaalen(formula_tmp, data_tmp, variance)
        }
      }

      # Remove time 0:
      # if(mod_tmp$coefficients[1,1]==0){
      #   mod_tmp$coefficients <- mod_tmp$coefficients[2:nrow(mod_tmp$coefficients),]
      # }

      # Save coefficients:
      coefficients[[i]] <- mod_tmp$coefficients
      coefficients.var[[i]] <- mod_tmp$coefficients.var

      all_times <- c(all_times, mod_tmp$coefficients[,1])

      if(constTRUE){
        if('gamma' %in% names(mod_tmp)){
          gamma[[i]] <- mod_tmp$gamma
          # if(variance){
          #   gamma.var[[i]] <- mod_tmp$gamma.var
          # }
        }
      }
    }
    names(coefficients) <- paste0('trans', strata_levels)
    if(variance) names(coefficients.var) <- paste0('trans', strata_levels)
    names(gamma) <- paste0('trans', strata_levels)
    # names(gamma.var) <- paste0('trans', strata_levels)

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
          if((nrow(coefficients.var[[i]]) != 0)){
            coefficients.var[[i]][,j] <- mstateNAfix(coefficients.var[[i]][,j], 0)
          }
        }
      }
    }

    names(relsurv.mods) <- split.transitions

    if(constTRUE){
      if(variance){
        out <- list(coefficients=coefficients, coefficients.var=coefficients.var,
                    gamma=gamma, #gamma.var=gamma.var,
                    split.transitions=split.transitions,
                    formula=formula, ratetable=ratetable, rmap=rmap,
                    relsurv.mods=relsurv.mods)
      } else{
        out <- list(coefficients=coefficients,
                    gamma=gamma,
                    split.transitions=split.transitions,
                    formula=formula, ratetable=ratetable, rmap=rmap,
                    relsurv.mods=relsurv.mods)
      }
    } else{
      if(variance){
        out <- list(coefficients=coefficients, coefficients.var=coefficients.var,
                    split.transitions=split.transitions,
                    formula=formula, ratetable=ratetable, rmap=rmap,
                    relsurv.mods=relsurv.mods)
      } else{
        out <- list(coefficients=coefficients,
                    split.transitions=split.transitions,
                    formula=formula, ratetable=ratetable, rmap=rmap,
                    relsurv.mods=relsurv.mods)
      }
    }

    class(out) <- 'aalen.model'

    return(out)

  } else{
    if(!missing('split.transitions')) warning('The split.transitions argument is ignored since you are not using a multi-state model.')

    # Save original data:
    data_orig <- data

    # Take care of rmap:
    if (!missing(rmap)) {
      rmap_tmp <- substitute(rmap)
      if(inherits(rmap_tmp, "call")){
        rmap <- rmap_tmp
      }
    }

    # Standard relsurv format:
    rform <- rformulate(formula_new, data, ratetable, stats::na.omit(),
                        rmap)

    # Objects:
    data <- rform$data
    xt <- rform$X
    ratetable <- rform$ratetable

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
    # dNt <- dNt(data, all_times)
    # Yt <- Yt(data, all_times)
    # dLambdaPs_0 <- dLambdaP0(data, all_times, event_times, ratetable, atts)

    atts <- attributes(ratetable)
    atts$cutpoints <- atts$cutpoints[lapply(atts$cutpoints,length)>0]
    # dLambdaPs <- dLambdaP(as.matrix(data), all_times, event_times, as.vector(ratetable), atts)

    # xx1 <- prepareX(Yt[[1]], as.matrix(xt))
    # fuu <- fitOLS2(xx1, dNt[[1]], Yt[[1]])

    # Find betas:
    # betas <- calculateBetas(data, as.matrix(xt), event_times)
    # betas <- calculateBetasRelsurv(data, as.matrix(xt), event_times, all_times, as.vector(ratetable), atts, as.matrix(data))
    # betas <- calculateBetasRelsurv2(data, as.matrix(xt), event_times, all_times, as.vector(ratetable), atts, as.matrix(data))

    # betas <- calculateBetasRelsurv22(data, as.matrix(xt), event_times, all_times, as.vector(ratetable), atts, as.matrix(data), 1)
    if(constTRUE){
      zt <- xt[,covs_const, drop=FALSE]
      xt <- xt[,!covs_const, drop=FALSE]

      # all_times_w0 <- all_times[2:length(all_times)]
      # betas0 <- calculateBetasGammasRelsurv22R(data, as.matrix(xt), as.matrix(zt), event_times, all_times_w0, ratetable, atts, as.matrix(data), 1)
      betas0 <- calculateBetasGammasRelsurv22R(data, as.matrix(xt), as.matrix(zt), event_times, all_times, ratetable, atts, as.matrix(data), 1)
    } else{
      betas0 <- calculateBetasRelsurv22R(data, as.matrix(xt), event_times, all_times, ratetable, atts, as.matrix(data), variance, var_estimator)
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
        betas_var2[,iu] <- cumsum(betas_var2[,iu])
      }
    }

    # Add times:
    if(constTRUE){
      # betas2 <- cbind(all_times_w0, betas2)
      # betas2 <- betas2[all_times_w0 %in% c(0, event_times),]
      betas2 <- cbind(all_times, betas2)
    } else{
      betas2 <- cbind(all_times, betas2)
    }

    # Add times:
    if(ncol(betas_var2)!=0){
      betas_var2 <- cbind(all_times, betas_var2)
    }

    # Add naming:
    colnames(betas2) <- c('time', '(Intercept)', colnames(xt))
    if(ncol(betas_var2)!=0){
      colnames(betas_var2) <- c('time', '(Intercept)', colnames(xt))
    }

    # Don't have double zeros:
    # if(betas2[1,1] == 0 & betas2[2,1] == 0){
    #   betas2 <- betas2[-1,]
    # }

    betas3 <- betas2[betas2[,1] %in% c(0, event_times), ]

    if(!(0 %in% betas3[,1])){
      betas3 <- rbind(matrix(0, nrow=1,ncol=ncol(betas3)), betas3)
    }

    if(ncol(betas_var2)!=0){
      betas_var2 <- betas_var2[betas_var2[,1] %in% c(0, event_times), ]
    }

    # Prepare final object:
    if(variance){
      var.obj <- betas_var2
      # var.obj[,2:ncol(var.obj)] <- abs(var.obj[,2:ncol(var.obj)])/4

      out <- list(coefficients=betas3, coefficients.var=var.obj, # prej je betas_var2
                  formula=formula, ratetable=ratetable, rmap=rmap)
    } else{
      out <- list(coefficients=betas3,
                  formula=formula, ratetable=ratetable, rmap=rmap)
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
# implement variances
# other left-truncation solutions? for now you put 0s
# time dependent covariates?

# multi-state:
# Popravi racunanje vo multi-state - za sea site vreminja gi vlecis, ne zemas vo predvid promeni vo L_P
