colVars <- function(x, na.rm = FALSE){
  f <- function(v, na.rm = na.rm) {
    if(is.numeric(v) || is.logical(v) || is.complex(v))
      stats::var(v, na.rm = na.rm)
    else NA
  }
  return(unlist(lapply(x, f, na.rm = na.rm)))
}

# Copied function from mstate:::NAfix.
mstateNAfix <- function (x, subst = -Inf)
{
  spec <- max(x[!is.na(x)]) + 1
  x <- c(spec, x)
  while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)] -
                                           1]
  x[x == spec] <- subst
  x <- x[-1]
  x
}

# Helper function:
nessie_spi <- function(formula = formula(data), data, ratetable = relsurv::slopop,
                       tis, starting.time, include.censoring=FALSE,
                       arg.example=FALSE, rmap){
  data_orig <- data
  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  na.action <- NA
  rform <- rformulate(formula, data, ratetable, na.action,
                      rmap)

  data <- rform$data
  data$Xs <- rep(1, nrow(data))
  n_rows <- nrow(data)

  # Fix demographic covariates:
  if(starting.time == "left.truncated"){
    rform$R[,"year"] <- rform$R[,"year"] - rform$R[,"age"]
    rform$R[,"age"] <- 0
  }
  if(include.censoring){
    # browser()
    wh <- which(rform$status==1)
    rform$Y[wh] <- max(rform$Y)

    if(arg.example){
      wh2 <- which(rform$status==1 & data$age==18262)
      rform$Y[wh2] <- 1826
    }
  }
  else{
    rform$Y <- rep(max(rform$Y), length(rform$Y))
    # status is not relevant in this case
  }

  out <- NULL
  out$yi <- NULL
  out$yidli <- NULL
  l_tis <- length(tis)
  temps <- lapply(1:n_rows, function(inx) {
    temp <- exp_prep(rform$R[inx, , drop = FALSE], rform$Y[inx], rform$ratetable,
                     rform$status[inx], times = tis, fast = TRUE, cmp=FALSE,ys=data$start[inx])
    s_pi <- exp(-cumsum(temp$yidli))

    s_pi_helper <- which.min(temp$yidli==0)-1
    if(s_pi_helper>1){ s_pi[1:s_pi_helper] <- 0}

    if(include.censoring){ s_pi[(s_pi_helper+1):l_tis] <- pmin(s_pi[(s_pi_helper+1):l_tis],
                                                               temp$yi[(s_pi_helper+1):l_tis])}

    c(s_pi, # s_pi
      temp$yidli*s_pi) # l_pi * s_pi
  })

  temps2 <- do.call("cbind", temps)

  temps2 <- rowSums(temps2)

  out$yi <- temps2[1:(length(temps2)/2)]
  out$yidli <- temps2[(length(temps2)/2+1):length(temps2)]
  return(out)
}

# Copied scales::trans_new:
# scales_trans_new <- function (name, transform, inverse, breaks = extended_breaks(),
#           minor_breaks = regular_minor_breaks(), format = format_format(),
#           domain = c(-Inf, Inf))
# {
#   if (is.character(transform))
#     transform <- match.fun(transform)
#   if (is.character(inverse))
#     inverse <- match.fun(inverse)
#   structure(list(name = name, transform = transform, inverse = inverse,
#                  breaks = breaks, minor_breaks = minor_breaks, format = format,
#                  domain = domain), class = "trans")
# }


#' Compute one of the life years measures
#'
#' Provides an estimate for one of the following measures: years lost (Andersen, 2013), years lost/saved (Andersen, 2017), or
#' life years difference (Manevski, Ruzic Gorenjec, Andersen, Pohar Perme, 2022).
#'
#' The life years difference (\code{measure='yd'}) is taken by default. If other
#' measures are of interest, use the \code{measure} argument.
#'
#' The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with the \code{rmap} argument. For example, if
#' age is in years in the data but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.
#'
#' Numerical integration is performed, argument
#' precision is set with argument \code{precision}, which defaults to 30-day
#' intervals for intergration. For higher accuracy take a smaller value (e.g. precision=1 makes
#' the integration on a daily basis).
#'
#' The observed curves are reported at event and censoring times. The
#' population curves are reported at all times used for the numerical integration.
#' Note that for the years lost (Andersen, 2013) measure, only the excess absolute risk is reported.
#'
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, \code{~1} specified on the right.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param measure choose which measure is used: 'yd' (life years difference; Manevski, Ruzic Gorenjec, Andersen, Pohar Perme, 2022), 'yl2017' (years lost/saved; Andersen 2017),
#' 'yl2013' (years lost/saved; Andersen 2013).
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @param var.estimator Choose the estimator for the variance ('none', 'bootstrap', 'greenwood'). Default is 'none'.
#' The 'greenwood' option is possible only for \code{measure='yd'}.
#' @param B if \code{var.estimator} is 'bootstrap'. The number of bootstrap replications. Default is 100.
#' @param precision precision for numerical integration of the population curve. Default is 30 (days).
#' The value may be decreased to get a
#' higher precision or increased to achieve a faster calculation.
#' @param add.times specific times at which the curves should be reported.
#' @param na.action a missing-data filter function. Default is \code{na.omit}.
#' @param conf.int the confidence level for a two-sided confidence interval. Default is 0.95.
#' @param timefix the timefix argument in survival::survfit.formula. Default is FALSE.
#' @param is.boot if TRUE, the function \code{years} has been called during a bootstrap replication.
#' @param first.boot if TRUE, this is the first bootstrap replication.
#' @return A list containing the years measure, the observed and population curves (or the excess curve for Andersen 2013).
#' The values are given as separate data.frames through time. Times are given in days, all areas are given in years.
#' For \code{measure='yl2017'} values are reported only at the last time point.
#' Functions \code{plot_f} and \code{plot_years} can be then used for plotting.
#' @seealso \code{\link{plot_f}}, \code{\link{plot_years}}
#' @examples
#'
#' library(relsurv)
#' # Estimate the life years difference for the rdata dataset.
#' mod <- years(Surv(time, cens)~1, data=rdata, measure='yd', ratetable=slopop,
#'              rmap=list(age=age*365.241), var.estimator = 'none')
#' # Plot the absolute risk (observed and population curve):
#' plot_f(mod)
#' # Plot the life years difference estimate:
#' plot_years(mod, conf.int=FALSE)
years <- function(
    formula=formula(data),
    data,
    measure=c('yd', 'yl2017', 'yl2013'),
    # estimator=c("F_P_final"),#, "F_P_Spi", "F_P_Spi2", "F_P", "F_P2", "all"),
    ratetable=relsurv::slopop,
    rmap,
    var.estimator=c('none', 'bootstrap', 'greenwood'),
    B=100,
    precision=30,
    add.times,
    na.action=stats::na.omit,
    conf.int=0.95,
    timefix=FALSE,
    # admin.cens,
    # cause.val,
    is.boot=FALSE,
    first.boot=FALSE
    # ,estimator.observed='Kaplan-Meier'
){

  # OLD ARGUMENTS:
  # F_P_Spi: Tako kot F_P_final, ignorira censoring. Ali pa vzame samo admin cens
  # F_P_Spi2: Vzame ves censoring
  # @param cause.val for competing risks, to be added.
  # @param admin.cens if a Date is supplied, administrative censoring is taken into account at that time
  # in the population curve. Works only if there's late entry, e.g. if the formula is \code{Surv(start,stop,event)~1}.

  ############ #
  # PREPARE OBJECTS:
  ############ #

  estimator=c("F_P_final") #  #' @param estimator which estimator should be used for calculating
  # estimator <- match.arg(estimator)
  arg.example <- FALSE # @param arg.example temporary argument, used for checking additionalities.

  Call <- match.call()
  if(!missing(rmap) & !is.boot & !first.boot)  rmap <- substitute(rmap)
  measure <- match.arg(measure)
  var.estimator <- match.arg(var.estimator)

  if(var.estimator=='bootstrap'){
    bootstrap <- TRUE
  } else if(var.estimator %in% c('none', 'greenwood')){
    bootstrap <- FALSE
  } else{
    stop('Incorrect value provided in argument var.estimator.')
  }

  if(!is.data.frame(data)) stop('Argument data is not a data.frame object.')
  data <- as.data.frame(data)

  out <- NULL
  late.values <- FALSE

  # These were arguments. To be deleted?
  exact.hazards <- FALSE # calculate hazards on a daily basis (to be checked)
  find.cond.time <- FALSE # if TRUE, return time at which there are at least 5 individuals in the at-risk set.

  # if(!missing(cause.val)){
  #   data$status <- ifelse(data$cause == cause.val, 1, 0)
  #   # Remove NAs:
  #   eniNAs <- which(is.na(data$status))
  #   if(length(eniNAs)>0) data <- data[-eniNAs,]
  # }

  # data$age <- round(data$age*365.241)
  # data$stop <- round(data$stop*365.241)

  # If Surv(start,stop, event) (possibly + mstate)
  if_start_stop <- length(as.character(formula[[2]])) %in% c(4,5)

  if(if_start_stop){
    start_col <- as.character(formula[[2]])[2]
    stop_col <- as.character(formula[[2]])[3]
    status_col <- as.character(formula[[2]])[4]
    starting_age <- as.vector(as.matrix(data[, start_col]))
  } else{
    stop_col <- as.character(formula[[2]])[2]
    if(!(stop_col %in% colnames(data))){
      stop(paste0('Instead of \'', stop_col, '\', please use a column from the data in the formula.'))
    }
  }

  # Check if no. at risk falls to zero at some point:
  if(if_start_stop){
    # Prepare at-risk matrix:
    find_tajms <- unique(sort(c(data[,start_col], data[,stop_col])))
    mat <- lapply(1:nrow(data), function(x) ifelse((data[x, start_col] < find_tajms) & (find_tajms <= data[x, stop_col]), 1, 0))
    mat2 <- matrix(unlist(mat), nrow = nrow(data), byrow = TRUE)
    # The sum of the individual at-risk processes:
    yi_left <- colSums(mat2)
    # If there's an empty at-risk at a later timepoint, censor the data:
    wh_yi <- which(yi_left==0)
    if(length(wh_yi)>1){
      if((!is.boot) & (!first.boot)){
        warning(paste0('In the time interval ', find_tajms[wh_yi[2]-1], '-', find_tajms[wh_yi[2]],
                       ' the at-risk sample is empty (nobody is followed). Survival cannot be estimated in this time interval.',
                       ' The data is censored at time ', find_tajms[wh_yi[2]-1], '.'))
      }
      # Censor data:
      data <- data[data[,start_col] <= find_tajms[wh_yi[2]-1], ]
      wh_cen <- which(data[, stop_col] > find_tajms[wh_yi[2]-1])
      data[wh_cen, stop_col] <- find_tajms[wh_yi[2]-1]
      data[wh_cen, status_col] <- 0
      if(!missing(add.times)){
        if(any(add.times > find_tajms[wh_yi[2]-1])) add.times <- add.times[add.times<=find_tajms[wh_yi[2]-1]]
      }
    }
    rm(mat,mat2)
  }
  data_orig <- data

  # if(starting.time=="left.truncated"){
  # if(!missing(admin.cens)){
  #   if(!inherits(admin.cens, 'Date')) warning('Object of class Date should be supplied to admin.cens.')
  #   end_date <- data$year+(data$stop-data$age)
  #   if(any(end_date > admin.cens)) warning('There are events that occur after the date of administrative censoring. Please check the values in arguments data and admin.cens.')
  #   id_admin_cens <- which(admin.cens==end_date)
  # }
  # }

  if(if_start_stop){
    starting.time <- 'left.truncated'
  } else{
    starting.time <- 'zero'
  }

  # Starting age
  starting_age <- rep(0,nrow(data))
  if(if_start_stop){
    starting_age <- as.vector(as.matrix(data[, start_col]))
  }
  starting_age <- as.numeric(starting_age)

  ############ #
  # YEARS ON DATA - GENERAL:
  ############ #

  surv_obj <- as.character(formula[[2]])

  if(missing(formula)){
    stop('Missing formula argument value.')
  } else{
    if('mstate' %in% surv_obj){
      juh <- 1:nrow(data)
      mod <- survival::survfit.formula(as.formula(Reduce(paste, deparse(formula))), data=data, timefix=timefix, id = juh, na.action=na.action)
    } else{
      mod <- survival::survfit.formula(formula, data=data, timefix=timefix, na.action=na.action)
    }
  }

  if('mstate' %in% surv_obj){
    surv_obj_new <- paste0(surv_obj[1], '(', surv_obj[2], ',', surv_obj[3])
    if(length(surv_obj)==5){
      surv_obj_new <- paste0(surv_obj_new, ',', surv_obj[4], ')')
    } else{
      surv_obj_new <- paste0(surv_obj_new, ')')
    }
    formula <- paste0(surv_obj_new, '~1')
  }
  status_obj <- surv_obj[length(surv_obj)]

  # if(!missing(cause.val)){
  #   mod$n.risk <- mod$n.risk[,1]
  #   mod$n.event <- mod$n.event[,cause.val+1]
  #   mod$surv <- 1-mod$pstate[,cause.val+1]
  #   mod$std.err <- mod$std.err[,cause.val+1]
  #   mod$cumhaz <- mod$cumhaz[,cause.val]
  # }

  if(!missing(add.times)){
    mod_sum <- summary(mod, times = sort(unique(c(mod$time, add.times))))
    if(any(!(add.times %in% mod_sum$time))){
      if(!is.boot){
        if(!first.boot){
          warning('Some values in add.times are after the last follow-up time. All measures are extrapolated up to these times. Please consider removing them.')
        }
        late.values <- TRUE

        miss_tajms <- add.times[!(add.times %in% mod_sum$time)]
        mod_sum$time <- c(mod_sum$time, miss_tajms)
        mod_sum$n.risk <- c(mod_sum$n.risk, rep(mod_sum$n.risk[length(mod_sum$n.risk)], length(miss_tajms)))
        mod_sum$n.event <- c(mod_sum$n.event, rep(0, length(miss_tajms)))
        mod_sum$surv <- c(mod_sum$surv, rep(mod_sum$surv[length(mod_sum$surv)], length(miss_tajms)))
        mod_sum$cumhaz <- c(mod_sum$cumhaz, rep(mod_sum$cumhaz[length(mod_sum$cumhaz)], length(miss_tajms)))

        # First fix std.err:
        if(is.nan(mod_sum$std.err[length(mod_sum$std.err)])){
          mod_sum$std.err[length(mod_sum$std.err)] <- mod_sum$std.err[length(mod_sum$std.err) - 1]
        }
        mod_sum$std.err <- c(mod_sum$std.err, rep(mod_sum$std.err[length(mod_sum$std.err)], length(miss_tajms)))
      }
    }
    mod$time <- mod_sum$time
    mod$n.risk <- mod_sum$n.risk
    mod$n.event <- mod_sum$n.event
    mod$surv <- mod_sum$surv
    mod$std.err <- mod_sum$std.err
    mod$cumhaz <- mod_sum$cumhaz
  }

  if(find.cond.time) return(mod$time[which.min(mod$n.risk<5)])

  # Calculate AUC:
  if(length(mod$time)>1){
    if(if_start_stop){
      survs <- c(1, mod$surv[1:(length(mod$surv)-1)])
      t_diff <- diff(c(mod$time[1], mod$time))
    } else{
      survs <- mod$surv
      t_diff <- diff(c(0, mod$time))
    }
    auc_data <- sum(t_diff*(1 - survs))
    auc_data_vec <- cumsum(t_diff*(1 - survs))
  } else{
    auc_data <- mod$time*mod$surv
    auc_data_vec <- auc_data
  }

  out$F_data <- 1-mod$surv
  out$auc_data <- auc_data/365.241
  out$auc_data_vec <- auc_data_vec/365.241

  # Exact hazards:
  if(exact.hazards){
    mod$time <- seq(min(mod$time), max(mod$time), by=1)
    mod$surv <- exp(-cumsum(rep(ratetable[1,1,1], max(mod$time)-min(mod$time)+1)))

    out$F_data <- 1-exp(-cumsum(c(0, rep(ratetable[1,1,1], max(mod$time)-min(mod$time)))))
    out$auc_data <- sum(out$F_data)/365.241
  }

  ############ #
  # SEPARATE YEARS FOR EVERY MEASURE:
  ############ #

  if(measure %in% c('yl2017', 'yl2013')){
    # YL_P preparation:
    data_yi <- data
    rform <- rformulate(formula, data, ratetable, na.action=na.action, rmap = rmap)
    data <- rform$data
    if(if_start_stop){
      if(!(start_col %in% colnames(data))){
        data[,start_col] <- data_orig[, start_col]
      }
    }

    # Check covariates:
    p <- rform$m
    if (p > 0) stop("There shouldn't be any covariates in the formula. This function gives non-parametric estimates of the hazards.")
    else data$Xs <- rep(1, nrow(data)) #if no covariates, just put 1

    out_n <- table(data$Xs) #table of strata
    out$time <- out$haz.excess <- out$haz.pop <- out$std.err <- out$strata <-  NULL

    kt <- 1 # the only stratum
    inx <- which(data$Xs == names(out_n)[kt]) #individuals within this stratum
    # tis <- sort(unique(rform$Y[inx])) #unique times
    if(!if_start_stop){
      tis <- rform$Y[inx] #unique times
      tis_seq <- seq(0, max(rform$Y[inx]), precision)
    } else{
      tis <- sort(unique(c(rform$Y[inx], data[, start_col]))) #unique times
      tis_seq <- seq(min(data[, start_col]), max(rform$Y[inx], data[, start_col]), precision)
    }

    if(!is.boot){
      tis <- sort(unique(c(tis, tis_seq)))
    }

    if(!missing(add.times)){
      tis <- sort(unique(c(tis, add.times)))
    }
    ltis <- length(tis)

    # Fix demographic covariates:
    if(if_start_stop){
      rform$R[,"year"] <- rform$R[,"year"] - rform$R[,"age"]
      rform$R[,"age"] <- 0
    }

    if(measure == 'yl2017'){
      # YL_O (used only for yl2017):
      if(if_start_stop){
        it_auc <- rep(NA, nrow(data_orig))
        mod_sum <- summary(mod, times=tis) # unique(sort(c(data_orig[,start_col], data_orig[,stop_col])))
        lsurv <- length(mod_sum$surv)
        val_mat <- matrix(0, nrow=nrow(data_orig), ncol=lsurv)
        for(it in 1:nrow(data_orig)){
          it_wh <- which(data_orig[it, start_col] == mod_sum$time)
          it_surv <- mod_sum$surv[it_wh:lsurv]/mod_sum$surv[it_wh]
          it_auc[it] <- sum(c(0, diff(mod_sum$time[it_wh:lsurv]))*(1 - it_surv))/365.241
          val_mat[it, it_wh:lsurv] <- cumsum(c(0, diff(mod_sum$time[it_wh:lsurv]))*(1 - it_surv))/365.241
        }
        # spodaj <- mod_sum$n.risk + cumsum(mod_sum$n.event) + cumsum(mod_sum$n.censor)
        YL_O_vec <- colMeans(val_mat) # colSums(val_mat)/spodaj
        YL_O <- mean(it_auc)
        F_O_time <- mod_sum$time
        F_O_ext <- data.frame(time=F_O_time, area=YL_O_vec)
        # Subset:
        F_O_ext2 <- subset(F_O_ext, time %in% mod$time)
        F_O_time <- F_O_ext2$time
        YL_O_vec <- F_O_ext2$area
      } else{
        YL_O_vec <- out$auc_data_vec
        YL_O <- out$auc_data
        F_O_time <- mod$time
        if(!(0 %in% F_O_time)){
          F_O_time <- c(0, F_O_time)
          YL_O_vec <- c(0, YL_O_vec)
        }
        # Prepare extended F_O object:
        if(0 %in% mod$time){
          F_O_temp <- data.frame(time=mod$time, surv=mod$surv)
        } else{
          F_O_temp <- data.frame(time=c(0, mod$time), surv=c(1, mod$surv))
        }
        F_O_ext <- data.frame(time=tis)
        F_O_ext <- merge(F_O_ext, F_O_temp, by='time', all.x=TRUE)
        F_O_ext$surv <- mstateNAfix(F_O_ext$surv, 0)

        tis_diff <- diff(c(0, F_O_ext$time))
        F_O_ext$area <- cumsum(tis_diff*(1 - F_O_ext$surv))/365.241
        F_O_ext <- F_O_ext[,c('time', 'area')]
      }

      F_O <- data.frame(time=F_O_time, area=YL_O_vec)

      ###
      # YL_P continue:
      it_auc_P <- rep(NA, nrow(data))
      it_auc_P_mat <- matrix(0, nrow=nrow(data), ncol=ltis)

      for(it in 1:nrow(data)){
        temp <- exp_prep(rform$R[it,,drop=FALSE],max(rform$Y),rform$ratetable,rform$status[it],times=tis,fast=FALSE, cmp=FALSE, ys=starting_age[it], netweiDM = FALSE)
        if(if_start_stop){
          it_wh <- which(data[it, start_col] == tis)
          hazs <- temp$yidli[it_wh:ltis]
          hazs[1] <- 0
          cumhazs <- cumsum(hazs)
          F_P <- 1 - exp(-cumhazs)
          it_auc_P[it] <- sum(c(tis[it_wh], diff(tis[it_wh:ltis]))*c(0, F_P[1:(length(F_P)-1)]))/365.241
          it_auc_P_mat[it,it_wh:ltis] <- cumsum(c(0, diff(tis[it_wh:ltis]))*c(0, F_P[1:(length(F_P)-1)]))/365.241
        } else{
          # it_wh <- which(data$age[it] == tis)
          hazs <- temp$yidli[1:ltis]
          hazs[1] <- 0
          cumhazs <- cumsum(hazs)
          F_P <- 1 - exp(-cumhazs)
          it_auc_P[it] <- sum(c(0, diff(tis))*c(0, F_P[1:(length(F_P)-1)]))/365.241
          it_auc_P_mat[it,] <- cumsum(c(0, diff(tis))*c(0, F_P[1:(length(F_P)-1)]))/365.241
        }

      }
      YL_P <- mean(it_auc_P)
      F_P <- data.frame(time=tis, area=colMeans(it_auc_P_mat))

      yd_curve <- data.frame(time=tis, est=F_O_ext$area - F_P$area)

      # Bootstrap:
      if(bootstrap){
        data_b <- data_orig
        data_b$id <- 1:nrow(data_b)
        yl_boot <- ylboot(theta=ylboot.iter, data=data_b, id="id",
                          B=B, verbose=0, #all_times = all_times,
                          ratetable=ratetable#, add.times=add.times
                          , starting.time, estimator, precision,
                          add.times = add.times,
                          formula = formula,
                          rmap = rmap, measure=measure
        )
        if(ncol(yl_boot[[2]])>nrow(F_O)){
          varsincol <- colVars(yl_boot[[2]], na.rm=TRUE)^(1/2)
          varsincol_df <- data.frame(time=yl_boot[[4]], area.se=varsincol)
          varsincol_df <- varsincol_df[varsincol_df$time %in% F_O$time,]
          F_O$area.se <- varsincol_df$area.se
        } else{
          F_O$area.se <- colVars(yl_boot[[2]], na.rm=TRUE)^(1/2)
        }
        F_P$area.se <- colVars(yl_boot[[3]], na.rm=TRUE)^(1/2)
        yl_boot <- as.data.frame(t(yl_boot[[1]]))
        yd_curve$est.se <- (colVars(yl_boot, na.rm=TRUE))^(1/2)
      }

      # Add CI:
      if((!is.boot) & (!first.boot)){
        if(!is.null(yd_curve$est.se)){
          yd_curve$lower <- yd_curve$est - yd_curve$est.se*stats::qnorm(0.5+conf.int/2)
          yd_curve$upper <- yd_curve$est + yd_curve$est.se*stats::qnorm(0.5+conf.int/2)
        }
      }
      # Values to be reported:
      if((!is.boot) & (!first.boot)){
        if(if_start_stop){
          # Report only at last time point - the values until this time are not suitable to report:
          out <- list(years=utils::tail(yd_curve,1), F_O=utils::tail(F_O,1), F_P=utils::tail(F_P,1), measure=measure)
        } else{
          # Report full measures:
          out <- list(years=yd_curve, F_O=F_O, F_P=F_P, measure=measure)
        }
      } else{
        out <- list(years=yd_curve, F_O=F_O, F_P=F_P, measure=measure)
      }
      return(out)

    } else{ # measure == 'yl2013'
      temp <- exp_prep(rform$R[,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status,times=tis, fast=TRUE, cmp=FALSE, ys=starting_age)
      temp$yi[temp$yi==0] <- Inf

      # Calculate measures:
      haz.pop <- temp$yidli/temp$yi
      mod_tis <- summary(mod, times = tis)
      F_E <- cumsum(mod_tis$surv*(mod_tis$n.event/mod_tis$n.risk - haz.pop))
      ltis <- length(tis)

      # To be checked, doesn't work ok
      # # Var as in Pavlic2018:
      # F_E_st <- sapply(1:ltis, function(s){
      #   (sum(mod_tis$surv[s:ltis]*(mod_tis$n.event[s:ltis]/mod_tis$n.risk[s:ltis] - haz.pop[s:ltis]))/mod_tis$surv[s]) # *c(0, diff(tis[s:ltis]))  /365.241
      # })
      # # Klemnova:
      # F_Ese <- (cumsum((mod_tis$surv)^2*(1 - F_E_st)^2*((mod_tis$n.event)/(mod_tis$n.risk^2))*c(0, diff(tis)))/365.241)^(1/2)
      # surv_int <- rev(cumsum(rev(c(0, diff(tis))*c(1, mod_tis$surv[1:(length(mod_tis$surv)-1)])))/365.241)
      #
      # # Moja:
      # F_E_int <- rev(cumsum(rev(c(0, diff(tis))*c(0, F_E[1:(length(F_E)-1)])))/365.241)
      # F_Ese <- (cumsum((surv_int)^2*(1 - F_E_st)^2*((mod_tis$n.event)/(mod_tis$n.risk^2))*c(0, diff(tis)))/365.241)^(1/2)
      #
      # # Observed:
      # F_Ese <- (cumsum(surv_int^2*((mod_tis$n.event)/(mod_tis$n.risk^2))*c(0, diff(tis)))/365.241)^(1/2)
      #
      # # Predlog glede na Andersen 2013:
      # F_Ese <- (cumsum((surv_int^2*(mod_tis$n.event - temp$yidli) + F_E_int^2*temp$yidli)/(mod_tis$n.risk^2)*c(0, diff(tis)))/365.241)^(1/2)

      # Calculate measures:
      YL <- cumsum(F_E*c(0, diff(tis)))/365.241
      F_E_area <- cumsum(c(0, diff(tis))*c(0, F_E[1:(length(F_E)-1)]))/365.241
      F_E_df <- data.frame(time=tis, prob=F_E, area=F_E_area) # , prob.se=F_Ese
      yd_curve <- data.frame(time=tis, est=YL)

      # Bootstrap:
      if(bootstrap){
        data_b <- data_orig
        data_b$id <- 1:nrow(data_b)
        yl_boot <- ylboot(theta=ylboot.iter, data=data_b, id="id",
                          B=B, verbose=0, #all_times = all_times,
                          ratetable=ratetable#, add.times=add.times
                          , starting.time, estimator, precision,
                          add.times = add.times,
                          formula = formula,
                          rmap = rmap, measure=measure
        )

        # Calculate area.se:
        area.se <- yl_boot[[2]]
        for(itar in 1:nrow(yl_boot[[2]])){
          prob_tmp <- as.vector(as.matrix(yl_boot[[2]][itar,]))
          area_tmp <- cumsum(c(0, diff(tis))*c(0, prob_tmp[1:(length(prob_tmp)-1)]))/365.241
          area.se[itar,] <- area_tmp
        }
        area.se <- as.vector(colVars(area.se, na.rm=TRUE))

        F_E_df$prob.se <- (colVars(yl_boot[[2]], na.rm=TRUE))^(1/2)
        F_E_df$area.se <- area.se
        yl_boot <- as.data.frame(t(yl_boot[[1]]))
        yd_curve$est.se <- (colVars(yl_boot, na.rm=TRUE))^(1/2)
      }

      if((!is.boot) & (!first.boot)){
        if(!is.null(yd_curve$est.se)){
          yd_curve$lower <- yd_curve$est - yd_curve$est.se*stats::qnorm(0.5+conf.int/2)
          yd_curve$upper <- yd_curve$est + yd_curve$est.se*stats::qnorm(0.5+conf.int/2)
        }
      }

      out <- list(years=yd_curve, F_E=F_E_df, measure=measure)
      return(out)
    }
  } else{ # measure == 'yd'

    ################################################### #
    # CIF on population:

    data_yi <- data
    rform <- rformulate(formula, data, ratetable, na.action=na.action, rmap = rmap)
    data <- rform$data
    if(if_start_stop){
      if(!(start_col %in% colnames(data))){
        data[,start_col] <- data_orig[, start_col]
      }
    }

    # Check covariates:
    p <- rform$m
    if (p > 0) stop("There shouldn't be any covariates in the formula. This function gives non-parametric estimates of the hazards.")
    else data$Xs <- rep(1, nrow(data)) #if no covariates, just put 1

    out_n <- table(data$Xs) #table of strata
    out$time <- out$haz.excess <- out$haz.pop <- out$std.err <- out$strata <-  NULL

    kt <- 1 # the only stratum
    inx <- which(data$Xs == names(out_n)[kt]) #individuals within this stratum
    if(!if_start_stop) tis <- sort(unique(c(rform$Y[inx], seq(0, max(rform$Y[inx]), precision)))) #unique times
    else tis <- sort(unique(c(rform$Y[inx], data[, start_col], seq(min(data[, start_col]), max(rform$Y[inx], data[, start_col]), precision)))) #unique times

    if(!missing(add.times)){
      tis <- sort(unique(c(tis, add.times)))
    }

    # Fix demographic covariates:
    if(if_start_stop){
      rform$R[,"year"] <- rform$R[,"year"] - rform$R[,"age"]
      rform$R[,"age"] <- 0
    }

    ### #
    # Greenwood Variance of area (not F):
    # First prepare objects:
    mod_gw <- summary(mod, times = tis)
    gw_df <- data.frame(time=mod_gw$time, surv=mod_gw$surv, n.risk=mod_gw$n.risk, n.event=mod_gw$n.event)
    # Then calculate:
    times_all2 <- c(0, diff(gw_df$time))/365.241
    surv_all <- c(1, gw_df$surv[1:(length(gw_df$surv)-1)])
    auc_all <- cumsum(times_all2*surv_all)
    area_var <- sapply(1:length(auc_all), function(x) {
      numer <- gw_df$n.risk[1:x]*(gw_df$n.risk[1:x] - gw_df$n.event[1:x])
      numer[numer==0] <- Inf
      sum(((auc_all[x] - auc_all[1:x])^2*gw_df$n.event[1:x])/numer)
    })
    if(is.nan(area_var[length(area_var)])){
      area_var[length(area_var)] <- area_var[length(area_var)-1]
    }
    ### #

    if(estimator=='F_P' | estimator=="all"){
      # Prepare at-risk matrix:
      # browser()
      # mat <- lapply(1:nrow(data), function(x) ifelse((data$start[x] < tis) & (tis <= data$Y[x]), 1, NA))
      # mat2 <- matrix(unlist(mat), nrow = nrow(data_yi), byrow = TRUE)
      # # The sum of the individual at-risk processes:
      # yi_left <- colSums(mat2)
      # yi_left[yi_left == 0] <- Inf
      #
      # mat3 <- lapply(1:nrow(data), function(x) data$age[x] + c(0, diff(tis)))

      if(any(rform$Y[inx]<=starting_age)) browser()

      temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=TRUE, cmp=FALSE, ys=starting_age)
      # Fix at-risk process, if needed:
      temp$yi[temp$yi==0] <- Inf

      out$time <- c(out$time, tis)						#add times

      # Calculate hazards:
      haz.pop <- temp$yidli/temp$yi
      out$haz.pop <- c(out$haz.pop,haz.pop)
      out$cum.haz.pop <- cumsum(out$haz.pop)
      out$F_P <- 1-exp(-out$cum.haz.pop)
      out$auc_pop <- sum(c(tis[1], diff(tis))*c(0, out$F_P[1:(length(out$F_P)-1)]))/365.241

    }

    data_spi2 <- data
    if(estimator=='F_P_Spi2' | estimator=="all"){

      if(any(data_spi2$start>=data_spi2$Y)) browser()
      # Take into account censoring:
      exp.surv2 <- nessie_spi(Surv(start, Y, stat)~1, data=data_spi2, ratetable=ratetable,
                              tis=tis, starting.time=starting.time, include.censoring = TRUE,
                              arg.example)
      out$haz.pop.spi2 <- exp.surv2$yidli/exp.surv2$yi
      out$cum.haz.pop.spi2 <- cumsum(out$haz.pop.spi2)
      out$F_P_Spi2 <- 1-exp(-out$cum.haz.pop.spi2)
      out$auc_pop_Spi2 <- sum(c(tis[1], diff(tis))*c(0, out$F_P_Spi2[1:(length(out$F_P_Spi2)-1)]))/365.241
    }
    if(estimator=='F_P_Spi' | estimator=="all"){
      if(TRUE){ # (!missing(admin.cens))    - tega nimamo vec
        data_spi2$stat <- 1
        # data_spi2$stat[id_admin_cens] <- 0  # - tole ni bilo zakomentirano, ko smo imeli admin.cens
        exp.surv <- nessie_spi(Surv(start, Y, stat)~1, data=data_spi2, ratetable=ratetable,
                               tis=tis, starting.time=starting.time, include.censoring = TRUE,
                               arg.example)
      } else{
        # Don't take into account censoring:
        exp.surv <- nessie_spi(Surv(start, Y, stat)~1, data=data_spi2, ratetable=ratetable,
                               tis=tis, starting.time=starting.time, include.censoring = FALSE,
                               arg.example)
      }

      out$haz.pop.spi <- exp.surv$yidli/exp.surv$yi
      out$cum.haz.pop.spi <- cumsum(out$haz.pop.spi)
      out$F_P_Spi <- 1-exp(-out$cum.haz.pop.spi)
      out$auc_pop_Spi <- sum(c(tis[1], diff(tis))*c(0, out$F_P_Spi[1:(length(out$F_P_Spi)-1)]))/365.241
    }
    if(estimator=='F_P_final'){
      # Shift all to the end:
      if(if_start_stop) data_yi[,stop_col] <- max(data_yi[,stop_col])

      rform2 <- rform
      rform <- rformulate(formula, data_yi, ratetable, na.action=na.action, rmap = rmap)

      # Shift all to the end:
      if(!if_start_stop){
        rform$Y <- rep(max(rform$Y), length(rform$Y))
        rform$data[,"Y"] <- rform$Y
      }
      data <- rform$data
      if(if_start_stop){
        if(!(start_col %in% colnames(data))){
          data[,start_col] <- data_orig[, start_col]
        }
      }

      # Check covariates:
      p <- rform$m
      if (p > 0) stop("There shouldn't be any covariates in the formula. This function gives non-parametric estimates of the hazards.")
      else data$Xs <- rep(1, nrow(data)) #if no covariates, just put 1

      out$haz.pop2 <- NULL

      kt <- 1 # the only stratum
      inx <- which(data$Xs == names(out_n)[kt]) #individuals within this stratum

      # Fix demographic covariates:
      if(if_start_stop){
        rform$R[,"year"] <- rform$R[,"year"] - rform$R[,"age"]
        rform$R[,"age"] <- 0
      }

      if(any(starting_age>=rform$Y[inx])) browser()
      temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=FALSE, cmp=FALSE, ys=starting_age, netweiDM = TRUE)
      temp$sidliD[1] <- 0
      # temp$sisD[1] <- 1
      temp$sisD[temp$sisD==0] <- Inf
      haz.pop2 <- temp$sidliD/temp$sisD

      out$haz.pop2 <- c(out$haz.pop2, haz.pop2)
      out$cum.haz.pop2 <- cumsum(out$haz.pop2)
      out$F_P2 <- 1-exp(-out$cum.haz.pop2)
      out$auc_pop2 <- sum(c(tis[1], diff(tis))*c(0, out$F_P2[1:(length(out$F_P2)-1)]))/365.241

      out$sidli <- temp$sidli
      out$sis <- temp$sis

      # DODATEK:
      haz.pop.ves.cas <- temp$sidli
      haz.pop.ves.cas[1] <- 0
      haz.pop.ves.cas <- haz.pop.ves.cas/temp$sis
      out$cum.haz.pop.ves.cas <- cumsum(haz.pop.ves.cas)
      out$F_P_ves_cas <- 1 - exp(-out$cum.haz.pop.ves.cas)
      out$auc_pop_ves_cas <- sum(c(tis[1], diff(tis))*c(0, out$F_P_ves_cas[1:(length(out$F_P_ves_cas)-1)]))/365.241
    }
    if(estimator=='F_P2' | estimator=="all"){
      # Shift all to the end:
      if(if_start_stop) data_yi[,stop_col] <- max(data_yi[,stop_col])

      rform2 <- rform

      rform <- rformulate(formula, data_yi, ratetable, na.action=na.action, rmap = rmap)

      # Shift all to the end:
      if(!if_start_stop){
        rform$Y <- rep(max(rform$Y), length(rform$Y))
        rform$data[,"Y"] <- rform$Y
      }
      data <- rform$data
      if(if_start_stop){
        if(!(start_col %in% colnames(data))){
          data[,start_col] <- data_orig[, start_col]
        }
      }

      # Check covariates:
      p <- rform$m
      if (p > 0) stop("There shouldn't be any covariates in the formula. This function gives non-parametric estimates of the hazards.")
      else data$Xs <- rep(1, nrow(data)) #if no covariates, just put 1

      out$haz.pop2 <- NULL

      kt <- 1 # the only stratum
      inx <- which(data$Xs == names(out_n)[kt]) #individuals within this stratum

      # Fix demographic covariates:
      if(if_start_stop){
        rform$R[,"year"] <- rform$R[,"year"] - rform$R[,"age"]
        rform$R[,"age"] <- 0
      }

      if(any(starting_age>=rform$Y[inx])) browser()

      # temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=TRUE, cmp=FALSE, ys=0)
      temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=TRUE, cmp=FALSE, ys=starting_age)

      # Fix at-risk process, if needed:
      temp$yi[temp$yi==0] <- Inf

      # Calculate hazards:
      haz.pop2 <- temp$yidli/temp$yi

      out$haz.pop2 <- c(out$haz.pop2, haz.pop2)
      out$cum.haz.pop2 <- cumsum(out$haz.pop2)
      out$F_P2 <- 1-exp(-out$cum.haz.pop2)
      # out$auc_pop2 <- sum(c(tis[1], diff(tis))*out$F_P2)/365.241
      out$auc_pop2 <- sum(c(tis[1], diff(tis))*c(0, out$F_P2[1:(length(out$F_P2)-1)]))/365.241
    }

    ###
    # Bootstrap:
    if(bootstrap){
      # browser()
      data_b <- data_orig
      data_b$id <- 1:nrow(data_b)
      yl_boot <- ylboot(theta=ylboot.iter, data=data_b, id="id",
                        B=B, verbose=0, #all_times = all_times,
                        ratetable=ratetable#, add.times=add.times
                        , starting.time, estimator, precision,
                        add.times = add.times,
                        formula = formula,
                        rmap = rmap, measure=measure
      )
      L_OP <- yl_boot[[3]]
      F_boot <- yl_boot[[2]]
      yl_boot <- as.data.frame(t(yl_boot[[1]]))
    }
    ###
    estimator.orig <- estimator
    if(estimator=='F_P_final') estimator = 'F_P2'

    out$strata <- c(out$strata, length(tis))				#number of times in this strata
    names(out$strata) <-  names(out_n)
    out$strata <-  NULL
    out$auc <- c(auc_data=out$auc_data, auc_pop=out$auc_pop, auc_pop2=out$auc_pop2, auc_pop_Spi=out$auc_pop_Spi, auc_pop_Spi2=out$auc_pop_Spi2)

    if(estimator=='all'){
      F_P_final <- data.frame(time=out$time,F_P=out$F_P, F_P2=out$F_P2, F_P_Spi=out$F_P_Spi, F_P_Spi2=out$F_P_Spi2)
    } else if(estimator=='F_P'){
      F_P_final <- data.frame(time=tis,prob=out$F_P)
    } else if(estimator=='F_P2'){
      F_P_final <- data.frame(time=tis,prob=out$F_P2)
    } else if(estimator=='F_P_Spi'){
      F_P_final <- data.frame(time=tis,prob=out$F_P_Spi)
    } else if(estimator=='F_P_Spi2'){
      F_P_final <- data.frame(time=tis,prob=out$F_P_Spi2)
    }

    # YD through time:
    F_data_yd <- data.frame(time=mod$time, F_data=out$F_data)
    pop.times <- F_P_final$time[!(F_P_final$time %in% mod$time)]
    if(length(pop.times) > 0){
      F_data_yd_tmp <- data.frame(time=pop.times, F_data=NA)
      F_data_yd <- rbind(F_data_yd, F_data_yd_tmp)
      F_data_yd <- F_data_yd[order(F_data_yd$time),]
      F_data_yd$F_data <- mstateNAfix(F_data_yd$F_data, 0)
    }
    F_data_yd$var <- area_var

    yd_data <- cumsum(c(F_data_yd$time[1], diff(F_data_yd$time))*c(0, F_data_yd$F_data[1:(nrow(F_data_yd)-1)]))/365.241

    # Population part:
    F_P_yd <- F_P_final
    yd_pop <- cumsum(c(F_P_yd$time[1], diff(F_P_yd$time))*c(0, F_P_yd$prob[1:(nrow(F_P_yd)-1)]))/365.241
    yd_curve <- data.frame(time=F_data_yd$time, yd=yd_data - yd_pop,
                           obs_var=F_data_yd$var,
                           # obs_var22=obs_var_time22,
                           yd_data=yd_data,
                           yd_pop=yd_pop
    )
    ###
    # Greenwood for prob:
    greenwood_est <- (mod$surv^2*cumsum(mod$n.event/((mod$n.risk - mod$n.event)*mod$n.risk)))^(1/2)
    # If Surv(t)=0 in the end, take the last var estimate:
    if(any(rev(mod$surv)==0)){
      greenwood_wh <- which(mod$surv==0)
      greenwood_est[greenwood_wh] <- greenwood_est[greenwood_wh[1]-1]
    }

    F_data_tmp <- data.frame(time=mod$time,
                             prob=out$F_data,
                             prob.se=greenwood_est,
                             area=NA,
                             area.se=NA)
    # Add values at time zero:
    F_tmp <- F_data_tmp[1,]
    F_tmp$time <- min(starting_age)
    F_tmp$prob <- 0
    F_tmp$prob.se <- 0
    if(!(F_tmp$time %in% F_data_tmp$time)) F_data_tmp <- rbind(F_tmp, F_data_tmp)

    if(!if_start_stop){
      F_P_final_tmp <- F_P_final[1,]
      F_P_final_tmp$time <- min(starting_age)
      F_P_final_tmp$prob <- 0
      if(!(F_P_final_tmp$time %in% F_P_final$time)) F_P_final <- rbind(F_P_final_tmp, F_P_final)
    }

    yd_curve_tmp <- yd_curve[1,]
    yd_curve_tmp$time <- min(starting_age)
    yd_curve_tmp[,2:ncol(yd_curve_tmp)] <- 0
    if(!(yd_curve_tmp$time %in% yd_curve$time)) yd_curve <- rbind(yd_curve_tmp, yd_curve)

    # Bootstrap:
    if(bootstrap){
      yd_curve$boot_var <- colVars(yl_boot, na.rm=TRUE)
      if(late.values){
        last_val <- utils::tail(yd_curve$boot_var[!is.na(yd_curve$boot_var)],1)
        yd_curve$boot_var[is.na(yd_curve$boot_var)] <- last_val
      }
      yl_sd_boot <- stats::sd(yl_boot[, ncol(yl_boot)], na.rm=TRUE)
    }

    # Add areas:
    F_data_tmp$area <- yd_curve$yd_data[yd_curve$time %in% F_data_tmp$time]
    F_P_final$area <- yd_curve$yd_pop#[yd_curve$time %in% F_P_final$time]
    F_data_tmp$area.se <- yd_curve$obs_var[yd_curve$time %in% F_data_tmp$time]^(1/2)

    # If, add boot variance:
    if(bootstrap & (!is.boot)){
      F_data_tmp$prob.se <- (F_boot$F_data[F_boot$time %in% F_data_tmp$time])^(1/2)
      F_P_final$prob.se <- (F_boot$F_P#[F_boot$time %in% F_P_final$time]
      )^(1/2)
      F_data_tmp$area.se <- L_OP$L_O[L_OP$time %in% F_data_tmp$time]^(1/2)
      F_P_final$area.se <- L_OP$L_P^(1/2)
    }

    # Column order:
    F_data_tmp <- F_data_tmp[, c('time', 'prob', 'area', 'prob.se', 'area.se')]

    # Choose relevant columns:
    if(bootstrap){
      yd_curve <- yd_curve[,c('time', 'yd', 'boot_var')]
    } else{
      yd_curve <- yd_curve[,c('time', 'yd', 'obs_var')]
    }
    yd_curve[,3] <- yd_curve[,3]^(1/2)
    colnames(yd_curve)[2:3] <- c('est', 'est.se')
    yd_curve$lower <- yd_curve$est - yd_curve$est.se*stats::qnorm(0.5+conf.int/2)
    yd_curve$upper <- yd_curve$est + yd_curve$est.se*stats::qnorm(0.5+conf.int/2)

    return_obj <- list(F_data=F_data_tmp,
                       F_P=F_P_final,
                       auc=out$auc,
                       yd_curve=yd_curve,
                       starting.time=starting.time,
                       estimator=estimator.orig,
                       out=out
    )

    if(bootstrap){
      return_obj[[length(return_obj)+1]] <- F_boot
      names(return_obj)[length(return_obj)] <- 'F_boot'
      return_obj[[length(return_obj)+1]] <- L_OP
      names(return_obj)[length(return_obj)] <- 'L_OP'
      return_obj <- append(return_obj, yl_sd_boot)
      names(return_obj)[length(return_obj)] <- 'yl_sd_boot'
    }

    return_short <- list(years=return_obj$yd_curve, F_O=return_obj$F_data, F_P=return_obj$F_P, measure=measure)

    if((bootstrap & (!is.boot)) #| ((!bootstrap) & (!is.boot))
    ){
      return_obj <- return_short
    }
    if((!bootstrap) & (!is.boot)){
      return_obj <- return_short
    }
    if(is.boot){
      return_obj <- return_short
    }

    if(var.estimator=='none'){
      return_obj$years <- return_obj$years[,1:2]
      find_cols <- (!grepl('.se', colnames(return_obj[[2]])))
      return_obj[[2]] <- return_obj[[2]][,find_cols]
      if(length(return_obj)==4){
        find_cols <- (!grepl('.se', colnames(return_obj[[3]])))
        return_obj[[3]] <- return_obj[[3]][,find_cols]
      }
    }

    return(return_obj)
  }
}

utils::globalVariables(c("time", "prob", "Curve", "est", "lower", "upper"))

# Bootstrap function:
ylboot <- function(theta, data, B = 5, id = "id", verbose = 0,
                   #all_times,
                   ratetable=relsurv::slopop, #add.times,
                   starting.time, estimator, precision,
                   add.times,
                   formula,
                   rmap, measure,
                   ...){
  ids <- unique(data[, id])
  n <- length(ids)
  if(!missing(add.times)){
    th <- ylboot.iter(formula, data, starting.time = starting.time, estimator = estimator, precision = precision,
                      ratetable=ratetable, first=TRUE, add.times = add.times, rmap = rmap, measure=measure, ...)
  } else{
    th <- ylboot.iter(formula, data, starting.time = starting.time, estimator = estimator, precision = precision,
                      ratetable=ratetable, first=TRUE, rmap = rmap, measure=measure, ...)
  }

  simple_par <- TRUE
  if(missing(add.times)) simple_par <- FALSE

  # Prepare objects:
  res <- data.frame(matrix(NA, nrow=B, ncol=nrow(th[[1]])))
  if(!missing(add.times)){
    add.times <- sort(unique(c(th[[1]]$time, add.times)))
  } else{
    add.times <- th[[1]]$time
  }

  Fdata <- data.frame(matrix(NA, nrow=B, ncol=length(add.times)))
  Fo <- data.frame(matrix(NA, nrow=B, ncol=nrow(th[[2]])))
  Fp <- data.frame(matrix(NA, nrow=B, ncol=length(add.times)))
  L_O <- data.frame(matrix(NA, nrow=B, ncol=length(add.times)))
  L_P <- data.frame(matrix(NA, nrow=B, ncol=length(add.times)))
  F_E <- data.frame(matrix(NA, nrow=B, ncol=length(add.times)))

  # Iteration:
  for (b in 1:B) {
    nek_obj <- ylboot.apply(formula, b, verbose, ids, data, id, add.times, starting.time, estimator, precision, ratetable, th, simple_par, rmap, measure, ...)
    res[b,1:length(nek_obj[[1]])] <- nek_obj[[1]]
    if(measure=='yl2013'){
      F_E[b,1:length(nek_obj[[2]])] <- nek_obj[[2]]
    }
    if(measure=='yl2017'){
      Fo[b,1:length(nek_obj[[2]])] <- nek_obj[[2]]
      Fp[b,1:length(nek_obj[[3]])] <- nek_obj[[3]]
    }
    if(measure=='yd'){
      subnek <- subset(nek_obj[[2]], time %in% add.times)

      sub_vec <- 1:nrow(subnek)
      Fdata[b,sub_vec] <- subnek$F_data
      Fp[b,sub_vec] <- subnek$F_P

      subnek2 <- subset(nek_obj[[3]], time %in% add.times)

      sub2_vec <- 1:nrow(subnek2)
      L_O[b,sub2_vec] <- subnek2$yd_data
      L_P[b,sub2_vec] <- subnek2$yd_pop
    }
  }
  res <- as.data.frame(t(res))
  if(measure == 'yl2013'){
    return(list(res, F_E))
  }
  if(measure == 'yl2017'){
    return(list(res, Fo, Fp, add.times))
  }
  else{
    if (verbose)
      cat("\n")

    F_obj <- data.frame(time=add.times,
                        F_data=colVars(Fdata, na.rm = TRUE),
                        F_P=colVars(Fp, na.rm = TRUE))
    L_OP <- data.frame(time=add.times,
                       L_O=colVars(L_O, na.rm = TRUE),
                       L_P=colVars(L_P, na.rm = TRUE))

    return(list(res, F_obj, L_OP))
  }
}

ylboot.apply <- function(formula, b, verbose, ids, data, id, add.times, starting.time, estimator, precision, ratetable, th, simple_par,
                         rmap, measure,
                         ...){

  if(starting.time=='left.truncated'){
    start_col <- as.character(formula[[2]])[2]
    stop_col <- as.character(formula[[2]])[3]
  } else{
    stop_col <- as.character(formula[[2]])[2]
  }

  if (verbose > 0) {
    cat("\nBootstrap replication", b, "\n")
  }
  bootdata <- NULL
  bids <- sample(ids, replace = TRUE)
  bidxs <- unlist(sapply(bids, function(x) which(x ==
                                                   data[, id])))
  bootdata <- data[bidxs, ]
  if (verbose > 0) {
    cat("applying theta ...")
  }

  if(length(unique(bootdata[,id]))==1){
    next
  }

  if(!missing(add.times) & simple_par){
    add.times.arg <- sort(unique(c(th[[1]]$time, add.times)))
  } else{
    add.times.arg <- th[[1]]$time
  }
  add.times.arg2 <- add.times.arg
  # Remove unnecessary times
  if(starting.time == 'left.truncated'){
    add.times.arg <- add.times.arg[add.times.arg<=max(bootdata[,stop_col])]
  } else{
    add.times.arg <- add.times.arg[add.times.arg<=max(bootdata[,stop_col])]# - bootdata[,start_col])]
  }

  thstar <- ylboot.iter(formula, bootdata, starting.time = starting.time, estimator = estimator, precision = precision,
                        ratetable=ratetable, add.times=add.times.arg, rmap=rmap, measure=measure, ...)
  if(measure == 'yl2013'){
    return(list(thstar[[1]]$est, thstar[[2]]$prob))
  }
  if(measure == 'yl2017'){
    FoO <- thstar[[2]]
    FpP <- thstar[[3]]
    thstar <- thstar[[1]]
    # if(nrow(th[[1]]) != nrow(thstar)) browser()

    if(nrow(FoO) < nrow(th[[2]])){
      mis.tajms <- th[[2]]$time[!(th[[2]]$time %in% FoO$time)]
      mis.tajms <- mis.tajms[mis.tajms <= max(FoO$time)]
      temp_df <- data.frame(time=mis.tajms, area=NA)
      FoO <- rbind(FoO, temp_df)
      FoO <- FoO[order(FoO$time),]
      FoO$area <- mstateNAfix(FoO$area, 0)
    }

    if(nrow(th[[1]]) < nrow(thstar)){
      thstar <- thstar[thstar$time %in% th[[1]]$time, ]
      FpP <- FpP[FpP$time %in% th[[1]]$time, ]
      foO <- foO[foO$time %in% th[[1]]$time, ]
    }

    if(length(th[[1]]$time[th[[1]]$time <= max(thstar$time)]) != length(thstar$time)) browser()
    pogoj <- any(th[[1]]$time[th[[1]]$time <= max(thstar$time)] != thstar$time)
    if(pogoj){
      missing_times <- th[[1]]$time[which(!(th[[1]]$time %in% thstar$time))]

      if(length(missing_times)>0){
        # There are times missing in thstar, add them:
        add_df <- thstar[1:length(missing_times),]
        add_df$time <- missing_times
        add_df$yd <- NA
        add_df$obs_var <- NA
        add_df$yd_data <- NA

        thstar <- rbind(thstar, add_df)
        thstar <- thstar[order(thstar$time),] # redundantno

        thstar$yd <- mstateNAfix(thstar$yd, 0)
        thstar$obs_var <- mstateNAfix(thstar$obs_var, 0)
        thstar$yd_data <- mstateNAfix(thstar$yd_data, 0)

        if(nrow(th[[1]]) < nrow(thstar)){
          thstar <- thstar[thstar$time %in% th[[1]]$time, ]
        }

        if(nrow(th[[1]]) != nrow(thstar)) browser()
      } else{
        # This means there's more times in thstar than needed. Remove unnecessary times:
        thstar <- thstar[-which(!(thstar$time %in% th[[1]]$time)),]
        FpP <- FpP[-which(!(FpP$time %in% th[[1]]$time)),]
        foO <- foO[-which(!(foO$time %in% th[[1]]$time)),]

        if(nrow(th[[1]]) != nrow(thstar)) browser()
      }
    }

    return(list(thstar$est, FoO$area, FpP$area))
  }

  L_OP <- thstar[[3]]
  Fobj <- thstar[[2]]
  thstar <- thstar[[1]]

  if(nrow(th[[1]]) < nrow(thstar)){
    thstar <- thstar[thstar$time %in% th[[1]]$time, ]
    L_OP <- L_OP[L_OP$time %in% th[[1]]$time, ]
    Fobj <- Fobj[Fobj$time %in% th[[1]]$time, ]
  }

  # Ali kaksne vrednosti manjkajo:
  if(length(th[[1]]$time[th[[1]]$time <= max(thstar$time)]) != length(thstar$time)) browser()
  pogoj <- any(th[[1]]$time[th[[1]]$time <= max(thstar$time)] != thstar$time)
  if(pogoj){

    missing_times <- th[[1]]$time[which(!(th[[1]]$time %in% thstar$time))]

    if(length(missing_times)>0){
      # There are times missing in thstar, add them:
      add_df <- thstar[1:length(missing_times),]
      add_df$time <- missing_times
      add_df$yd <- NA
      add_df$obs_var <- NA
      add_df$yd_data <- NA

      thstar <- rbind(thstar, add_df)
      thstar <- thstar[order(thstar$time),] # redundantno

      thstar$yd <- mstateNAfix(thstar$yd, 0)
      thstar$obs_var <- mstateNAfix(thstar$obs_var, 0)
      thstar$yd_data <- mstateNAfix(thstar$yd_data, 0)

      if(nrow(th[[1]]) < nrow(thstar)){
        thstar <- thstar[thstar$time %in% th[[1]]$time, ]
      }

      if(nrow(th[[1]]) != nrow(thstar)) browser()
    } else{
      # This means there's more times in thstar than needed. Remove unnecessary times:
      thstar <- thstar[-which(!(thstar$time %in% th[[1]]$time)),]
      L_OP <- L_OP[-which(!(L_OP$time %in% th[[1]]$time)),]
      Fobj <- Fobj[-which(!(Fobj$time %in% th[[1]]$time)),]

      if(nrow(th[[1]]) != nrow(thstar)) browser()
    }
  }

  # thstar$b <- b
  # Save result:
  # res[b,] <-

  list(thstar$est, Fobj, L_OP)
}

ylboot.iter <- function(formula, data, #all_times,
                        starting.time, estimator, precision,
                        ratetable=relsurv::slopop,
                        first=FALSE, add.times,
                        rmap, measure
){
  if(!missing(rmap))  rmap <- as.call(rmap)
  if(first){
    is.boot <- FALSE
    first.boot <- TRUE
  } else{
    is.boot <- TRUE
    first.boot <- FALSE
  }

  # Round, if needed:
  tolerance <- 1e-15

  if(missing(add.times)){
    object <- years(formula = formula, data = data, ratetable = ratetable,
                    precision=precision, var.estimator='greenwood', is.boot=is.boot, first.boot = first.boot, rmap = rmap, measure=measure)
    # estimator = estimator,
  } else{
    object <- years(formula = formula, data = data, ratetable = ratetable,
                    precision=precision, var.estimator='greenwood', add.times=add.times, is.boot=is.boot, first.boot = first.boot, rmap = rmap, measure=measure)
    # estimator = estimator,
  }
  if(measure=='yd'){
    if(first) return(list(object$years, object$F_O))
    else{
      # return(object$yd_curve)
      Fobj <- merge(object$F_P[,c('time','prob')], object$F_O[,c('time','prob')], by='time', all.x=TRUE)
      Fobj <- Fobj[,c(1,3,2)]
      colnames(Fobj)[2:3] <- c('F_data','F_P')

      L_OP <- merge(object$F_P[,c('time','area')], object$F_O[,c('time','area')], by='time', all.x = TRUE)
      L_OP <- L_OP[,c(1,3,2)]
      colnames(L_OP)[2:3] <- c('yd_data', 'yd_pop')

      return(list(object$years,
                  Fobj,
                  L_OP))
    }
  } else if(measure=='yl2013'){
    return(list(object$years, object$F_E))
  } else{
    return(list(object$years, object$F_O, object$F_P))
  }

}

plot_helper_years <- function(years, obj){

  df_poly <- data.frame(time=years[[obj]]$time/365.241,
                        prob=years[[obj]]$prob)
  df_st <- df_poly[1,]
  df_st$prob <- 0
  df_end <- df_poly[nrow(df_poly),]
  df_end$prob <- 0
  df_poly <- rbind(df_st, df_poly, df_end)

  df_poly
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Plot the absolute risk (observed and population curve)
#'
#' Plots the estimated observed and population curve for the
#' life years difference (Manevski, Ruzic Gorenjec, Andersen, Pohar Perme, 2022).
#'
#' A ggplot2 implementation for plotting the observed and population curves. The type of curves is
#' dependent upon the measure calculated using \code{years} function (argument \code{measure}).

#' @param years the object obtained using function \code{years}.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param xbreak the breaks on the x axis (this is supplied to \code{scale_x_continuous}).
#' @param ybreak the breaks on the y axis (this is supplied to \code{scale_y_continuous}).
#' @param xlimits define the limits on the x axis (this is supplied to \code{scale_x_continuous}).
#' @param ylimits define the limits on the y axis (this is supplied to \code{scale_y_continuous}).
#' @param show.legend if TRUE, the legend is shown on the graph.
#' @return A ggplot object
#' @seealso \code{\link{years}}, \code{\link{plot_years}}
#'
plot_f <- function(years, xlab='Time interval', ylab='Absolute risk', xbreak, ybreak, xlimits, ylimits, show.legend=TRUE){
  # years: object given from the years() function
  # xlab: define xlab
  # ylab: define ylab
  # xbreak: The breaks on x axis
  # ybreak: The breaks on y axis
  # xlimits: Define the limits on the x axis
  # ylimits: Define the limits on the y axis
  # show.legend: TRUE by default (shows the legend)

  # Checks:
  if(years$measure != 'yd'){
    stop("The plot_f function is available only for the YD measure (argument measure='yd' in the years function).")
  }

  out <- rbind(
    cbind(years$F_O[,c('time', 'prob')], Curve='Observed'),
    cbind(years$F_P[,c('time', 'prob')], Curve='Population')
  )

  if(missing(xlimits)){
    xlimits <- c(min(out$time), max(out$time))/365.241
  }
  if(missing(ylimits)){
    ylimits <- c(0,max(out$prob))*1.1
  }

  colorji <- gg_color_hue(3)
  colorji <- colorji[c(1,3)]

  g <- ggplot2::ggplot(out)+
    ggplot2::geom_step(ggplot2::aes(time/365.241, prob, color=Curve)#, size=1.001
    )+
    ggplot2::scale_color_manual(values=colorji)+
    ggplot2::xlab(xlab)+
    ggplot2::ylab(ylab)

  poly_data <- plot_helper_years(years, 'F_O')
  poly_P <- plot_helper_years(years, 'F_P')

  g <- g+
    pammtools::geom_stepribbon(ggplot2::aes(x=time/365.241, ymin=0, ymax=prob, fill=Curve), alpha=0.3, linetype='dashed')+
    ggplot2::scale_fill_manual(values = colorji)

  if(!missing(xbreak)){
    g <- g +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlimits, breaks = xbreak)
  } else{
    g <- g +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlimits)
  }

  if(!missing(ybreak)){
    g <- g +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits=ylimits, breaks = ybreak)
  } else{
    g <- g +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits=ylimits)
  }

  g <- g +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = 'bottom',
                   legend.title = ggplot2::element_blank())+
    ggplot2::theme(text = ggplot2::element_text(size=14))+
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.minor.x = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.major.y = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.minor.y = ggplot2::element_line(linetype='dashed', colour = 'grey85'))

  if(!show.legend){
    g <- g +
      ggplot2::theme(legend.position = 'none')
  }

  g
}

#' Plot the years measure
#'
#' Plot the years measure obtained from the \code{years} function.
#'
#' A ggplot2 implementation for plotting the years measure. The type of curve is
#' dependent upon the measure calculated using the \code{years} function (argument \code{measure}).

#' @param years the object obtained using function \code{years}.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param xbreak the breaks on the x axis (this is supplied to \code{scale_x_continuous}).
#' @param ybreak the breaks on the y axis (this is supplied to \code{scale_y_continuous}).
#' @param xlimits define the limits on the x axis (this is supplied to \code{scale_x_continuous}).
#' @param ylimits define the limits on the y axis (this is supplied to \code{scale_y_continuous}).
#' @param conf.int if TRUE, the confidence interval is plotted.
#' @param ymirror mirror the y values (w.r.t. the x axis).
#' @param yminus use function y -> -y when plotting.
#' @return A ggplot object
#' @seealso \code{\link{years}}, \code{\link{plot_f}}
#'
plot_years <- function(years, xlab='Time interval', ylab='Years', xbreak, ybreak, xlimits, ylimits, conf.int=FALSE, ymirror=FALSE, yminus=FALSE){

  out <- years$years

  if(conf.int){
    if(is.null(out$lower)){
      stop('Confidence intervals not present in the years object. Please set conf.int=FALSE or use the var.estimator argument in the years function.')
    }
  }

  if(years$measure=='yl2017' & nrow(out)==1){
    stop('The years measure is reported at the end of follow-up thus it is not plotted.')
  }

  if(yminus){
    out$est <- -out$est
    if(!is.null(out$lower)){
      tmp_lower <- out$lower
      out$lower <- -out$upper
      out$upper <- -tmp_lower
    }
  }

  if(missing(xlimits)){
    xlimits <- c(min(out$time[1]), max(out$time))/365.241
  }
  if(missing(ylimits)){
    tmp_vec <- out$est
    if(!is.null(out$lower)) tmp_vec <- c(out$est, out$lower, out$upper)
    ymax <- max(tmp_vec)
    ymin <- min(tmp_vec)
    ylimits <- c(ymin,ymax)*1.1
  }

  g <- ggplot2::ggplot(out)+
    ggplot2::geom_step(ggplot2::aes(time/365.241, est)#, size=1.001
                       )

  if(conf.int){
    g <- g+
      ggplot2::geom_step(ggplot2::aes(time/365.241, lower), linetype='dashed')+
      ggplot2::geom_step(ggplot2::aes(time/365.241, upper), linetype='dashed')
  }

  g <- g+
    ggplot2::xlab(xlab)+
    ggplot2::ylab(ylab)

  if(!missing(xbreak)){
    g <- g+
      ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlimits, breaks = xbreak)
  } else{
    g <- g+
      ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlimits)
  }

  # Helper:
  trans <- function(x) -x
  inv <- function(x) -x
  reverse_fun <- scales::trans_new(name = "reverse_new",
                                   transform = trans,
                                   inverse = inv
  )

  if(!missing(ybreak)){
    g <- g +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylimits, breaks = ybreak)
  } else{
    g <- g +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylimits)
  }
  if(ymirror){
    g <- g +
      ggplot2::coord_trans(y=reverse_fun)
  }

  g <- g +
    ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size=14))+
    ggplot2::expand_limits(y = 0)+
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.minor.x = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.major.y = ggplot2::element_line(linetype='dashed', colour = 'grey85'),
      panel.grid.minor.y = ggplot2::element_line(linetype='dashed', colour = 'grey85'))

  g

}

