#' Subject-specific prediction for an excess hazards model
#'
#' A function that returns subject-specific predictions from an rsadd object.
#' @param object An rsadd object
#' @param newdata A data.frame with one row (add covariate values in columns)
#' @param type By default ('Haz') returns cumulative hazards. If 'prob' it calculates CIFs. If 'both' a list is returned with Haz and prob.
#' @param add.times specific times at which predictions should be obtained.
#' @param ... Not used for now
#' @return A data.frame containing the times, excess and population components.
#'
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @examples
#' data(rdata)
#' fit <- rsadd(Surv(time,cens)~sex+age,rmap=list(age=age*365.241),
#'              ratetable=slopop,data=rdata,int=5,method="EM")
#' predict(fit, data.frame(sex=2,age=68, year=as.Date('1982-06-24')))
predict.rsadd <- function(object, newdata, type='Haz', add.times,
                          ...){
  # Zaenkrat delamo z object$times. Treba bo dodati se add.times. Ta add.times bo verjetno v rsadd

  type <- match.arg(type, c("Haz", "prob", "both"))

  # lin.pred <- 0
  # for(cov in names(object$coefficients)){
  #   lin.pred <- lin.pred + sum(newdata[1,cov]*object$coefficients[cov])
  # }
  lin.pred <- sum(model.matrix(object, newdata[1,])*object$coefficients)

  # haz_function(object$formula, rdata, object$ratetable, rmap=list(age=age*365.241), add.times=0, include.all.times = FALSE)

  kar <- deparse(object$formula[[2]])
  autkam <- gsub(' ', '', strsplit(substr(kar, start = 6, stop=nchar(kar)-1), ',')[[1]])
  newdata_2 <- data.frame(matrix(1, nrow=1, ncol=length(autkam)))
  if(length(autkam)==3) newdata_2[,1] <- 0
  colnames(newdata_2) <- autkam
  newdata_2 <- cbind(newdata, newdata_2)

  if(is.null(object$rmap)){
    rform <- suppressWarnings(rformulate(object$formula,
                                         newdata_2,object$ratetable))
  } else{
    rform <- suppressWarnings(rformulate(object$formula,
                                         newdata_2,object$ratetable, rmap = object$rmap))
  }

  fk <- (attributes(rform$ratetable)$factor != 1)
  nfk <- length(fk)

  # Include add.times:
  all_times <- object$times
  if(!missing(add.times)){
    add.times <- add.times[add.times<=all_times[length(all_times)]]
    add.times <- add.times[!(add.times %in% all_times)]
    all_times <- sort(unique(c(add.times, all_times)))
  }

  pop.surv <- sapply(1:length(all_times), function(i) exp_prep(rform$data[1, 4:(nfk + 3), drop = FALSE], all_times[i], object$ratetable))
  Haz.p <- -log(pop.surv)

  # Expand Lambda0:
  if (!missing(add.times)) {
    df_tmp <- data.frame(time = object$times, Lambda0 = object$Lambda0)
    df_novo <- data.frame(time = add.times, Lambda0 = NA)
    df_add <- rbind(df_tmp, df_novo)
    df_add <- df_add[order(df_add$time), ]
    df_add$Lambda0 <- mstateNAfix(df_add$Lambda0,
                                            0)
    rownames(df_add) <- as.character(1:nrow(df_add))

    Haz.e <- df_add$Lambda0 * exp(lin.pred)
  } else{
    Haz.e <- object$Lambda0 * exp(lin.pred)
  }

  if(type=='Haz'){
    data.frame(time=all_times, Haz.e=Haz.e, Haz.p=Haz.p)
  } else if(type=='prob'){
    S_O <- exp(-(Haz.e + Haz.p))

    prob.e <- cumsum(S_O*diff(c(0, Haz.e)))
    prob.p <- cumsum(S_O*diff(c(0, Haz.p)))

    data.frame(time=all_times, prob.e=prob.e, prob.p=prob.p)
  } else{
    out <- list()
    out$Haz <- data.frame(time=all_times, Haz.e=Haz.e, Haz.p=Haz.p)

    S_O <- exp(-(Haz.e + Haz.p))

    prob.e <- cumsum(S_O*diff(c(0, Haz.e)))
    prob.p <- cumsum(S_O*diff(c(0, Haz.p)))
    out$prob <- data.frame(time=all_times, prob.e=prob.e, prob.p=prob.p)
    out
  }
}
