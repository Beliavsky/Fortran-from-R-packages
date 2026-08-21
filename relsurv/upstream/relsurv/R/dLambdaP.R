dLambdaPR <- function(data, all_times, event_times, ratetable, atts){

  # Objects for exp_prep / expc:
  fk <- (atts$factor != 1)
  nfk <- length(fk)
  cuts <- atts$cutpoints

  ltype <- length(atts$type)
  rfac <- rep(0, ltype)
  for(i in 1:ltype){
    if(atts$type[i]==1) rfac[i] <- 1
  }
  rfac <- as.integer(rfac)
  # rfac <- ifelse(atts$type == 1, 1, 0)
  adim <- atts$dim
  acuts <- unlist(cuts)

  Yt_all <- Yt(data, all_times)

  ltimes <- length(all_times)
  nr <- nrow(data)

  data_m <- as.matrix(data)

  # outcome <- vector("list", length = ltimes)
  outcome <- matrix(0, nrow=nr, ncol=ltimes)

  for(i in 1:ltimes){
    if(i==1){
      if(all_times[i]==0) next

      tstart <- 0
    } else{
      tstart <- all_times[i-1]
    }

    tstop <- all_times[i]

    # ord_id <- order(data$Y)
    # data <- data[ord_id, ]

    wh_at_risk <- (Yt_all[[i]] == 1)

    data_tmp <- data_m[wh_at_risk, 4:(nfk + 3), drop = FALSE]

    data_tmp[,fk] <- data_tmp[,fk]+tstart

    # pop_survs <- exp_prep(data_tmp,
    #                       rep(tstop-tstart, sum(wh_at_risk)),
    #                       ratetable)

    times <- rep(tstop-tstart, sum(wh_at_risk))

    temp <- .Call("expc", rfac, adim,
                  acuts, ratetable, data_tmp, times,
                  PACKAGE = "relsurv")
    pop_survs <- temp$surv

    # for(j in 1:nr){
    #   if(Yt_all[[i]][j] == 1){
    #     outcome[j,i] <- 1 # pop hazard v tem intervalu
    #     xx <- exp_prep(data[j, 4:(nfk + 3), drop = FALSE], data$Y -
    #                      data$start, ratetable)
    #   }
    # }

    outcome[wh_at_risk, i] <- -log(pop_survs)
  }

  # Cum hazards:
  outcome <- t(apply(outcome, 1, cumsum))

  # Hazards at event times only:
  # whe <- which(all_times %in% event_times)
  # outcome <- outcome[, whe]

  # Add additional zeros for time 0 (for diffs in hazard):
  outcome <- cbind(rep(0, nrow(outcome)), outcome)

  return(outcome)
}

# exp_prep(x,y,ratetable,status,times,fast=FALSE,ys,prec,cmp=F,netweiDM=FALSE)
#
# jee <- exp_prep(rform$R, rform$Y, rform$ratetable, rform$status, all_times, fast=FALSE, data$start, 1, FALSE, netweiDM=FALSE)
#
# jee <- exp_prep(rform$R, rform$Y, rform$ratetable, rform$status, all_times[1:2], fast=FALSE, data$start, 1, FALSE, netweiDM=FALSE)

