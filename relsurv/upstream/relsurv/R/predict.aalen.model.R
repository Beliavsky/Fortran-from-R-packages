#' Subject-specific prediction for an Aalen model
#'
#' A function that returns subject-specific predictions from a survaalen or rsaalen object.
#' @param object A survaalen or rsaalen object
#' @param newdata A data.frame with one row (add covariate values in columns)
#' @param type By default ('Haz') returns cumulative hazards. If 'prob' it calculates CDF/CIFs. If 'both' a list is returned with Haz and prob.
#' @param ... Not used for now
#' @return A data.frame containing the times and survival (or excess/population) estimates.
#'
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @examples
#' data(rdata)
#' mod <- survaalen(Surv(time, cens)~sex+age, data=rdata)
#' predict(mod, data.frame(age=68, sex=2))
#'
#' mod <- rsaalen(Surv(time, cens)~sex+age, data=rdata, ratetable=slopop,
#'                 rmap=list(age=age*365.241))
#' predict(mod, data.frame(age=68, sex=2, year=as.Date('1982-06-24')))
predict.aalen.model <- function(object, newdata, type='Haz', ...){
  # Zaenkrat delamo z object$times. Treba bo dodati se add.times.

  type <- match.arg(type, c("Haz", "prob", "both"))

  # lin.pred <- object$coefficients[,2]
  # if(ncol(object$coefficients) >= 3){
  #   for(cov in colnames(object$coefficients)[3:length(colnames(object$coefficients))]){
  #     lin.pred <- lin.pred + newdata[1,cov]*object$coefficients[,cov]
  #   }
  # }

  kar <- deparse(object$formula[[2]])
  autkam <- gsub(' ', '', strsplit(substr(kar, start = 6, stop=nchar(kar)-1), ',')[[1]])
  newdata[autkam] <- 1

  # Calculate linear predictor:
  lin.pred <- as.vector(model.matrix(object, newdata[1,])%*%t(object$coefficients[,2:ncol(object$coefficients)]))


  # rsaalen:
  if('ratetable' %in% names(object)){
    newdata_2 <- data.frame(matrix(1, nrow=1, ncol=length(autkam)))
    if(length(autkam)==3) newdata_2[,1] <- 0
    colnames(newdata_2) <- autkam
    newdata_2 <- cbind(newdata, newdata_2)

    rform <- suppressWarnings(rformulate(object$formula,
                                         newdata_2, object$ratetable, stats::na.omit(), rmap = object$rmap))


    fk <- (attributes(rform$ratetable)$factor != 1)
    nfk <- length(fk)
    pop.surv <- sapply(1:nrow(object$coefficients), function(i) exp_prep(rform$data[1, 4:(nfk + 3), drop = FALSE], object$coefficients[i,1], object$ratetable))
    Haz.p <- -log(pop.surv)

    if(type=='Haz'){
      data.frame(time=object$coefficients[,1], Haz.e=lin.pred, Haz.p=Haz.p)
    } else if(type=='prob'){
      S_O <- exp(-(lin.pred + Haz.p))

      prob.e <- cumsum(S_O*diff(c(0, lin.pred)))
      prob.p <- cumsum(S_O*diff(c(0, Haz.p)))

      data.frame(time=object$coefficients[,1], prob.e=prob.e, prob.p=prob.p)
    } else{
      out <- list()
      out$Haz <- data.frame(time=object$coefficients[,1], Haz.e=lin.pred, Haz.p=Haz.p)

      # copied from up:
      S_O <- exp(-(lin.pred + Haz.p))

      prob.e <- cumsum(S_O*diff(c(0, lin.pred)))
      prob.p <- cumsum(S_O*diff(c(0, Haz.p)))
      out$prob <- data.frame(time=object$coefficients[,1], prob.e=prob.e, prob.p=prob.p)
      out

    }
  } else{
    # survaalen:
    if(type=='Haz'){
      data.frame(time=object$coefficients[,1], Haz=lin.pred)
    } else if(type=='prob'){
      data.frame(time=object$coefficients[,1], prob=1-exp(-lin.pred))
    } else{
      out <- list()
      out$Haz <- data.frame(time=object$coefficients[,1], Haz=lin.pred)
      out$prob <- data.frame(time=object$coefficients[,1], prob=1-exp(-lin.pred))
      out
    }

  }
}

