#' bond_future_price
#'
#' @param call_prices a vector of call prices on a bond futures, in numeric format
#' @param call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param put_prices a vector of put prices on a bond futures, in numeric format
#' @param put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param nb_log a number for the number of component densities in the lognormal mixture to model the bond futures' price, either 2 or 3, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param day_count_conv a number for the day count convention, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot a number for the options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param ctd_matu a date for the maturity date of the current Cheapest-to-Deliver Bond of the bond futures, in Date format
#' @param fut_price a number for the bond futures' price on calibration date, in numeric format
#' @param fut_matu a date for the bond futures' maturity date, in Date format
#' @param option_matu a date for the maturity date of the options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param nationality a character for the nationality of the issuer of the bond underlying the futures contract, in character format (NA by default)
#' @param currency a character for the currency in which the futures contract and the options are traded, in character format (NA by default)
#'
#' @returns the mean and standard deviation of each component lognormal density and the weight on the first density (for a mixture of 2) or on the first 2 densities (for a mixture of 3) in numeric format,  a discretized domain of the futures price at option's maturity by increments of 0.001 units of the currency of the futures contract, in numeric format, the probability density for each value in the discretized domain in numeric format, the cumulative density for each value in the discretized domain in numeric format, the type of convergence in the optimization in numeric format (0 indicating successful convergence), the mean, standard deviation, skewness and kurtosis of the distribution of the futures price in numeric format, quantiles of order 0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 90%, 95%, 99%, 99.5% and 99.9% of the distribution in numeric format, the mode of the distribution in numeric format, the options' prices predicted by the model in numeric format, a plot of the RND and a plot of the CDF of the futures price
#' @export
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm
#' @importFrom utils head tail
#' @import dplyr
#' @import lubridate
#' @import zoo
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' bond_future_price(c(10.39,9.92,9.46,9.00,8.55,8.10,7.66,7.23,
#' 6.81,6.39,5.98,5.58,5.20,4.82,4.46,4.10,3.76,3.44,3.13,2.83,2.56,
#' 2.29,2.05,1.82,1.61,1.42,1.25,1.09,0.95,0.82,0.71,0.61,0.53,0.45,
#' 0.38,0.33,0.28,0.23,0.20,0.17,0.14,0.12,0.10,0.08),
#' seq(106, 127.5, 0.5),
#' c(0.22,0.25,0.29,0.33,0.38,0.43,0.49,0.56,0.64,0.72,0.81,0.91,
#' 1.03,1.15,1.29, 1.43,1.59,1.77,1.96,2.16,2.39,2.62,2.88,3.15,
#' 3.44,3.75,4.08, 4.42,4.78,5.15,5.54,5.94,6.36,6.78,7.21,7.66,
#' 8.11,8.56,9.03, 9.50,9.97,10.45,10.93,11.41),
#' seq(106, 127.5, 0.5),
#' 2,
#' 0.0344,
#' 1,
#' 3,
#' as.Date("2033-11-01"),
#' 116.17,
#' as.Date("2024-12-10"),
#' as.Date("2024-11-22"),
#' as.Date("2024-06-14"),
#' "Italian",
#' "EUR")
#' }
#'

bond_future_price <- function(call_prices, call_strikes, put_prices, put_strikes, nb_log, r, day_count_conv, cot,
                              ctd_matu, fut_price, fut_matu, option_matu, start_date, nationality = NA, currency = NA){

  if(length(nb_log) == 1 & length(r) == 1 & length(day_count_conv) == 1 & length(cot) == 1 & length(ctd_matu) == 1 &
     length(fut_price) == 1 & length(fut_matu) == 1 & length(option_matu) == 1 & length(start_date) == 1 &
     length(nationality) == 1 & length(currency) == 1 & length(call_prices) > 1 & length(call_strikes) > 1 &
     length(put_prices) > 1 & length(put_strikes) > 1){

    bond_fut <- data.frame(ctd_matu, fut_price, option_matu, start_date, fut_matu, nationality, currency) %>%
      rename_with(~c("ctd_matu", "fut_price", "option_matu", "start_date", "fut_matu", "nationality", "currency"))

    if(bond_fut$start_date < bond_fut$option_matu & bond_fut$option_matu <= bond_fut$fut_matu){

      if(day_count_conv == 1){
        bond_fut <- bond_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                          as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")))
      } else if(day_count_conv == 2){
        bond_fut <- bond_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/360)
      } else if(day_count_conv == 3){
        bond_fut <- bond_fut %>% mutate(option_term =as.numeric(option_matu - start_date)/365)
      } else {
        bond_fut <- bond_fut %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                        stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                        plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                           ceiling_date(start_date, "months")))/30),
                                        option_term = (stub_1 + stub_2 + plain_months*30)/360)}

      call <- function(x, KC){
        d1_C <- (x[1] + x[2]^2 - log(KC))/x[2]
        d2_C <- d1_C - x[2]
        if(cot %in%c(1, 2)){call <- exp(-r*T)*(exp(x[1] + (x[2]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))
        } else(call <- exp(x[1] + (x[2]^2/2))*pnorm(d1_C) - KC*pnorm(d2_C))
      }

      call_mix <- function(x, KC){
        if(cot %in%c(1, 3)){
          ifelse(length(x) == 5, return(x[5]*call(x[c(1, 3)], KC) + (1 - x[5])*call(x[c(2, 4)], KC) ),
                 return(x[7]*call(x[c(1, 4)], KC) + x[8]*call(x[c(2, 5)], KC) + (1 - sum(x[7:8]))*call(x[c(3, 6)], KC)))
        } else{ ifelse(length(x) == 7, return(x[5]*call(x[c(1, 3)], KC) + (1 - x[5])*call(x[c(2, 4)], KC) ),
                       return(x[7]*call(x[c(1, 4)], KC) + x[8]*call(x[c(2, 5)], KC) + (1 - sum(x[7:8]))*call(x[c(3, 6)], KC))) }
      }

      esp <- function(x){exp(x[1] + (x[2]^2/2))}

      esp_mix <- function(x){
        if(cot %in%c(1, 3)){
          ifelse(length(x) == 5, x[5]*esp(x[c(1, 3)]) + (1 - x[5])*esp(x[c(2, 4)]),
                 x[7]*esp(x[c(1, 4)]) + x[8]*esp(x[c(2, 5)]) + (1 - sum(x[7:8]))*esp(x[c(3, 6)]))
        } else{ ifelse(length(x) == 7, x[5]*esp(x[c(1, 3)]) + (1 - x[5])*esp(x[c(2, 4)]),
                       x[7]*esp(x[c(1, 4)]) + x[8]*esp(x[c(2, 5)]) + (1 - sum(x[7:8]))*esp(x[c(3, 6)])) }
      }
      put <- function(x, KP){
        d1_C <- (x[1] + x[2]^2 - log(KP))/x[2]
        d2_C <- d1_C - x[2]
        if(cot %in%c(1, 2)){put <- exp(-r*T)*( -exp(x[1] + (x[2]^2/2))*pnorm(-d1_C) + KP*pnorm(-d2_C))
        } else(put <- -exp(x[1] + (x[2]^2/2))*pnorm(-d1_C) + KP*pnorm(-d2_C))
      }

      put_mix <- function(x, KP){
        if(cot %in%c(1, 3)){
          ifelse(length(x) == 5, return(x[5]*put(x[c(1, 3)], KP) + (1 - x[5])*put(x[c(2, 4)], KP) ),
                 return(x[7]*put(x[c(1, 4)], KP) + x[8]*put(x[c(2, 5)], KP) + (1 - sum(x[7:8]))*put(x[c(3, 6)], KP)))
        } else {  ifelse(length(x) == 7, return(x[5]*put(x[c(1, 3)], KP) + (1 - x[5])*put(x[c(2, 4)], KP) ),
                         return(x[7]*put(x[c(1, 4)], KP) + x[8]*put(x[c(2, 5)], KP) + (1 - sum(x[7:8]))*put(x[c(3, 6)], KP)))}
      }

      if(nb_log == 2) {PR <- seq(0.01, 0.49, 0.01)
      } else {PR <- seq(0.01, 1, 0.01)
      PR <- expand.grid(c(rep(list(PR), 2)))
      PR <- PR[rowSums(PR) < 1, ]}

      suppressWarnings({

        if(cot %in%c(1, 3)){
          model_prices <- function(x){
            return(list(model_call_price = call_mix(x, KC), model_put_price = put_mix(x, KP)))}
        } else {model_prices <- function(x){
          C_INF <- pmax(esp_mix(x) - KC, call_mix(x, KC))
          C_SUP <- exp(r*T)*call_mix(x, KC)
          P_INF <- pmax(KP - esp_mix(x), put_mix(x, KP))
          P_SUP <- exp(r*T)*put_mix(x, KP)
          itm_fwd_call <- as.numeric(KC <= esp_mix(x))
          itm_fwd_put <- as.numeric(KP >= esp_mix(x))
          w_call <- itm_fwd_call*first(tail(x, 2)) + (1 - itm_fwd_call)*last(x)
          w_put <- itm_fwd_put*first(tail(x, 2)) + (1 - itm_fwd_put)*last(x)
          CALL <- w_call*C_INF + (1 - w_call)*C_SUP
          PUT <- w_put*P_INF + (1 - w_put)*P_SUP
          return(list(model_call_price = CALL, model_put_price = PUT))}
        }

        MSE_mix <- function(x){
          MSE_mix <- sum((C - model_prices(x)$model_call_price)^2, na.rm = T) +
            sum((P - model_prices(x)$model_put_price)^2, na.rm = T) + (FWD - esp_mix(x))^2
          return(MSE_mix)}
      })

      objective <- function(x){
        if(cot %in% c(1, 3)){
          ifelse( length(PR) !=2,
                  MSE_mix( c(x[1:4], PR[i])),
                  MSE_mix( c(x[1:6], PR[i, 1], PR[i, 2] )))
        } else{ ifelse( length(PR) !=2,
                        MSE_mix( c(x[1:4], PR[i], x[6:7])),
                        MSE_mix( c(x[1:6], PR[i, 1], PR[i, 2], x[9:10] ))) }
      }

      C <- call_prices
      P <- put_prices
      KC <- call_strikes
      KP <- put_strikes
      T <- bond_fut$option_term
      FWD <- bond_fut$fut_price

      m1 <- m2 <- m3 <- s1 <- s2 <- s3 <- SCE <- NA

      if(nb_log == 2){PARA <- as.matrix(data.frame(m1, m2, s1, s2, pi1 = PR, w1 = 0.5, w2 = 0.5, SCE))
      } else {PARA <- as.matrix(data.frame(m1, m2, m3, s1, s2, s3, pi1 = PR[, 1], pi2 = PR[, 2],
                                           w1 = 0.5, w2 = 0.5, SCE))}

      start <- c(rep(c(log(FWD), 0.1), each = nb_log), rep(0.5, 2) )

      if(FWD != 1){
        lower <- c(rep(c( (sign(1 - FWD)*0.5 + 1)*log(FWD), 1e-6), each = nb_log), rep(0, 2) )
        upper <- c(rep(c( (sign(FWD - 1)*0.5 + 1)*log(FWD), 0.8), each = nb_log), rep(1, 2) )
      } else {
        lower <- c(rep(c( -1, 1e-6), each = nb_log), rep(0, 2) )
        upper <- c(rep(c( 1, 0.8), each = nb_log), rep(1, 2) ) }

      suppressWarnings({
        for (i in 1:length(PR)){
          if(cot %in%c(1,3)){
            sol <- nlminb(start = start[1:(length(start) - 2)], objective = objective,
                          lower = lower[1:(length(lower) - 2)], upper = upper[1:(length(upper) - 2)],
                          control = list(iter.max = 500))
            PARA[i, grep( paste( c("m", "s"), collapse = "|"), colnames(PARA))] <- sol$par[1:(2*nb_log)]
          } else{
            sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, control = list(iter.max = 500))
            PARA[i, grep( paste( c("m", "s", "w"), collapse = "|"), colnames(PARA))] <- sol$par[1:(2*nb_log + 2)] }
          PARA[i, "SCE"] <- sol$objective

        }
      })


      if(length(which(PARA[, ncol(PARA)]!="Inf")) != 0){
        PARA <- PARA
      } else{

        if(nb_log == 2){
          if(cot %in%c(1,3)){
            PARA <- data.frame(m1 = start[1], m2 = start[2], s1 = start[3], s2 = start[4],
                               pi1 = 0.5, SCE = NA)
          } else {PARA <- data.frame(m1 = start[1], m2 = start[2], s1 = start[3], s2 = start[4],
                                     pi1 = 0.5, w1 = 0.5, w2 = 0.5, SCE = NA)}
        }else{
          if(cot %in%c(1,3)){
            PARA <- data.frame(m1 = start[1], m2 = start[2], m3 = start[3], s1 = start[4],
                               s2 = start[5], s3 = start[6], pi1 = 0.33, pi2 = 0.33, SCE = NA)
          } else {data.frame(m1 = start[1], m2 = start[2], m3 = start[3], s1 = start[4],
                             s2 = start[5], s3 = start[6], pi1 = 0.33, pi2 = 0.33,
                             w1 = 0.5, w2 = 0.5, SCE = NA) }
        }
      }


      PARA <- PARA[ !is.na(PARA[, "m1"]), ]
      if(nrow(PARA) > 1){
        param <- PARA[which.min(PARA[, "SCE"]), -ncol(PARA)]
      } else {param <- as.matrix(PARA[, -ncol(PARA)], nrow = nrow(PARA))}

      param[param == 0] <- 1e-6

      L <- U <- rep(0, length(param))
      L[sign(param) == -1] <- 1.5*param[sign(param) == -1]
      L[sign(param) == 1] <- 1e-2*param[sign(param) == 1]
      L[(length(L) - 1):length(L)] <- 0
      if(cot%in%c(1, 3)){
        L <- L[1: (length(L)-2)]
      } else{L <- L }

      U[sign(param) == -1] <- 1e-2*param[sign(param) == -1]
      U[sign(param) == 1] <- 1.5*param[sign(param) == 1]
      U[(length(U) - 1):length(U)] <- 1
      if(cot%in%c(1, 3)){
        U <- U[1: (length(U)-2)]
      } else{U <- U }

      CI <- c(L, -U)
      UI <- rbind(diag(length(L)), -diag(length(L)))

      if(cot%in%c(1, 3)){
        param <- param[1: (length(param)-2)]
      } else{param <- param }

      param <- matrix(param)

      suppressWarnings({
        solu <- constrOptim(param, MSE_mix, NULL, ui = UI, ci = CI, mu = 1e-05, method = "Nelder-Mead")
      })

      model_p <- model_prices(c(solu$par))

      if(cot%in%c(1,3) ){ params <- solu$par
      } else {params <- solu$par[1: (length(solu$par) - 2)]}

      names(params) <- colnames(PARA)[1:(3*nb_log - 1)]

      range_px <- range(c(KP, KC))
      PX <- Reduce(seq, 1e3*range_px)*1e-3

      sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
      PDF <- function(x, y){
        ifelse(length(params) == 5,
               return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
               return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1 - sum(x[7:8])), y))) }

      sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }
      CDF <- function(x, y){
        ifelse(length(params) == 5,
               return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1 - x[5]), y) ),
               return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1 - sum( x[7:8])), y)) ) }

      DNR <- PDF(params, PX)

      if(sum(rollmean(PDF(params, PX), 2)*diff(PX), na.rm = T) < 1){

        x_axis <- 1e-2
        PX_2 <- PX
        range_px_2 <- range_px
        ratio <- 1
        while (ratio > 1e-5){
          integral <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
          range_px_2 <- c(1 - x_axis, 1)*range_px_2
          PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3
          integral_2 <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
          ratio <- integral_2 - integral}

        ratio <- 1
        while (ratio > 1e-5){
          integral <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
          range_px_2 <- c(1, 1 + x_axis)*range_px_2
          PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3
          integral_2 <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
          ratio <- integral_2 - integral}

        while (sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T) < 0.9991){
          range_px_2 <- c(1 - x_axis, 1 + x_axis)*range_px_2
          PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3}

        extension <- diff(range(PX_2))/diff(range(PX))

        if(extension <= 10){
          DNR_2 <- PDF(params, PX_2)

          NCDF <- CDF(params, PX_2)

          if(DNR_2[1] < DNR_2[2] &
             DNR_2[length(DNR_2) - 1] > DNR_2[length(DNR_2)] &
             min(DNR_2)%in%DNR_2[c(1, length(DNR_2))]){

            E <- sum(rollmean(PX_2*DNR_2, 2)*diff(PX_2))
            moments <- function(x){ return(sum(rollmean(DNR_2*(PX_2 - E)^x , 2)*diff(PX_2)))}
            SD <- sqrt(moments(2))
            SK <- moments(3)/SD^3
            KU <- moments(4)/SD^4
            moments <- c(mean = E, stddev = SD, skewness = SK, kurtosis = KU)

            df <- data.frame(price = PX_2, density = DNR_2)
            cdf <- data.frame(price = PX_2, cdf = NCDF)

            thres <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.995, 0.999)

            if(length(which(cdf$cdf > last(thres))) > 0 & length(which(cdf$cdf < first(thres))) > 0){
              quantiles <- list()
              for (j in 1:length(thres)){
                quantiles[[j]] <- mean(df$price[c(min(which(cdf$cdf > thres[j] - 1e-3)), max(which(cdf$cdf < thres[j] + 1e-3)))])}

              qt <- data.frame(quantiles) %>% rename_with(~paste0("q", 100*thres))

              mode <- PX_2[which.max(DNR_2)]

              graph <- PX_2 >= qt$q0.1 & PX_2 <= qt$q99.9
              PX_graph <- PX_2[graph]
              DNR_graph <- DNR_2[graph]
              NCDF_graph <- cdf$cdf[graph]

              df_graph <- data.frame(price = PX_graph, density = DNR_graph)
              cdf_graph <- data.frame(price = PX_graph, cdf = NCDF_graph)

              raw_bond_term <- as.numeric(bond_fut$ctd_matu - bond_fut$fut_matu)/365
              if(raw_bond_term <  2.5){ bond_term <- round(raw_bond_term/2)*2
              } else if (raw_bond_term < 5.5) { bond_term <- round(raw_bond_term/5)*5
              } else {bond_term <- round(raw_bond_term/10)*10}

              pdf <- ggplot() + geom_line(data = df_graph, aes(x = price, y = density)) +
                labs(x = paste0("futures price (", bond_fut$currency, ") of maturity ", bond_fut$fut_matu),
                     y = "probability density") + theme_bw() +
                theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                labs(title = paste0(bond_term, "-y ", bond_fut$nationality, " bond futures price on ",
                                    bond_fut$option_matu, " as of ", bond_fut$start_date),
                     subtitle = paste0("Risk Neutral Probability Density for a mixture of ", nb_log, " lognormals"))

              ncdf <- ggplot() + geom_line(data = cdf_graph, aes(x = price, y = cdf)) +
                labs(x = paste0("futures price (", bond_fut$currency, ") of maturity ", bond_fut$fut_matu),
                     y = "probability density") + theme_bw() +
                theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                labs(title = paste0(bond_term, "-y ", bond_fut$nationality, " bond futures price on ",
                                    bond_fut$option_matu, " as of ", bond_fut$start_date),
                     subtitle = paste0("Risk Neutral Cumulative Probability for a mixture of ", nb_log, " lognormals"))

              bond_future_price <- list(params, df$price, df$density, cdf$cdf, solu$convergence, moments, qt, mode, model_p, pdf, ncdf)
              names(bond_future_price) <- c("params", "domain", "rnd", "cdf", "CV", "moments", "quantiles", "mode", "model_prices", "rnd_plot", "cdf_plot")

              return(bond_future_price)

            } else {message(paste0("A mixture of ", nb_log, " lognormal distributions is not consistent with this options's style"))}
          } else {message(paste0("A mixture of ", nb_log, " lognormal distributions is not consistent with this options's style"))}
        } else {message("impossible to retrieve a density")}
      } else {message("impossible to retrieve a density")}
    } else {message("input dates are not consistent")}
  } else {message("inputs do not have the required length")}
}
