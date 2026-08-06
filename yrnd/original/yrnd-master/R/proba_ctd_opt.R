#' proba_ctd_opt
#'
#' @param call_prices a vector of call prices on a bond futures, in numeric format
#' @param call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param put_prices a vector of put prices on a bond futures, in numeric format
#' @param put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param nb_log a number for the number of component densities in the lognormal mixture to model the bond futures' price, either 2 or 3, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param r_2 a number for the riskfree spot rate whose maturity is equal to the futures' maturity, in numeric format
#' @param day_count_conv a number for the day count convention, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot a number for the options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param ctd_matu a date for the maturity date of the current Cheapest-to-Deliver Bond of the bond futures, in Date format
#' @param fut_price a number for the bond futures' price on calibration date, in numeric format
#' @param fut_matu a date for the bond futures' maturity date, in Date format
#' @param option_matu a date for the maturity date of the options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param bond_ISIN a vector of ISIN codes for the bonds in the delivery basket of the bond futures, in character format
#' @param bond_coupon a vector of the corresponding coupon rates for the bonds in the delivery basket, in numeric format
#' @param bond_cp_f a vector of the corresponding frequencies of coupon payment, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_matu a vector of the corresponding maturity dates for the bonds in the basket of deliverable bonds, in Date format
#' @param bond_nomi a single number for the nominal of the bonds (100 by default) in numeric format
#' @param bond_conv_factor a vector of the corresponding conversion factors for the bonds in the delivery basket, in numeric format
#' @param bond_ytm a vector of the corresponding yield to maturities at observation date for the bonds in the delivery basket, in numeric format
#' @param sett a number for the number of days between the ex-coupon date and the coupon payment date of the bonds in the delivery basket, in numeric format
#'
#' @returns for the bonds in the delivery basket, their ISIN in character format and their probability of being the CtD bond at options' maturity, in numeric format
#' @export
#'
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm
#' @importFrom utils head tail
#' @import dplyr
#' @import lubridate
#' @import zoo
#' @import ggplot2
#' @import tvm
#'
#' @examples
#' \donttest{
#' proba_ctd_opt( c(24.10, 23.10, 22.12, 21.12, 20.12, 19.14, 18.14, 17.16, 16.18, 15.20,
#' 14.22, 13.24, 12.28, 11.32, 10.36, 9.44, 8.50, 7.60, 6.72, 5.86, 5.04, 4.28,
#' 3.56, 2.88, 2.30, 1.78,  1.36, 1.02, 0.76, 0.56, 0.42, 0.30, 0.22, 0.18, 0.14,
#' 0.10,  0.08, 0.06, 0.06, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
#' seq(85, 131),
#' c(0.02,  0.02, 0.02, 0.02, 0.02, 0.04, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.18, 0.22,
#' 0.26, 0.34,  0.40, 0.50, 0.62, 0.76, 0.94, 1.18, 1.46, 1.78, 2.20, 2.68, 3.26, 3.92,
#' 4.66,  5.46, 6.32, 7.20, 8.12, 9.08, 10.04, 11.00, 11.98, 12.96, 13.96, 14.94, 15.94,
#' 16.92, 17.92, 18.92, 19.92, 20.92, 21.92),
#' seq(85, 131),
#' 2,
#' 0.0187,
#' 0.019,
#' 1,
#' 3,
#' as.Date("2054-08-15"),
#' 109.1,
#' as.Date("2026-09-08"),
#' as.Date("2026-08-21"),
#' as.Date("2026-05-28"),
#' c("DE0001102572", "DE0001102614", "DE0001030757", "DE000BU2D004", "DE000BU2D012"),
#' c(0.000, 0.018, 0.018, 0.025, 0.029),
#' rep(1, 5),
#' as.Date(c("2052-08-15", "2053-08-15", "2053-08-15", "2054-08-15", "2056-08-15")),
#' 100,
#' c(0.361698, 0.641260, 0.641260, 0.750372, 0.809987),
#' c(0.03500, 0.03507, 0.03492, 0.03510, 0.03514),
#' 2)
#' }
#'

proba_ctd_opt <- function(call_prices, call_strikes, put_prices, put_strikes, nb_log, r, r_2, day_count_conv, cot,
                      ctd_matu, fut_price, fut_matu, option_matu, start_date, bond_ISIN, bond_coupon,
                      bond_cp_f, bond_matu, bond_nomi = 100, bond_conv_factor, bond_ytm, sett){

  if(length(nb_log) == 1 & length(r) == 1 & length(day_count_conv) == 1 & length(cot) == 1 & length(ctd_matu) == 1 &
     length(fut_price) == 1 & length(fut_matu) == 1 & length(option_matu) == 1 & length(start_date) == 1 &
     length(call_prices) > 1 & length(call_strikes) > 1 & length(put_prices) > 1 & length(put_strikes) > 1 &
     length(bond_ISIN) > 1 & identical(length(bond_ISIN), length(bond_coupon), length(bond_cp_f), length(bond_matu),
                                       length(bond_conv_factor), length(bond_ytm))){

    if(start_date < option_matu & option_matu <= fut_matu &
       length(which(as.Date(fut_matu) - as.Date(bond_matu) > 0)) == 0 ){

      deliverables <- data.frame(bond_ISIN, bond_coupon, bond_matu, bond_conv_factor, bond_ytm, bond_cp_f) %>%
        rename_with(~c("ISIN", "coupon", "matu", "conv_factor", "ytm", "cp_freq")) %>% mutate_at("matu", as.Date) %>%
        mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"), "-", format(matu, "%m-%d"))))

      Nomi <- bond_nomi
      sett <- sett

      deliv_bonds <- bind_rows(deliverables[deliverables$cp_freq == 1,] %>%
                                 mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu < ., . %m-% years(1), .))) %>%
                                 mutate(curr_cp_dt = prev_cp_dt %m+% years(1), next_cp_dt = curr_cp_dt %m+% years(1)),
                               deliverables[deliverables$cp_freq == 2, ] %>%
                                 mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu - . < - months(6), . %m-% years(1),
                                                                         ifelse(option_matu - . < 0, . %m-% months(6), .)))) %>%
                                 mutate(curr_cp_dt = prev_cp_dt %m+% months(6), next_cp_dt = curr_cp_dt %m+% months(6)))

      terms <- data.frame(fut_matu, option_matu) %>% mutate_all(as.Date)

      if(day_count_conv == 1){
        terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                    as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")),
                                  res_term = as.numeric(fut_matu - option_matu)/
                                    as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
      } else if(day_count_conv == 2){
        terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/360,
                                  res_term = as.numeric(fut_matu - option_matu)/360)
      } else if(day_count_conv == 3){
        terms <- terms %>% mutate(option_term =as.numeric(option_matu - start_date)/365,
                                  res_term = as.numeric(fut_matu - option_matu)/365)
      } else {
        terms <- terms %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                  stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                  plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                     ceiling_date(start_date, "months")))/30),
                                  option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                  stub_1_res = max(0, 30 - as.numeric(format(option_matu + sett, "%d"))),
                                  stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                  plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                        ceiling_date(option_matu + sett, "months") )/30),
                                  res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

      terms <- terms %>% mutate(res_term_2 = 0)

      rate_table <- data.frame(term = terms$option_term + c(0, terms$res_term), rates = c(r, r_2)) %>%
        mutate(d_fact = term*rates)

      fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

      bef <- which(terms$fut_matu < deliv_bonds$curr_cp_dt)
      aft <- which(terms$fut_matu >= deliv_bonds$curr_cp_dt)

      if(day_count_conv == 1){
        deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett)/
                                                                 as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") )),
                                 deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett)/
                                                                 as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") ),
                                                               acc_matu = Nomi*coupon/cp_freq*res_term_2 ))
      } else if(day_count_conv == 2) {
        deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett)/360),
                                 deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett)/360,
                                                               acc_matu = Nomi*coupon/cp_freq*res_term_2))
      } else if(day_count_conv == 3){
        deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett)/365),
                                 deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett)/365,
                                                               acc_matu = Nomi*coupon/cp_freq*res_term_2))
      } else{
        deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(prev_cp_dt + sett, "%d"))),
                                                               stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                               plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                 ceiling_date(prev_cp_dt + sett, "months") )/30),
                                                               acc_matu = Nomi*coupon/cp_freq/360*(stub_2 + stub_1 + max(0, plain_months)*30)),
                                 deliv_bonds[aft, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(curr_cp_dt + sett, "%d"))),
                                                               stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                               plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                 ceiling_date(curr_cp_dt + sett, "months") )/30),
                                                               res_term_2 = (stub_1 + stub_2 + max(0, plain_months)*30)/360,
                                                               acc_matu = Nomi*coupon/cp_freq*res_term_2))
      }

      deliv_bonds <- deliv_bonds[match(deliverables$ISIN, deliv_bonds$ISIN), ]

      rate_table_bond <- fwd_2 <- list()
      for (i in aft){
        rate_table_bond[[i]] <- rate_table %>%
          add_row(term = terms$res_term + terms$option_term - deliv_bonds$res_term_2[i])
        rate_table_bond[[i]]$rates[3] <- approx(rate_table_bond[[i]]$term[1:2],
                                                rate_table_bond[[i]]$rates[1:2],
                                                xout = rate_table_bond[[i]]$term[3],
                                                method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = F)$y
        rate_table_bond[[i]]$d_fact[3] <- rate_table_bond[[i]]$term[3]*rate_table_bond[[i]]$rates[3]
        fwd_2[[i]] <- diff(rate_table_bond[[i]]$d_fact[-1])/diff(rate_table_bond[[i]]$term[-1])
      }

      true_cp_dt <- stub <- cp_dt_2 <- list()

      for (i in 1:nrow(deliv_bonds)){

        if(deliv_bonds$cp_freq[i] == 1){
          true_cp_dt[[i]] <- seq(from = deliv_bonds$curr_cp_dt[i], to = deliv_bonds$matu[i], by = "year")
        } else { true_cp_dt[[i]] <- seq(from = deliv_bonds$curr_cp_dt[i], to = deliv_bonds$matu[i], by = "quarter")
        true_cp_dt[[i]] <- true_cp_dt[[i]][seq(1, length(true_cp_dt[[i]]), by = 2)] }

        true_cp_dt[[i]] <- c(head(true_cp_dt[[i]], -1) + sett, tail(true_cp_dt[[i]], 1))
        true_cp_dt[[i]] <- c(terms$option_matu, true_cp_dt[[i]])

        if(day_count_conv == 1){
          if(deliv_bonds$cp_freq[i] == 1){ stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
            as.numeric(deliv_bonds$curr_cp_dt[i] - deliv_bonds$prev_cp_dt[i])
          } else {stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
            as.numeric(ceiling_date(deliv_bonds$prev_cp_dt[i] + sett, "year") -
                         floor_date(deliv_bonds$prev_cp_dt[i] + sett, "year")) }
          cp_dt_2[[i]] <- stub[[i]] + c(0, seq(1, length(true_cp_dt[[i]]) - 2)/deliv_bonds$cp_freq[i])
        } else if(day_count_conv == 2){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/360
        } else if(day_count_conv == 3){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/365
        } else {stub[[i]] <- (max(0, 30 - as.numeric(format(true_cp_dt[[i]][1], "%d"))) +
                                as.numeric(round((floor_date(true_cp_dt[[i]][2], "months") -
                                                    ceiling_date(true_cp_dt[[i]][1], "months"))/30))*30 +
                                min(30, as.numeric(format(true_cp_dt[[i]][2], "%d"))))/360
        cp_dt_2[[i]] <- stub[[i]] + c(0, seq(1, length(true_cp_dt[[i]]) - 2)/deliv_bonds$cp_freq[i])}

        cp_dt_2[[i]] <- list(cp_dt_2[[i]])
      }

      cf_matu <- Nomi*(1 + deliv_bonds$coupon/deliv_bonds$cp_freq)
      cf_other <- split(rep(deliv_bonds$coupon/deliv_bonds$cp_freq*Nomi,
                            sapply(cp_dt_2, lengths) - 1),
                        rep(seq_along(cp_dt_2), sapply(cp_dt_2, lengths) - 1))


      tri <- function(x){
        if(terms$fut_matu < deliv_bonds$curr_cp_dt[i]){
          tri <- mapply(xirr, cf = mapply(c, -(x*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] )*exp( -fwd_1*terms$res_term),
                                          cf_other[i], cf_matu[i], SIMPLIFY = F),
                        tau = mapply(c, 0, mapply(unlist, cp_dt_2[i], SIMPLIFY = F), SIMPLIFY = F))
        } else{tri <- mapply(xirr, cf = mapply(c, -(x*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] +
                                                      deliv_bonds$coupon[i]/deliv_bonds$cp_freq[i]*Nomi*exp(fwd_2[[i]]*deliv_bonds$res_term_2[i]) )*exp( -fwd_1*terms$res_term),
                                               cf_other[i], cf_matu[i], SIMPLIFY = F),
                             tau = mapply(c, 0, mapply(unlist, cp_dt_2[i], SIMPLIFY = F), SIMPLIFY = F) )}
      }

      bond_fut <- bond_future_price(call_prices, call_strikes, put_prices, put_strikes, nb_log, r,
                                    day_count_conv, cot, ctd_matu, fut_price, fut_matu, option_matu, start_date)

      if(length(bond_fut) > 0 ){

        ytm <- list()
        for (i in 1:nrow(deliv_bonds)){
          ytm[[i]] <- tri(bond_fut$domain)}

        ytm <- do.call(cbind, ytm) %>% data.frame %>% rename_with(~c(deliv_bonds$ISIN))

        ctd <- list()
        for (i in 1:nrow(deliv_bonds)){
          ctd[[i]] <- replicate(nrow(deliv_bonds), ytm[, i] - deliv_bonds$ytm[i]) +
            t(replicate(nrow(ytm), deliv_bonds$ytm)) }

        dirty <- function(x){
          dcf <- mapply("/", list(c(unlist(cf_other[i]), cf_matu[i])),
                        mapply("^", 1 + x, list(unlist(cp_dt_2[[i]])),  SIMPLIFY = F), SIMPLIFY = F)
          dirty <- unlist(lapply(dcf, sum))}

        net_basis <- ctd_conf <- ctd_conf_2 <- list()

        for (k in 1:length(ctd)){
          net_basis[[k]] <- list()
          for (i in which(terms$fut_matu < deliv_bonds$curr_cp_dt)){
            net_basis[[k]][[i]] <- dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
              (bond_fut$domain*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i]  )}
          for (i in which(terms$fut_matu >= deliv_bonds$curr_cp_dt)){
            net_basis[[k]][[i]] <-  dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
              (bond_fut$domain*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] +
                 deliv_bonds$coupon[i]/deliv_bonds$cp_freq[i]*Nomi*exp(fwd_2[[i]]*deliv_bonds$res_term_2[i]) ) }
          net_basis[[k]] <- do.call(cbind, net_basis[[k]])
          ctd_conf_2[[k]] <- ctd_conf[[k]] <- apply(net_basis[[k]], 1, which.min)
          ctd_conf_2[[k]][ctd_conf[[k]] != k] <- 0
        }

        ctd_conf_2 <- do.call(cbind, ctd_conf_2)

        for (i in 1:length(ctd)){
          ctd_conf_2[ctd_conf_2[, i] == i, -i] <- 0}

        prb <- rowSums(ctd_conf_2)

        if(length(which(prb == 0)) > 0){
          prb[prb == 0] <- ctd_conf[[length(ctd)]][which(prb == 0)]
        }

        ctd_pot <- unique(prb)

        prob <- list()
        for (i in 1:length(ctd_pot)){
          prob[[i]] <- sum(bond_fut$rnd[prb == ctd_pot[i]])*first(diff(bond_fut$domain))}

        prob <- round(unlist(prob), 3)

        probas <- data.frame(ISIN = c(deliv_bonds$ISIN[ctd_pot], deliv_bonds$ISIN[-ctd_pot]),
                             proba_ctd_matu = c(prob, rep(0, length(deliv_bonds$ISIN[-ctd_pot])))) %>%
          arrange(desc(proba_ctd_matu))

        return(probas)
      } else {message("impossible to retrieve the probabilities to be CtD")}
    } else {message("input dates are not consistent")}
  } else {message("inputs do not have the required length")}
}
