#' bond_yield_spread
#'
#' @param call_prices_1 a vector of call prices on the first bond futures, in numeric format
#' @param call_strikes_1 a vector of call strikes attached to the call prices, in numeric format
#' @param put_prices_1 a vector of put prices on the first bond futures, in numeric format
#' @param put_strikes_1 a vector of put strikes attached to the put prices, in numeric format
#' @param call_prices_2 a vector of call prices on the second bond futures, in numeric format
#' @param call_strikes_2 a vector of call strikes attached to the call prices, in numeric format
#' @param put_prices_2 a vector of put prices on the second bond futures, in numeric format
#' @param put_strikes_2 a vector of put strikes attached to the put prices, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param r_2 a number for the riskfree spot rate whose maturity is equal to the futures contracts' maturity, in numeric format
#' @param day_count_conv_1 a number for the day count convention of the first bond futures, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param day_count_conv_2 a number for the day count convention of the second bond futures, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot_1 a number for style of the options on the first bond futures, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param cot_2 a number for style of the options on the second bond futures, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param ctd_matu_1 a date for the maturity date of the current Cheapest-to-Deliver Bond of the first bond futures, in Date format
#' @param ctd_matu_2 a date for the maturity date of the current Cheapest-to-Deliver Bond of the second bond futures, in Date format
#' @param fut_price_1 a number for the first bond futures' price on calibration date, in numeric format
#' @param fut_price_2 a number for the second bond futures' price on calibration date, in numeric format
#' @param fut_matu a date for the bond futures' maturity date, in Date format
#' @param option_matu a date for the maturity date of the options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param bond_ISIN_1 a vector of ISIN codes for the bonds in the delivery basket of the first bond futures, in character format
#' @param bond_coupon_1 a vector of the corresponding coupon rates for the bonds in the delivery basket, in numeric format
#' @param bond_cp_f_1 a vector of the corresponding frequencies of coupon payment, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_matu_1 a vector of the corresponding maturity dates for the bonds in the basket of deliverable bonds, in Date format
#' @param bond_conv_factor_1 a vector of the corresponding conversion factors for the bonds in the delivery basket, in numeric format
#' @param bond_ytm_1 a vector of the corresponding yield to maturities at observation date for the bonds in the delivery basket, in numeric format
#' @param sett_1 a number for the number of days between the ex-coupon date and the coupon payment date of the bonds in the delivery basket, in numeric format
#' @param bond_ISIN_2 a vector of ISIN codes for the bonds in the delivery basket of the second bond futures, in character format
#' @param bond_coupon_2 a vector of the corresponding coupon rates for the bonds in the delivery basket, in numeric format
#' @param bond_cp_f_2 a vector of the corresponding frequencies of coupon payment, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_matu_2 a vector of the corresponding maturity dates for the bonds in the basket of deliverable bonds, in Date format
#' @param bond_conv_factor_2 a vector of the corresponding conversion factors for the bonds in the delivery basket, in numeric format
#' @param bond_ytm_2 a vector of the corresponding yield to maturities at observation date for the bonds in the delivery basket, in numeric format
#' @param sett_2 a number for the number of days between the ex-coupon date and the coupon payment date of the bonds in the delivery basket, in numeric format
#' @param bond_nomi a single number for the nominal of the bonds (100 by default) in numeric format
#' @param corr a number for the correlation coefficient between the returns on the two bond futures prices with identical maturity, in numeric format
#' @param ctry_1 a character for the nationality of the issuer of the bond in the first bond futures, in character format (NA by default)
#' @param ctry_2 a character for the nationality of the issuer of the bond in the second bond futures, in character format (NA by default)
#'
#' @returns 10,000 realizations of the bond yield spread between two different issuers or two issues of the same issuer with different maturities, using options on two bond futures with the same maturity, and a density plot of the spread, in bps
#' @export
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm qlnorm sd
#' @importFrom utils head tail
#' @importFrom MASS mvrnorm
#' @import dplyr
#' @import lubridate
#' @import tvm
#' @import zoo
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' bond_yield_spread(
#' c(12.64, 12.14, 11.65, 11.15, 10.65, 10.15,  9.65,  9.16, 8.66,
#' 8.16, 7.67,  7.17, 6.68, 6.19, 5.70, 5.21, 4.73, 4.25, 3.78, 3.32,
#' 2.87,  2.45, 2.05, 1.67, 1.33, 1.03, 0.77, 0.57, 0.41, 0.29,
#' 0.21,  0.15, 0.10, 0.07, 0.06, 0.04, 0.04, 0.03, 0.02, 0.02, 0.02,
#' 0.02,  0.01, 0.01, 0.01, 0.01, 0.01),
#' c(seq(114, 136.5, 0.5), 137.5),
#' c(0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.04, 0.04,
#' 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.15, 0.19, 0.24, 0.32, 0.42, 0.54,
#' 0.70, 0.90, 1.14, 1.44, 1.78, 2.16, 2.58, 3.02, 3.47, 3.94, 4.43, 4.91,
#' 5.41, 5.90, 6.39, 6.89, 7.39, 7.89, 8.38, 8.88, 9.38, 9.88, 10.88),
#' c(seq(114, 136.5, 0.5), 137.5),
#' c(13.46, 12.96, 12.46, 11.96, 11.47, 10.97, 10.47, 9.98, 9.48, 8.99, 8.49,
#' 8.00,  7.51, 7.02, 6.53, 6.04, 5.56, 5.08, 4.61, 4.15, 3.69, 3.25, 2.82,
#' 2.41,  2.02, 1.67, 1.34, 1.05, 0.80, 0.59, 0.43, 0.30, 0.21, 0.14, 0.09,
#' 0.06, 0.04,  0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
#' seq(105.5, 129, 0.5),
#' c(0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06,
#' 0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.17, 0.21, 0.25, 0.31, 0.38, 0.47,
#' 0.58, 0.73, 0.90, 1.11, 1.36, 1.65, 1.99, 2.36, 2.77, 3.20, 3.65, 4.12,
#' 4.60, 5.08, 5.57, 6.07, 6.57, 7.06, 7.56, 8.06, 8.56, 9.06, 9.56, 10.06),
#' seq(105.5, 129, 0.5),
#' 0.0187,
#' 0.019,
#' 1,
#' 1,
#' 3,
#' 3,
#' as.Date("2035-08-15"),
#' as.Date("2035-08-01"),
#' 126.66,
#' 119.04,
#' as.Date("2026-09-08"),
#' as.Date("2026-08-21"),
#' as.Date("2026-06-17"),
#' c("DE000BU2Z056", "DE000BU2Z064"),
#' c(0.026, 0.029),
#' rep(1, 2),
#' as.Date(c("2035-08-15", "2036-02-15")),
#' c(0.770088, 0.781249),
#' c(0.02898, 0.02926),
#' 2,
#' c("IT0005631590", "IT0005648149", "IT0005676504", "IT0005402117",
#' "IT0005706285", "IT0005433195"),
#' c(0.0365, 0.0360, 0.0345, 0.0145, 0.0380, 0.0095),
#' rep(2, 6),
#' as.Date(c("2035-08-01", "2035-10-01", "2036-02-01", "2036-03-01",
#' "2036-07-01", "2037-03-01")),
#' c(0.845212, 0.839556, 0.824415, 0.679587, 0.844354, 0.616457),
#' c(0.03567, 0.03587, 0.03620, 0.03614, 0.03665, 0.03712),
#' 2,
#' 100,
#' 0.4,
#' "Germany",
#' "Italy")
#' }
#'

bond_yield_spread <- function(call_prices_1, call_strikes_1, put_prices_1, put_strikes_1,
                              call_prices_2, call_strikes_2, put_prices_2, put_strikes_2,
                              r, r_2, day_count_conv_1, day_count_conv_2, cot_1, cot_2,
                              ctd_matu_1, ctd_matu_2, fut_price_1, fut_price_2, fut_matu, option_matu, start_date,
                              bond_ISIN_1, bond_coupon_1, bond_cp_f_1, bond_matu_1, bond_conv_factor_1, bond_ytm_1,
                              sett_1, bond_ISIN_2, bond_coupon_2, bond_cp_f_2, bond_matu_2, bond_conv_factor_2, bond_ytm_2,
                              sett_2, bond_nomi = 100, corr, ctry_1 = NA, ctry_2 = NA){

  if(length(fut_matu) == 1 & length(option_matu) == 1 & length(start_date) == 1 & length(corr) == 1 &
     length(r) == 1 & length(r_2) == 1 & length(day_count_conv_1) == 1 &
     length(cot_1) == 1 & length(ctd_matu_1) == 1 & length(fut_price_1) == 1 & length(sett_1 == 1) &
     length(call_prices_1) > 1 & length(call_strikes_1) > 1 & length(put_prices_1) > 1 & length(put_strikes_1) > 1 &
     length(bond_ISIN_1) > 1 &
     identical(length(bond_ISIN_1), length(bond_coupon_1), length(bond_cp_f_1), length(bond_matu_1),
               length(bond_conv_factor_1), length(bond_ytm_1)) &
     length(day_count_conv_2) == 1 & length(cot_2) == 1 & length(ctd_matu_2) == 1 &
     length(fut_price_2) == 1 & length(sett_2 == 1) &
     length(call_prices_2) > 1 & length(call_strikes_2) > 1 & length(put_prices_2) > 1 & length(put_strikes_2) > 1 &
     length(bond_ISIN_2) > 1 &
     identical(length(bond_ISIN_2), length(bond_coupon_2), length(bond_cp_f_2), length(bond_matu_2),
               length(bond_conv_factor_2), length(bond_ytm_2))){

    if(as.Date(start_date) < as.Date(option_matu) & as.Date(option_matu) <= as.Date(fut_matu) &
       length(which(as.Date(fut_matu) - as.Date(c(bond_matu_1, bond_matu_2)) > 0)) == 0 ){

      rnd_1 <- bond_future_price(call_prices_1, call_strikes_1, put_prices_1, put_strikes_1, 2, r,
                                 day_count_conv_1, cot_1, ctd_matu_1, fut_price_1, fut_matu, option_matu, start_date)

      rnd_2 <- bond_future_price(call_prices_2, call_strikes_2, put_prices_2, put_strikes_2, 2, r,
                                 day_count_conv_2, cot_2, ctd_matu_2, fut_price_2, fut_matu, option_matu, start_date)

      if(length(rnd_1) > 0 & length(rnd_2) > 0){

        suppressMessages({

          params <- bind_cols(rnd_1$params, rnd_2$params) %>% t %>% data.frame %>%
            bind_cols(c("ctry_1", "ctry_2")) %>% bind_cols(c(fut_price_1, fut_price_2)) %>%
            rename_with(~c(names(rnd_1$params), "ctry", "F_0")) %>% relocate(c(ctry, F_0))

          corr_matrix <- matrix(c(1, rep(corr, nrow(params)) , 1), nrow = nrow(params),
                                dimnames = list(params$ctry, params$ctry) )

          simulate_mixture <- function(x, U) {
            U1 <- pmin(U/x[5], 1)
            U2 <- pmax((U - x[5])/(1 - x[5]), 0)
            ifelse(U < x[5], qlnorm(U1, meanlog = x[1], sdlog = x[3]),
                   qlnorm(U2, meanlog = x[2], sdlog = x[4]) )
          }

          simulate_correlated_returns <- function(n, params, corr_matrix) {
            Z <- mvrnorm(n, mu = rep(0, nrow(params)), Sigma = corr_matrix, tol = 1e-06, empirical = FALSE)
            U <- pnorm(Z)
            simR <- data.frame(matrix(nrow = n, ncol = nrow(params)))  %>% rename_with(~params$ctry)
            for (i in 1:nrow(params)) {
              p <- params %>% filter(ctry == params$ctry[i]) %>% dplyr::select(-c(ctry, F_0)) %>% unlist
              F_T <- simulate_mixture(p, U[ ,i])
              simR[, i] <- (F_T/params$F_0[i]) - 1 }
            simR
          }


          joint <- simulate_correlated_returns(10000, params, corr_matrix)
          joint <- 1 + joint
          joint <- joint*t(replicate(nrow(joint), params$F_0))

          tri <- function(x){
            if(terms$fut_matu < deliv_bonds$curr_cp_dt[i]){
              tri <- mapply(xirr, cf = mapply(c, -(x*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] )*exp( -fwd_1*terms$res_term),
                                              cf_other[i], cf_matu[i], SIMPLIFY = F),
                            tau = mapply(c, 0, mapply(unlist, cp_dt_2[i], SIMPLIFY = F), SIMPLIFY = F) )
            } else{tri <- mapply(xirr, cf = mapply(c, -(x*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] +
                                                          deliv_bonds$coupon[i]/deliv_bonds$cp_freq[i]*Nomi*exp(fwd_2[[i]]*deliv_bonds$res_term_2[i]) )*exp( -fwd_1*terms$res_term),
                                                   cf_other[i], cf_matu[i], SIMPLIFY = F),
                                 tau = mapply(c, 0, mapply(unlist, cp_dt_2[i], SIMPLIFY = F), SIMPLIFY = F) )}
          }

          dirty <- function(x){
            dcf <- mapply("/", list(c(unlist(cf_other[i]), cf_matu[i])),
                          mapply("^", 1 + x, list(unlist(cp_dt_2[[i]])), SIMPLIFY = F), SIMPLIFY = F)
            dirty <- unlist(lapply(dcf, sum))}


          deliverables_1 <- data.frame(bond_ISIN_1, bond_coupon_1, bond_matu_1, bond_conv_factor_1, bond_ytm_1, bond_cp_f_1) %>%
            rename_with(~c("ISIN", "coupon", "matu", "conv_factor", "ytm", "cp_freq")) %>% mutate_at("matu", as.Date) %>%
            mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"), "-", format(matu, "%m-%d"))))

          deliverables_2 <- data.frame(bond_ISIN_2, bond_coupon_2, bond_matu_2, bond_conv_factor_2, bond_ytm_2, bond_cp_f_2) %>%
            rename_with(~c("ISIN", "coupon", "matu", "conv_factor", "ytm", "cp_freq")) %>% mutate_at("matu", as.Date) %>%
            mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"), "-", format(matu, "%m-%d"))))

          Nomi <- bond_nomi
          sett_1 <- sett_1
          sett_2 <- sett_2

          terms <- data.frame(fut_matu, option_matu) %>% mutate_all(as.Date)


          deliv_bonds <- bind_rows(deliverables_1[deliverables_1$cp_freq == 1,] %>%
                                     mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu < ., . %m-% years(1), .))) %>%
                                     mutate(curr_cp_dt = prev_cp_dt %m+% years(1)),
                                   deliverables_1[deliverables_1$cp_freq == 2, ] %>%
                                     mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu - . < - months(6), . %m-% years(1),
                                                                             ifelse(option_matu - . < 0, . %m-% months(6), .)))) %>%
                                     mutate(curr_cp_dt = prev_cp_dt %m+% months(6)))

          if(day_count_conv_1 == 1){
            terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                        as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")),
                                      res_term = as.numeric(fut_matu - option_matu)/
                                        as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
          } else if(day_count_conv_1 == 2){
            terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/360,
                                      res_term = as.numeric(fut_matu - option_matu)/360)
          } else if(day_count_conv_1 == 3){
            terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/365,
                                      res_term = as.numeric(fut_matu - option_matu)/365)
          } else {
            terms <- terms %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                      stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                      plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                         ceiling_date(start_date, "months")))/30),
                                      option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                      stub_1_res = max(0, 30 - as.numeric(format(option_matu + sett_1, "%d"))),
                                      stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                      plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                            ceiling_date(option_matu + sett_1, "months") )/30),
                                      res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

          terms <- terms %>% mutate(res_term_2 = 0)

          rate_table <- data.frame(term = terms$option_term + c(0, terms$res_term), rates = c(r, r_2)) %>%
            mutate(d_fact = term*rates)

          fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

          bef <- which(terms$fut_matu < deliv_bonds$curr_cp_dt)
          aft <- which(terms$fut_matu >= deliv_bonds$curr_cp_dt)

          if(day_count_conv_1 == 1){
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_1)/
                                                                     as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") )),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_1)/
                                                                     as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") ),
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2 ))
          } else if(day_count_conv_1 == 2) {
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_1)/360),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_1)/360,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          } else if(day_count_conv_1 == 3){
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_1)/365),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_1)/365,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          } else{
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(prev_cp_dt + sett_1, "%d"))),
                                                                   stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                                   plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                     ceiling_date(prev_cp_dt + sett_1, "months") )/30),
                                                                   acc_matu = Nomi*coupon/cp_freq/360*(stub_2 + stub_1 + max(0, plain_months)*30)),
                                     deliv_bonds[aft, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(curr_cp_dt + sett_1, "%d"))),
                                                                   stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                                   plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                     ceiling_date(curr_cp_dt + sett_1, "months") )/30),
                                                                   res_term_2 = (stub_1 + stub_2 + max(0, plain_months)*30)/360,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          }

          deliv_bonds <- deliv_bonds[match(deliverables_1$ISIN, deliv_bonds$ISIN), ]

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

            true_cp_dt[[i]] <- c(head(true_cp_dt[[i]], -1) + sett_1, tail(true_cp_dt[[i]], 1))
            true_cp_dt[[i]] <- c(terms$option_matu, true_cp_dt[[i]])

            if(day_count_conv_1 == 1){
              if(deliv_bonds$cp_freq[i] == 1){ stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
                as.numeric(deliv_bonds$curr_cp_dt[i] - deliv_bonds$prev_cp_dt[i])
              } else {stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
                as.numeric(ceiling_date(deliv_bonds$prev_cp_dt[i] + sett_1, "year") -
                             floor_date(deliv_bonds$prev_cp_dt[i] + sett_1, "year")) }
              cp_dt_2[[i]] <- stub[[i]] + c(0, seq(1, length(true_cp_dt[[i]]) - 2)/deliv_bonds$cp_freq[i])
            } else if(day_count_conv_1 == 2){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/360
            } else if(day_count_conv_1 == 3){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/365
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

          ytm <- list()
          for (i in 1:nrow(deliv_bonds)){
            ytm[[i]] <- tri(joint$ctry_1)}

          ytm <- do.call(cbind, ytm) %>% data.frame %>% rename_with(~c(deliv_bonds$ISIN))


          ctd <- list()
          for (i in 1:nrow(deliv_bonds)){
            ctd[[i]] <- replicate(nrow(deliv_bonds), ytm[, i] - deliv_bonds$ytm[i]) +
              t(replicate(nrow(ytm), deliv_bonds$ytm)) }

          net_basis <- ctd_conf <- ctd_conf_2 <- list()

          for (k in 1:length(ctd)){
            net_basis[[k]] <- list()
            for (i in which(terms$fut_matu < deliv_bonds$curr_cp_dt)){
              net_basis[[k]][[i]] <- dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
                (joint$ctry_1*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i]  )}
            for (i in which(terms$fut_matu >= deliv_bonds$curr_cp_dt)){
              net_basis[[k]][[i]] <-  dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
                (joint$ctry_1*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] +
                   deliv_bonds$coupon[i]/deliv_bonds$cp_freq[i]*Nomi*exp(fwd_2[[i]]*deliv_bonds$res_term_2[i]) ) }
            net_basis[[k]] <- do.call(cbind, net_basis[[k]])
            ctd_conf_2[[k]] <- ctd_conf[[k]] <- apply(net_basis[[k]], 1, which.min)
            ctd_conf_2[[k]][ctd_conf[[k]] != k] <- 0
          }

          ctd_conf_2 <- do.call(cbind, ctd_conf_2)

          for (i in 1:length(ctd)){
            ctd_conf_2[ctd_conf_2[, i] == i, -i] <- 0}

          ctd_conf_2[ctd_conf_2 > 1] <- 1

          ytm_f <- ytm*ctd_conf_2
          ytm_f <- rowSums(ytm_f)


          deliv_bonds <- bind_rows(deliverables_2[deliverables_2$cp_freq == 1,] %>%
                                     mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu < ., . %m-% years(1), .))) %>%
                                     mutate(curr_cp_dt = prev_cp_dt %m+% years(1)),
                                   deliverables_2[deliverables_2$cp_freq == 2, ] %>%
                                     mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu - . < - months(6), . %m-% years(1),
                                                                             ifelse(option_matu - . < 0, . %m-% months(6), .)))) %>%
                                     mutate(curr_cp_dt = prev_cp_dt %m+% months(6)))

          if(day_count_conv_2 == 1){
            terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                        as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")),
                                      res_term = as.numeric(fut_matu - option_matu)/
                                        as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
          } else if(day_count_conv_2 == 2){
            terms <- terms %>% mutate(option_term = as.numeric(option_matu - start_date)/360,
                                      res_term = as.numeric(fut_matu - option_matu)/360)
          } else if(day_count_conv_2 == 3){
            terms <- terms %>% mutate(option_term =as.numeric(option_matu - start_date)/365,
                                      res_term = as.numeric(fut_matu - option_matu)/365)
          } else {
            terms <- terms %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                      stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                      plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                         ceiling_date(start_date, "months")))/30),
                                      option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                      stub_1_res = max(0, 30 - as.numeric(format(option_matu + sett_2, "%d"))),
                                      stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                      plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                            ceiling_date(option_matu + sett_2, "months") )/30),
                                      res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

          terms <- terms %>% mutate(res_term_2 = 0)

          rate_table <- data.frame(term = terms$option_term + c(0, terms$res_term), rates = c(r, r_2)) %>%
            mutate(d_fact = term*rates)

          fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

          bef <- which(terms$fut_matu < deliv_bonds$curr_cp_dt)
          aft <- which(terms$fut_matu >= deliv_bonds$curr_cp_dt)

          if(day_count_conv_2 == 1){
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_2)/
                                                                     as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") )),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_2)/
                                                                     as.numeric(ceiling_date(terms$fut_matu, "year") - floor_date(terms$fut_matu, "year") ),
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2 ))
          } else if(day_count_conv_2 == 2) {
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_2)/360),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_2)/360,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          } else if(day_count_conv_2 == 3){
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(acc_matu = Nomi*coupon/cp_freq*as.numeric(terms$fut_matu - prev_cp_dt - sett_2)/365),
                                     deliv_bonds[aft, ] %>% mutate(res_term_2 = as.numeric(terms$fut_matu - curr_cp_dt - sett_2)/365,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          } else{
            deliv_bonds <- bind_rows(deliv_bonds[bef, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(prev_cp_dt + sett_2, "%d"))),
                                                                   stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                                   plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                     ceiling_date(prev_cp_dt + sett_2, "months") )/30),
                                                                   acc_matu = Nomi*coupon/cp_freq/360*(stub_2 + stub_1 + max(0, plain_months)*30)),
                                     deliv_bonds[aft, ] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(curr_cp_dt + sett_2, "%d"))),
                                                                   stub_2 = min(30, as.numeric(format(terms$fut_matu, "%d"))),
                                                                   plain_months = round(as.numeric(floor_date(terms$fut_matu, "months") -
                                                                                                     ceiling_date(curr_cp_dt + sett_2, "months") )/30),
                                                                   res_term_2 = (stub_1 + stub_2 + max(0, plain_months)*30)/360,
                                                                   acc_matu = Nomi*coupon/cp_freq*res_term_2))
          }

          deliv_bonds <- deliv_bonds[match(deliverables_2$ISIN, deliv_bonds$ISIN), ]

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

            true_cp_dt[[i]] <- c(head(true_cp_dt[[i]], -1) + sett_2, tail(true_cp_dt[[i]], 1))
            true_cp_dt[[i]] <- c(terms$option_matu, true_cp_dt[[i]])

            if(day_count_conv_2 == 1){
              if(deliv_bonds$cp_freq[i] == 1){ stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
                as.numeric(deliv_bonds$curr_cp_dt[i] - deliv_bonds$prev_cp_dt[i])
              } else {stub[[i]] <- as.numeric(diff(head(true_cp_dt[[i]], 2)))/
                as.numeric(ceiling_date(deliv_bonds$prev_cp_dt[i] + sett_2, "year") -
                             floor_date(deliv_bonds$prev_cp_dt[i] + sett_2, "year")) }
              cp_dt_2[[i]] <- stub[[i]] + c(0, seq(1, length(true_cp_dt[[i]]) - 2)/deliv_bonds$cp_freq[i])
            } else if(day_count_conv_2 == 2){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/360
            } else if(day_count_conv_2 == 3){ cp_dt_2[[i]] <- tail(as.numeric(true_cp_dt[[i]] - first(true_cp_dt[[i]])), -1)/365
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

          ytm <- list()
          for (i in 1:nrow(deliv_bonds)){
            ytm[[i]] <- tri(joint$ctry_2)}

          ytm <- do.call(cbind, ytm) %>% data.frame %>% rename_with(~c(deliv_bonds$ISIN))

          ctd <- list()
          for (i in 1:nrow(deliv_bonds)){
            ctd[[i]] <- replicate(nrow(deliv_bonds), ytm[, i] - deliv_bonds$ytm[i]) +
              t(replicate(nrow(ytm), deliv_bonds$ytm)) }

          net_basis <- ctd_conf <- ctd_conf_2 <- list()

          for (k in 1:length(ctd)){
            net_basis[[k]] <- list()
            for (i in which(terms$fut_matu < deliv_bonds$curr_cp_dt)){
              net_basis[[k]][[i]] <- dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
                (joint$ctry_2*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i]  )}
            for (i in which(terms$fut_matu >= deliv_bonds$curr_cp_dt)){
              net_basis[[k]][[i]] <-  dirty(ctd[[k]][, i])*exp(fwd_1*terms$res_term) -
                (joint$ctry_2*deliv_bonds$conv_factor[i] + deliv_bonds$acc_matu[i] +
                   deliv_bonds$coupon[i]/deliv_bonds$cp_freq[i]*Nomi*exp(fwd_2[[i]]*deliv_bonds$res_term_2[i]) ) }
            net_basis[[k]] <- do.call(cbind, net_basis[[k]])
            ctd_conf_2[[k]] <- ctd_conf[[k]] <- apply(net_basis[[k]], 1, which.min)
            ctd_conf_2[[k]][ctd_conf[[k]] != k] <- 0
          }

          ctd_conf_2 <- do.call(cbind, ctd_conf_2)

          for (i in 1:length(ctd)){
            ctd_conf_2[ctd_conf_2[, i] == i, -i] <- 0}

          ctd_conf_2[ctd_conf_2 > 1] <- 1

          ytm_f_2 <- ytm*ctd_conf_2
          ytm_f_2 <- rowSums(ytm_f_2)

          joint_f <- data.frame(ctry_1 = ytm_f, ctry_2 = ytm_f_2) %>%
            mutate(spread = (ctry_2 - ctry_1)*1e4)

          dens_ex <- data.frame(spread = joint_f$spread)

          bond_term <- list()
          raw_bond_term <- as.numeric(as.Date(c(ctd_matu_1, ctd_matu_2)) - as.Date(fut_matu))/365
          for (i in 1:length(raw_bond_term)){
            if(raw_bond_term[i] <  2.5){ bond_term[[i]] <- round(raw_bond_term[i]/2)*2
            } else if (raw_bond_term[i] < 5.5) { bond_term[[i]] <- round(raw_bond_term[i]/5)*5
            } else {bond_term[[i]] <- round(raw_bond_term[i]/10)*10}
          }

          spread_plot <- ggplot(dens_ex, aes(x = spread)) + geom_density() +
            labs(title = paste0(ctry_2, "-", ctry_1, " bond yield spread on ", option_matu,
                                " as of ", start_date), subtitle = "Risk Neutral Probability Density",
                 x = paste0("spread (bp) between ", ctry_2, " ", round(bond_term[[2]]), "-y bond and ",
                            ctry_1, " ", round(bond_term[[1]]), "-y bond yields"), y = "probability density") +
            theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm")) +
            theme_bw()

          return(list(spread_distrib = dens_ex, spread_plot = spread_plot))

        })

      } else {message("input dates are not consistent")}
    } else {message("inputs do not have the required length")}
  } else{message("impossible to retrieve densities for both ctrys")}
}
