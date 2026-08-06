#' deliv_bonds_charac_bbg
#'
#' @param bbg_ticker a Bloomberg ticker of a bond futures contract, in character format
#' @param date a date for the recovery of the characteristics of the deliverable bonds in the futures contract, in Date format
#'
#' @returns For the bonds in the delivery basket, their ISIN code in character format, their coupon rate and frequency of coupon payment in numeric format (1 for annual and 2 for semi-annual), their maturity date in Date format, their conversion factor, yield to maturity and day count convention in numeric format (1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360). Provided physical delivery of the contract, the same information is also displayed for the current Cheapest-to-Deliver Bond in a second table, with in addition its net basis in currency units of the futures contract, as calculated by Bloomberg, in numeric format
#' @export
#' @import dplyr
#' @import Rblpapi
#' @import tibble
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' deliv_bonds_charac_bbg("UBU6", as.Date("2026-05-19"))
#' }
#'

deliv_bonds_charac_bbg <- function(bbg_ticker, date){

  blpConnect()

  if (length(bbg_ticker) > 1){ message("please enter one ticker at a time")
  } else{

    bbg_ticker <- paste0(bbg_ticker, " Comdty")
    basket <- bds(bbg_ticker, "FUT_DLVRBLE_BNDS_ISINS")

    underlying <- bdp(bbg_ticker, "FUT_NOTL_BOND") %>% filter(FUT_NOTL_BOND != "")

    if(nrow(underlying) == 0){message("please enter a bond futures Bloomberg ticker")
    } else{

      basket <- basket %>% rename_with(~c("conv_factor", "ISIN")) %>%
        mutate_at("ISIN", ~gsub("   Corp", " Corp", .))

      if(nrow(basket) == 1){
        px <- data.frame(bdh(basket$ISIN, "YLD_YTM_MID", start.date = date, end.date = date))
        rownames(px) <- basket$ISIN
      } else {
        px <- do.call(rbind, bdh(basket$ISIN, "YLD_YTM_MID", start.date = date, end.date = date)) }
      if(nrow(px) == 0){
        message("please enter a trading day")
      } else{
        px <- px %>% dplyr::select(-date) %>% rename_with(~"ytm") %>%
          mutate_at("ytm", ~./100) %>% rownames_to_column("ISIN")

        corres <- data.frame(day_count_conv_prov = c(35, 20, 9, 10, 4, 3, 2, 1),
                             day_count_conv = c(4, 4, 3, 1, 4, 3, 2, 1))

        charac <- bdp(basket$ISIN, c("CPN", "CPN_FREQ", "DAY_CNT", "MATURITY", "ID_CUSIP")) %>%
          rename_with(~c("coupon", "cp_freq", "day_count_conv_prov", "matu", "cusip")) %>%
          mutate_at("coupon", ~./100) %>% mutate_at(c("cp_freq", "day_count_conv_prov"), as.numeric) %>%
          rownames_to_column("ISIN") %>% left_join(basket, by = "ISIN") %>% left_join(px, by = "ISIN") %>%
          mutate_at("ISIN", ~gsub(" Corp", "", .)) %>% arrange(matu) %>%
          left_join(corres, by = "day_count_conv_prov") %>%
          dplyr::select(-day_count_conv_prov) %>%
          mutate_at("cusip", ~substr(., 1, nchar(.) - 1))

        ctd <- bdh(bbg_ticker, "FUT_CTD_CUSIP", start.date = date, end.date = date) %>%
          rename_at(2, ~"cusip")

        if (nrow(ctd) == 0){ message("the contract is not physically delivered - no real CtD")
          return(charac)
        } else {
          if(!ctd$cusip%in%charac$cusip){
            ctd_isin <- bdp(bbg_ticker, "FUT_CTD_ISIN") %>% rename_all(~"ISIN")
            ctd <- charac %>% filter(ISIN == ctd_isin$ISIN) %>% dplyr::select(-cusip)
            charac <- charac %>% dplyr::select(-cusip) %>% na.omit()
          } else {
            ctd <- charac %>% filter(cusip == ctd$cusip)
            charac <- charac %>% dplyr::select(-cusip) %>% na.omit()
            ctd <- ctd %>% dplyr::select(-cusip)
          }
          net_basis_ctd <- bdh(bbg_ticker, "FUT_CTD_NET_BASIS", start.date = date, end.date = date) %>%
            dplyr::select(-date) %>% rename_all(~"net_basis")

          ctd <- data.frame(ctd, net_basis_ctd)
        }

        return(list(basket = charac, ctd = ctd))
      }
    }
  }
}
