#' Master calculator for all the analytics of volatility smiles required for a date
#'
#' This is all-in-one calculator. The inputs are symbol, date (YYYY-MM-DD), 
#' and quartic config file location, 
#' and the optional external data directory. The data structure and documentation here are really rough.
#' They are used to calcuate teh data needed for the quartic paper.
#' They need to be polished and refined after the quartic paper is released.
#'
#' @param symbol character pointing to the standard option data file
#' @param date_str character in the form of YYYY-MM-DD
#' @param config_file character, config file from the quarter optimx fit
#' @param extdata_dir character, external data directory
#' @param int_rate numeric, the interest rate used to calculate BS implied volatility from market data
#' @param div_yield numeric, the vididend yield used to calculate BS implied volatility from market data
#' @param idx integer, indicating the index of the option chain
#' @param df_day data frame for the day
#' @param master the list structure from the output of \code{ecop.term_master_calculator}
#' @param otype character, option type of p or c
#' @param opt the list structure from the output of \code{ecop.smile_data_calculator}
#'
#' @return The nested list containing all analytics of volatility smiles for a date.
#'         The first level keys are the date strings.
#'         The first level attributes are \code{quartic.config} which is a data frame, 
#'         lists of \code{days, volumes, classes}, and values of \code{undl_price, max_idx}.
#'
#' @keywords term structure
#'
#' @export ecop.term_master_calculator
#' @export ecop.smile_data_calculator
#' @export ecop.term_atm
#'
### <======================================================================>
"ecop.term_master_calculator" <- function(symbol, date_str, int_rate=0, div_yield=0, config_file=NULL, extdata_dir=NULL)
{
    master = list()
    if (!is.null(config_file)) {
        master$quartic.config <- read.csv(config_file)
    }
    df0 <- ecop.read_csv_by_symbol(symbol, extdata_dir=extdata_dir)
    #
    Date = days = NULL # Simply to avoid R CMD check complaining
    
    df <- subset(df0, Date == as.Date(date_str))
    
    master$symbol <- symbol
    master$date_str <- date_str
    master$days <- sort(unique(df[["days"]]))
    master$volumes <- sapply(master$days, function(d) sum(subset(df, days==d)[["VOL"]]))
    master$classes <- as.character(sapply(master$days, function(d) unique(subset(df, days==d)[["CLS"]])))
    master$max_idx <- length(master$days)
    master$undl_price <- head(unique(df$UNDL_PRC),1)
    
    idx_range = seq(1, master$max_idx) 
    dt = list()
    for (idx in idx_range) { 
        call <- ecop.smile_data_calculator(idx, df, master, int_rate, div_yield, otype="c")
        put <- ecop.smile_data_calculator(idx, df, master, int_rate, div_yield, otype="p")
        dt[[idx]] <- list(call=call, put=put, int_rate=int_rate, div_yield=div_yield)
    }
    master$dt = dt
        
    return(master)
}
### <---------------------------------------------------------------------->
#' @rdname ecop.term_master_calculator
"ecop.smile_data_calculator" <- function(idx, df_day, master, int_rate, div_yield, otype)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    d <- master$days[idx]
    cls <- as.character(master$classes[idx])
    stopifnot (cls %in% c("SPX", "SPXW"))
    
    ttm <- if (cls=="SPXW") d/365 else (d-1/3)/365 # 1/3 of day shorter for AM settled

    # --------------------------------------------
    Date = days = PC = L_ASK = NULL # Simply to avoid R CMD check complaining

    df <- subset(df_day, days==d & PC==toupper(otype) & L_ASK > 0)
    K_raw <- df$STRK_PRC
    IV_raw <- ecop.bs_implied_volatility(df$L_ASK, K_raw, master$undl_price, ttm=ttm, 
                                     int_rate=int_rate, div_yield=div_yield, otype=otype)
    IV2_raw <- ecop.bs_implied_volatility(df$L_BID, K_raw, master$undl_price, ttm=ttm, 
                                      int_rate=int_rate, div_yield=div_yield, otype=otype)

    k_raw <- log(K_raw/master$undl_price)
    valid <- which(IV_raw > 0.0)
    valid2 <- which(IV2_raw > 0.0)
    
    IV <- IV_raw[valid]
    IV2 <- IV2_raw[valid2]
    k <- k_raw[valid]
    k2 <- k_raw[valid2]
    
    opt = list(df=df, undl_price=master$undl_price, 
               K_raw=K_raw, k_raw=k_raw,
               IV_raw=IV_raw, IV2_raw, 
               IV=IV, IV2=IV2, 
               k=k, k2=k2,
               days=d,
               ttm=ttm,
               class=cls,
               int_rate=int_rate, 
               div_yield=div_yield)

    return(opt)
}
### <---------------------------------------------------------------------->
#' @rdname ecop.term_master_calculator
"ecop.term_atm" <- function(opt) {

    # --------------------------------------------
    k = IV = NULL # Simply to avoid R CMD check complaining

    atm <- which(opt$k >= -0.02 & opt$k <= 0.0)
    if (length(atm) <= 4) {
        atm <- which(opt$k >= -0.05 & opt$k <= 0.0)
    }
    
    k_atm <- opt$k[atm]
    IV_atm <- opt$IV[atm]
    IV2_atm <- opt$IV2[atm]
    N <- length(IV_atm)
    atm_skew <- (IV_atm[1]-IV_atm[N])/(k_atm[1]-k_atm[N])
    atm_vol <- IV_atm[N] 
    min_vol = min(IV[k > 0])
    
    dat = list(k=k_atm, IV=IV_atm, IV2=IV2_atm, N=N, 
               atm_skew=atm_skew, 
               atm_vol=atm_vol, 
               min_vol=min_vol)
    return(dat)
    
}
### <---------------------------------------------------------------------->
