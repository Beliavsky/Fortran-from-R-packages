## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE,        # CRAN fails handling external data too often. That's why we have to prevent running the code here. 
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
#  library(treasuryTR)
#  
#  yield_1y <- get_yields("DGS1")
#  Sys.sleep(1)
#  yield_10y <- get_yields("DGS10")
#  Sys.sleep(1)
#  yield_20y <- get_yields("DGS20")
#  Sys.sleep(1)
#  
#  tr_1y <- total_return(yield_1y, maturity = 1, scale = 261)
#  tr_10y <- total_return(yield_10y, maturity = 10, scale = 261)
#  tr_20y <- total_return(yield_20y, maturity = 20, scale = 261)
#  
#  head(cbind.xts(tr_1y, tr_10y, tr_20y))
#  
#  library(PerformanceAnalytics)
#  table.AnnualizedReturns(cbind.xts(tr_1y, tr_10y, tr_20y), Rf = tr_1y, scale = 262)

## ---- fig.width=5-------------------------------------------------------------
#  yield_10y_monthly <- yield_10y[endpoints(yield_10y, on = "months", k = 1)]
#  
#  tr_10y_monthly <- total_return(yield_10y_monthly, 10, scale = 12)
#  
#  performance_10y <- cumprod(1+tr_10y_monthly[-1])-1
#  
#  plot(performance_10y)
#  
#  table.AnnualizedReturns(tr_10y_monthly, scale = 12)

## -----------------------------------------------------------------------------
#  library(dplyr)
#  
#  yield_10y_df <- get_yields("DGS10", format_out = "tibble")
#  
#  tr_10y_df <- yield_10y_df %>%
#    mutate(TR = total_return(DGS10, maturity = 10))
#  
#  tr_10y_df %>%
#    filter(!is.na(TR)) %>%
#    summarise(mu = mean(TR)*262,
#              sigma = sd(TR)*sqrt(262))

## -----------------------------------------------------------------------------
#  tr_10y_df_stepbystep <- yield_10y_df %>%
#    mutate(mod_duration = mod_duration(DGS10, 10),
#           convexity = convexity(DGS10, 10),
#           TR = total_return(DGS10, maturity = 10,
#                             mdur = mod_duration,
#                             convex = convexity))
#  
#  tail(tr_10y_df_stepbystep)

## ---- fig.width=5-------------------------------------------------------------
#  library(dataseries)
#  library(tidyr)
#  library(ggplot2)
#  
#  swiss_yields <- ds(c("ch_snb_rendoblim.1j",
#                       "ch_snb_rendoblim.10j",
#                       "ch_snb_rendoblim.20j"))
#  
#  swiss_tr <- swiss_yields %>%
#    mutate(TR1 = total_return(ch_snb_rendoblim.1j/100, maturity = 1, scale = 12),
#           TR10 = total_return(ch_snb_rendoblim.10j/100, maturity = 10, scale = 12),
#           TR20 = total_return(ch_snb_rendoblim.20j/100, maturity = 20, scale = 12)) %>%
#    select(time, starts_with("TR")) %>%
#    pivot_longer(cols = -time) %>%
#    filter(!is.na(value)) %>%
#    arrange(name, time)
#  
#  swiss_tr %>%
#    group_by(name) %>%
#    mutate(performance = cumprod(1+value)-1) %>%
#    ggplot(aes(x = time, y = performance*100, color = name)) +
#    geom_line() +
#    scale_y_continuous() +
#    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
#    labs(title = "Cumulative performance since 1962 of Swiss Confederation Bonds",
#         x = "", y = "%", color = "") +
#    theme_classic() +
#    theme(legend.position = "top",
#          plot.title = element_text(hjust = 0.5),
#          axis.text.x = element_text(angle = 90, vjust = 0.5))
#  
#  swiss_tr %>%
#    group_by(name, year = format(time, "%Y")) %>%
#    summarise(TR = prod(1+value)-1) %>%
#    ggplot(aes(x = year, y = TR*100, fill = name)) +
#    geom_col(position = position_dodge2(), alpha = 0.9) +
#    scale_y_continuous() +
#    labs(title = "TR per calendar year of Swiss Confederation Bonds",
#         x = "", y = "%", fill = "") +
#    theme_classic() +
#    theme(legend.position = "top",
#          plot.title = element_text(hjust = 0.5),
#          axis.text.x = element_text(angle = 90, vjust = 0.5),
#          panel.grid.major.x = element_line())

