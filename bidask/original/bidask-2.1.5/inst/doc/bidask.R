## ---- include = FALSE---------------------------------------------------------
data.table::setDTthreads(1)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  out.width="100%",
  dpi = 300,
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(bidask)

## ---- results='hide'----------------------------------------------------------
library(dplyr)
library(crypto2)
df <- crypto_list(only_active=TRUE) %>%
  filter(symbol %in% c("BTC", "ETH")) %>%
  crypto_history(start_date = "20200101", end_date = "20221231")

## -----------------------------------------------------------------------------
head(df)

## -----------------------------------------------------------------------------
df %>%
  mutate(yyyy = format(timestamp, "%Y")) %>%
  group_by(symbol, yyyy) %>%
  arrange(timestamp) %>%
  summarise("EDGE" = edge(open, high, low, close))

## -----------------------------------------------------------------------------
library(ggplot2)
df %>%
  group_by(symbol) %>%
  arrange(timestamp) %>%
  mutate("EDGE (rolling)" = edge_rolling(open, high, low, close, width = 30)) %>%
  ggplot(aes(x = timestamp, y = `EDGE (rolling)`, color = symbol)) +
  geom_line() +
  theme_minimal()

## -----------------------------------------------------------------------------
df %>%
  group_by(symbol) %>%
  arrange(timestamp) %>%
  mutate("EDGE (expanding)" = edge_expanding(open, high, low, close)) %>%
  ggplot(aes(x = timestamp, y = `EDGE (expanding)`, color = symbol)) +
  geom_line() +
  theme_minimal()

## -----------------------------------------------------------------------------
library(quantmod)
x <- try(getSymbols("MSFT", auto.assign = FALSE, start = "2019-01-01", end = "2022-12-31"), silent = TRUE)
if("try-error" %in% class(x)){
  print("Error in getSymbols; using syntentic data instead")
  x <- sim(5000, units = "day")
}

## -----------------------------------------------------------------------------
class(x)

## -----------------------------------------------------------------------------
head(x)

## -----------------------------------------------------------------------------
spread(x)

## -----------------------------------------------------------------------------
edge(open = x[,1], high = x[,2], low = x[,3], close = x[,4])

## -----------------------------------------------------------------------------
sp <- spread(x, width = endpoints(x, on = "months"))
plot(sp)

## -----------------------------------------------------------------------------
sp <- spread(x, width = 21)
plot(sp)

## -----------------------------------------------------------------------------
setDefaults(getSymbols.av, api.key = "<API-KEY>")

## ---- include=FALSE-----------------------------------------------------------
x <- read.csv(system.file("extdata", "msft.csv", package = "bidask"))
x <- xts(x[,-1], order.by = as.POSIXct(x[,1]))

## -----------------------------------------------------------------------------
x <- x["T09:30/T16:00"]
head(x)

## -----------------------------------------------------------------------------
sp <- spread(x, width = endpoints(x, on = "day"))
plot(sp, type = "b")

## -----------------------------------------------------------------------------
sp <- spread(x, width = endpoints(x, on = "day"), method = c("EDGE", "AR", "CS", "ROLL"))
plot(sp, type = "b", legend.loc = "topright")

