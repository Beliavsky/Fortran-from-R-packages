## code to prepare dataset for Candelon and Moura (2024)
rm(list = ls())
cat("\014")
library(dplyr)
library(readxl)
library(seasonal)
library(lubridate)
library(YieldCurve)
library(tibble)
library(xts)

####### 1) Trade Flows #######
Economies <- c("China", "Brazil", "Mexico", "Uruguay")
TradeData <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "TradeData")

years <- 1948:2019
C <- length(Economies)
T_dim <- length(years)

TradeFlows <- list()

for (j in 1:C) {
  # Create a regular data.frame instead of tibble
  Tradeflows_CS <- data.frame(matrix(NA_real_, nrow = C, ncol = T_dim))
  colnames(Tradeflows_CS) <- years
  rownames(Tradeflows_CS) <- Economies

  for (i in 1:C) {
    if (Economies[j] == Economies[i]) {
      Tradeflows_Bilateral <- rep(0, times = T_dim)
    } else {
      Tradeflows_Bilateral <- TradeData %>%
        filter(
          `Country Name` == Economies[j],
          `Counterpart Country Name` == Economies[i]
        ) %>%
        mutate(across(all_of(as.character(years)), as.numeric)) %>%
        summarise(across(
          all_of(as.character(years)),
          ~ if (any(is.na(.))) NA_real_ else sum(.)
        )) %>%
        unlist(use.names = FALSE)
    }

    # Assign numeric vector directly (works in data.frame)
    Tradeflows_CS[i, ] <- as.numeric(Tradeflows_Bilateral)
  }

  # Convert to tibble at the end
  TradeFlows[[Economies[j]]] <- as_tibble(Tradeflows_CS)
  rownames(TradeFlows[[Economies[j]]]) <- Economies
}


####### 2) Global unspanned factors #######
GlobalMacroData <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "MacroGlobal")

Gl_Eco_Act <- GlobalMacroData$Glob_BusCyc / 100000000
Gl_Inflation <- (GlobalMacroData$Global_CPI_OECD /
  dplyr::lag(GlobalMacroData$Global_CPI_OECD, 12) - 1) * 100

GlobalMacro <- GlobalMacro[, colSums(is.na(GlobalMacro)) < nrow(GlobalMacro)]
GlobalMacroData$`Sample Date` <- format(as.Date(GlobalMacroData$`Sample Date`), "%d-%m-%Y")
colnames(GlobalMacro) <- GlobalMacroData$`Sample Date`[-(1:12)]


####### 3) Domestic unspanned factors #######
DomesticMacroData <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "MacroDomestic")

# a) Economic activity (remove NAs and adjust sample dates)
BR_Eco_Act <- na.omit(DomesticMacroData$BR_Eco_Act)[-(1:12)]
CH_Eco_Act <- na.omit(DomesticMacroData$CH_Eco_Act)[-(1:12)]
ME_Eco_Act <- na.omit(DomesticMacroData$ME_Eco_Act)[-(1:12)]

# a1) Get Uruguayan data
UR_DataGDP_Quarterly <- DomesticMacroData$UR_Eco_Act
TimeSpan <- DomesticMacroData$`Sample Date`
# GDP
GDPquarterly <- log(na.omit(UR_DataGDP_Quarterly))
InitialYear <- 2002
GDP <- ts(GDPquarterly, start = c(InitialYear, 1), frequency = 4) # GDP series starts in 2002Q1
GDPquarterlySA <- seas(GDP)$data # Seasonally adjusted series
GDPquarterlySA <- GDPquarterlySA[, 3] # Seasonally adjusted series (final)
# Interpolation
LogMonGDP <- spline(1:length(GDP), GDPquarterlySA, n = length(TimeSpan))$y
# YoY growth
YoY_growth <- c()
for (j in 1:(length(LogMonGDP) - 12)) {
  YoY_growth[j + 12] <- LogMonGDP[j + 12] - LogMonGDP[j]
}

UR_GDP_YoY_Monthly <- data.frame(TimeSpan, YoY_growth)
UR_GDP_YoY_Monthly_Adj <- subset(UR_GDP_YoY_Monthly, TimeSpan >= as.Date("2004-06-01") & TimeSpan <= as.Date("2020-01-01"))
UR_Eco_Act <- 100 * UR_GDP_YoY_Monthly_Adj$YoY_growth


# b) Inflation (YoY growth rate)
BR_Inflation <- (DomesticMacroData$BR_CPI /
  dplyr::lag(DomesticMacroData$BR_CPI, 12) - 1) * 100

CH_Inflation <- (DomesticMacroData$CH_CPI /
  dplyr::lag(DomesticMacroData$CH_CPI, 12) - 1) * 100

ME_Inflation <- (DomesticMacroData$ME_CPI /
  dplyr::lag(DomesticMacroData$ME_CPI, 12) - 1) * 100

UR_Inflation <- (DomesticMacroData$UR_CPI /
  dplyr::lag(DomesticMacroData$UR_CPI, 12) - 1) * 100

# Remove NAs and adjust sample dates
BR_Inflation <- na.omit(BR_Inflation)
CH_Inflation <- na.omit(CH_Inflation)
ME_Inflation <- na.omit(ME_Inflation)
UR_Inflation <- na.omit(UR_Inflation)

# Adjust outputs
DomMacro <- rbind(CH_Eco_Act, CH_Inflation, BR_Eco_Act, BR_Inflation, ME_Eco_Act, ME_Inflation, UR_Eco_Act, UR_Inflation)
row.names(DomMacro) <- c(
  "Eco_Act China", "Inflation China", "Eco_Act Brazil", "Inflation Brazil", "Eco_Act Mexico",
  "Inflation Mexico", "Eco_Act Uruguay", "Inflation Uruguay"
)
colnames(DomMacro) <- GlobalMacroData$`Sample Date`[-(1:12)]


####### 4) Bond yields #######

# a) China and Mexico
Yields_CH_all <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "Yields_China")
Yields_ME_all <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "Yields_Mexico")

# Clean output to export
Yields_CH <- Yields_CH_all[, -which(names(Yields_CH_all) == "Period")]
Yields_CH <- t(Yields_CH)
rownames(Yields_CH) <- paste0(rownames(Yields_CH), "_China")

Yields_ME <- Yields_ME_all[, -which(names(Yields_ME_all) == "Period")]
Yields_ME <- t(Yields_ME)
rownames(Yields_ME) <- paste0(rownames(Yields_ME), "_Mexico")

# b) Brazil
Yields_BR_all <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "Yields_Brazil")

# Use a Nelson-Siegel model to fit the missing data: Y120M from June/04 to Ago/05
t0 <- "2004/6/01" # First date of the sample

# Load Data
sheet_names <- excel_sheets("inst/extdata/YieldsBrazilJun04Aug05.xlsx")
list_all <- lapply(sheet_names, function(x) {
  as.data.frame(read_excel("inst/extdata/YieldsBrazilJun04Aug05.xlsx", sheet = x))
})
names(list_all) <- sheet_names

T_dim <- length(list_all)

y <- data.frame(matrix(NA, nrow = T_dim, ncol = 2))
colnames(y) <- c("Time", "Y120M")
y$Time <- seq(as.Date(t0), by = "month", length.out = T_dim)

for (i in 1:T_dim) {
  # Define variables
  rate.10y <- as.matrix(list_all[[i]]$Yields / 100)
  rate.10y <- t(rate.10y)
  colnames(rate.10y) <- list_all[[i]]$MatMonts
  maturity.10y <- as.numeric(colnames(rate.10y)) / 12

  # Estimate and fit the NS model
  NSParameters <- Nelson.Siegel(rate = rate.10y, maturity = maturity.10y)
  NSParameters <- xts(NSParameters, order.by = Sys.Date())

  # store fitted value
  y$Y120M[i] <- as.numeric(NSrates(NSParameters, 120))
}

Yields_BR_all$Y120M[1:T_dim] <- y$Y120M
Yields_BR <- Yields_BR_all[, -which(names(Yields_BR_all) == "Period")]
Yields_BR <- t(Yields_BR)
rownames(Yields_BR) <- paste0(rownames(Yields_BR), "_Brazil")


# c) Uruguay
YieldsDaily_UR <- read_excel("inst/extdata/Dataset_CM_2024.xlsx", sheet = "Yields_Uruguay")

YieldsMonthly_UR <- YieldsDaily_UR %>%
  mutate(Dates = as.Date(Dates)) %>%
  group_by(year = year(Dates), month = month(Dates)) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(-year, -month)

Yields_UR <- YieldsMonthly_UR[, -which(names(YieldsMonthly_UR) == "Dates")]
Yields_UR <- t(Yields_UR)
rownames(Yields_UR) <- paste0(rownames(Yields_UR), "_Uruguay")


# Gather yields to export
Yields <- rbind(Yields_CH, Yields_BR, Yields_ME, Yields_UR)
colnames(Yields) <- colnames(GlobalMacro)


####### 5) RiskFacFull data #######
# Compute country-specific spanned factors
N <- 3
PP <- Spanned_Factors(Yields, Economies, N)

RiskFacFull <- rbind(
  GlobalMacro, # Global factors
  DomMacro[1:2, ], PP[1:3, ], # China's factors
  DomMacro[3:4, ], PP[4:6, ], # Brazil's factors
  DomMacro[5:6, ], PP[7:9, ], # Mexico's factors
  DomMacro[7:8, ], PP[10:12, ]
) # Uruguay's factors


####### 6) GVARFactors data #######
t0 <- "2004-06-01"
tF <- "2020-01-01"
ModelType <- "GVAR multi"

GlobalVar <- c("Gl_Eco_Act", "Gl_Inflation")
DomVar <- c("Eco_Act", "Inflation")
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)

# Transition matrix
Wgvar <- Transition_Matrix(t_First = "2004", t_Last = "2019", Economies, type = "Sample Mean", TradeFlows)

# Build list of variables
MacroDataList <- list()
Period <- as.POSIXct(as.Date(colnames(GlobalMacro), format = "%d-%m-%Y"))
MacroDataList$Global <- as_tibble(data.frame(Period, t(GlobalMacro)))
colnames(MacroDataList$Global)[-1] <- GlobalVar
MacroDataList$China <- as_tibble(data.frame(Period, t(DomMacro[1:2, ])))
colnames(MacroDataList$China)[-1] <- DomVar
MacroDataList$Brazil <- as_tibble(data.frame(Period, t(DomMacro[3:4, ])))
colnames(MacroDataList$Brazil)[-1] <- DomVar
MacroDataList$Mexico <- as_tibble(data.frame(Period, t(DomMacro[5:6, ])))
colnames(MacroDataList$Mexico)[-1] <- DomVar
MacroDataList$Uruguay <- as_tibble(data.frame(Period, t(DomMacro[7:8, ])))
colnames(MacroDataList$Uruguay)[-1] <- DomVar


YieldsDataList <- list()
YieldsDataList$China <- as_tibble(data.frame(Period, t(Yields[1:6, ])))
colnames(YieldsDataList$China)[-1] <- c("Y3M", "Y6M", "Y12M", "Y36M", "Y60M", "Y120M")
YieldsDataList$Brazil <- as_tibble(data.frame(Period, t(Yields[7:12, ])))
colnames(YieldsDataList$Brazil)[-1] <- c("Y3M", "Y6M", "Y12M", "Y36M", "Y60M", "Y120M")
YieldsDataList$Mexico <- as_tibble(data.frame(Period, t(Yields[13:18, ])))
colnames(YieldsDataList$Mexico)[-1] <- c("Y3M", "Y6M", "Y12M", "Y36M", "Y60M", "Y120M")
YieldsDataList$Uruguay <- as_tibble(data.frame(Period, t(Yields[19:24, ])))
colnames(YieldsDataList$Uruguay)[-1] <- c("Y3M", "Y6M", "Y12M", "Y36M", "Y60M", "Y120M")

GVARFactors <- DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelType, MacroDataList, YieldsDataList, Wgvar)

# Output to export
usethis::use_data(TradeFlows, overwrite = TRUE)
usethis::use_data(GlobalMacro, overwrite = TRUE)
usethis::use_data(DomMacro, overwrite = TRUE)
usethis::use_data(Yields, overwrite = TRUE)
usethis::use_data(RiskFacFull, overwrite = TRUE)
usethis::use_data(GVARFactors, overwrite = TRUE)

