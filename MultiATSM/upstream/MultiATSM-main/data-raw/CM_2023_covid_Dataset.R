## code to prepare dataset for Candelon and Moura (2023)
rm(list = ls())
cat("\014")
library(dplyr)
library(readxl)
library(mFilter)
source(paste(getwd(), "/data-raw/Aux_functions_CM2023.R", sep = ""))

# Common inputs
Economies <- c("Brazil", "India", "Russia", "Mexico")

####### 1) Trade Flows #######
TradeData <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "TradeData")

years <- 1948:2020
C <- length(Economies)
T_dim <- length(years)

TradeFlows_covid <- list()

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
  TradeFlows_covid[[Economies[j]]] <- as_tibble(Tradeflows_CS)
  rownames(TradeFlows_covid[[Economies[j]]]) <- Economies
}

####### 2) Global factors #######
GlobalData <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "GlobalFactors")

# a)US growth
lambdaHP <- 6.25 * (52)^4
US_Output_growth <- hpfilter(GlobalData$`US_gdp growth`, freq = lambdaHP, type = "lambda")$cycle
# Adjust sample span
US_Output_growth <- US_Output_growth[1:80]
SP500 <- GlobalData$SP500[1:80]

# b) China growth
ChinaMonthlyData <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "China_OECD_data")
ChinaEcoAct <- Monthly2Weekly(ChinaMonthlyData, Economy = "China", EcoVar = "EcoAct")
ChinaEcoAct_covid <- subset(ChinaEcoAct, Date >= as.Date("2020-03-22") & Date <= as.Date("2021-09-26"))
China_Output_growth <- (ChinaEcoAct_covid$YoY_EcoAct_China) * 100

GlobalMacro_covid <- rbind(US_Output_growth, China_Output_growth, SP500)
colnames(GlobalMacro_covid) <- format(as.Date(GlobalData$Date), "%d-%m-%Y")[1:80]

####### 3) Domestic factors #######
# a) Inflation
Inflation <- list()

CPIMonthlyData <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "CPI")

for (i in 1:length(Economies)) {
  Series_CS <- Monthly2Weekly(CPIMonthlyData, Economy = Economies[i], EcoVar = "CPI")
  # Before, applying the HP filter, adjust sample to be compatible with the other variables
  Series_CS_covid <- subset(Series_CS, Date >= as.Date("2020-03-22") & Date <= as.Date("2022-01-09"))
  Inflation[[Economies[i]]] <- 100 * (hpfilter(Series_CS_covid[[3]], freq = lambdaHP, type = "lambda")$cycle)
  # Adjust sample of the paper
  Inflation[[Economies[i]]] <- Inflation[[Economies[i]]][1:80]
}

# b) Economic activity
EcoAct_allCountries <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "EcoAct")
Output_growth <- list()

for (i in 1:length(Economies)) {
  EcoAct_CS <- EcoAct_allCountries[[Economies[i]]]
  Output_growth[[Economies[i]]] <- hpfilter(EcoAct_CS, freq = lambdaHP, type = "lambda")$cycle
  Output_growth[[Economies[i]]] <- Output_growth[[Economies[i]]][1:80]
}

# c) R-rate
Rrate_allCountries <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "R_rate")
R_rate <- list()
for (i in 1:length(Economies)) {
  R_rate_CS <- Rrate_allCountries[[Economies[i]]]
  R_rate[[Economies[i]]] <- hpfilter(R_rate_CS, freq = lambdaHP, type = "lambda")$cycle
  R_rate[[Economies[i]]] <- 100 * (R_rate[[Economies[i]]][1:80])
}

# d) CDS
CDS_all <- read_excel("inst/extdata/Dataset_CM_2023.xlsx", sheet = "CDS")

# Prepare output to save
DomMacro_covid <- rbind(
  Inflation$Brazil, Output_growth$Brazil, CDS_all$Brazil, R_rate$Brazil,
  Inflation$India, Output_growth$India, CDS_all$India, R_rate$India,
  Inflation$Russia, Output_growth$Russia, CDS_all$Russia, R_rate$Russia,
  Inflation$Mexico, Output_growth$Mexico, CDS_all$Mexico, R_rate$Mexico
)

VarLabels <- c("Inflation", "Output_growth", "CDS", "COVID")
Labels_all <- as.vector(outer(VarLabels, Economies, paste))
rownames(DomMacro_covid) <- Labels_all
colnames(DomMacro_covid) <- colnames(GlobalMacro_covid)

# Output to export
usethis::use_data(TradeFlows_covid, overwrite = TRUE)
usethis::use_data(GlobalMacro_covid, overwrite = TRUE)
usethis::use_data(DomMacro_covid, overwrite = TRUE)
