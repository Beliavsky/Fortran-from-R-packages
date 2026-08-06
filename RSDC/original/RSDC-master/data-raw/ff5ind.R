# Build data/ff5ind.rda: daily returns of five Fama-French industry
# portfolios with the MCCC index and the VIX, aligned on trading days,
# 2005-01-03 -> 2025-06-30 (the MCCC 2025 release ends June 2025).
#
# Sources (all downloaded here; a local cache is used when present):
#  - Industry portfolios: Kenneth R. French Data Library, "10 Industry
#    Portfolios [Daily]" (value-weighted returns in percent); industry
#    definitions from Fama & French (1997). Data (c) Fama and French.
#  - MCCC: Sentometrics Media Climate Change Concerns workbook (2025
#    release), monthly Aggregate series, forward-filled across the trading
#    days of each month (the package's convention, see ?mccc); Ardia,
#    Bluteau, Boudt & Inghelbrecht (2023).
#  - VIX: CBOE Volatility Index, daily close, via FRED (series VIXCLS).
#
# Run from the package root: source("data-raw/ff5ind.R")

INDUSTRIES <- c("Manuf", "Enrgy", "HiTec", "Hlth", "Utils")
FIRST <- as.Date("2005-01-01")
LAST  <- as.Date("2025-06-30")

## ---- 1. Fama-French 10 industry daily returns --------------------------------
ff_cache <- "data-raw/cache_ff10_daily.csv"
if (!file.exists(ff_cache)) {
  tmp <- tempfile(fileext = ".zip")
  download.file(paste0("https://mba.tuck.dartmouth.edu/pages/faculty/",
                       "ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip"), tmp)
  file.copy(unzip(tmp, exdir = tempdir()), ff_cache)
}
raw <- readLines(ff_cache)
start <- grep("^\\s*,NoDur", raw)[1]
end   <- grep("Average Equal Weighted", raw)[1]
ff <- read.csv(text = paste(raw[start:(end - 2)], collapse = "\n"))
names(ff)[1] <- "DATE"
ff$DATE <- as.Date(as.character(ff$DATE), "%Y%m%d")
ff <- ff[!is.na(ff$DATE), c("DATE", INDUSTRIES)]
ff[, -1][ff[, -1] <= -99] <- NA          # French codes missing values as -99

## ---- 2. VIX (FRED) ------------------------------------------------------------
vix_cache <- "data-raw/cache_vixcls.csv"
if (!file.exists(vix_cache))
  download.file("https://fred.stlouisfed.org/graph/fredgraph.csv?id=VIXCLS",
                vix_cache)
vix <- read.csv(vix_cache)
names(vix) <- c("DATE", "vix")
vix$DATE <- as.Date(vix$DATE)
vix$vix  <- suppressWarnings(as.numeric(vix$vix))   # holidays are "."

## ---- 3. MCCC monthly Aggregate (Sentometrics 2025 release) --------------------
mccc_cache <- "data-raw/cache_mccc_2025.xlsx"
if (!file.exists(mccc_cache))
  download.file(paste0("https://www.dropbox.com/scl/fi/uucc6401uje293ofc3ahq/",
                       "Sentometrics_US_Media_Climate_Change_Index.xlsx",
                       "?dl=1&rlkey=jvgb6xg9w4ctdz5cdl6qun5md"),
                mccc_cache, mode = "wb")
mm <- suppressMessages(readxl::read_excel(mccc_cache,
                                          sheet = "2025 update monthly",
                                          skip = 6, col_names = FALSE))
mccc_m <- data.frame(
  month = as.Date(suppressWarnings(as.numeric(mm[[1]])), origin = "1899-12-30"),
  mccc  = suppressWarnings(as.numeric(mm[[2]])))
mccc_m <- mccc_m[stats::complete.cases(mccc_m), ]
stopifnot(max(mccc_m$month) >= as.Date("2025-06-01"))

## ---- 4. Merge on trading days, window, forward-fill MCCC ----------------------
m <- merge(ff, vix, by = "DATE")
m <- m[m$DATE >= FIRST & m$DATE <= LAST & stats::complete.cases(m), ]
idx <- match(as.Date(format(m$DATE, "%Y-%m-01")), mccc_m$month)
stopifnot(!anyNA(idx))                    # every trading day has its month
m$mccc <- mccc_m$mccc[idx]
ff5ind <- m[, c("DATE", INDUSTRIES, "mccc", "vix")]
rownames(ff5ind) <- NULL

stopifnot(nrow(ff5ind) == 5155,           # matches the frozen K = 5 analyses
          !anyNA(ff5ind),
          identical(range(ff5ind$DATE), as.Date(c("2005-01-03", "2025-06-30"))))
save(ff5ind, file = "data/ff5ind.rda", compress = "xz")
cat("data/ff5ind.rda written:", nrow(ff5ind), "rows x", ncol(ff5ind), "cols\n")
