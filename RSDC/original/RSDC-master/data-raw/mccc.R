# Build the packaged dataset "mccc" from the monthly source in inst/extdata.
#
# Source: Media Climate Change Concerns (MCCC) index, "Aggregate" series,
# monthly frequency (Ardia, Bluteau, Boudt, Inghelbrecht, 2023). The full
# multi-cluster workbook (MCCC.xlsx) is large and is NOT shipped; only the
# monthly Aggregate column is retained here as a small CSV.
#
# The package convention is a *daily* trading calendar (see greenbrown), so the
# monthly index is forward-filled across the calendar days of each month and
# then aligned to the greenbrown trading dates. This yields a daily covariate
# that is row-for-row mergeable with greenbrown for the TVTP estimation
# workflow (rsdc_estimate(method = "tvtp", X = ...)).

monthly <- utils::read.csv(
  system.file("extdata", "mccc-monthly.csv", package = "RSDC"),
  stringsAsFactors = FALSE
)
monthly$Date <- as.Date(monthly$Date)
monthly <- monthly[is.finite(monthly$Aggregate) & !is.na(monthly$Date), ]
monthly <- monthly[order(monthly$Date), ]

# Forward-fill each monthly value across every calendar day of its month
daily <- do.call(rbind, lapply(seq_len(nrow(monthly)), function(i) {
  a <- as.Date(cut(monthly$Date[i], "month"))            # first day of month
  b <- seq(a, by = "month", length.out = 2L)[2L] - 1L    # last day of month
  data.frame(DATE = seq(a, b, by = "day"), mccc = monthly$Aggregate[i])
}))
daily <- daily[!duplicated(daily$DATE), ]

# Align to the greenbrown trading calendar (inner join on DATE)
gb <- get(load("data/greenbrown.rda"))
mccc <- merge(data.frame(DATE = gb$DATE), daily, by = "DATE")
mccc <- mccc[order(mccc$DATE), ]
rownames(mccc) <- NULL
mccc <- as.data.frame(mccc)

# Save as compressed .rda in data/
usethis::use_data(mccc, overwrite = TRUE, compress = "xz")
