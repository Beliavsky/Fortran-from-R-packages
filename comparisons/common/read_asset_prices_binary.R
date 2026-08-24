read_asset_prices_binary <- function(filename) {
  con <- file(filename, open="rb")
  on.exit(close(con))

  magic <- readChar(con, nchars=8L, useBytes=TRUE)
  if (!identical(magic, "APRICE01")) stop("invalid asset-price binary signature")
  nrow <- readBin(con, integer(), n=1L, size=4L, endian="little")
  ncol <- readBin(con, integer(), n=1L, size=4L, endian="little")
  symbol_width <- readBin(con, integer(), n=1L, size=4L, endian="little")
  if (nrow < 2L || ncol < 1L || symbol_width != 16L) stop("invalid binary dimensions")

  symbols <- character(ncol)
  for (j in seq_len(ncol)) {
    symbols[j] <- trimws(readChar(con, nchars=symbol_width, useBytes=TRUE))
  }
  date_values <- readBin(con, integer(), n=3L*nrow, size=4L, endian="little")
  date_values <- matrix(date_values, ncol=3L, byrow=TRUE)
  dates <- as.Date(sprintf("%04d-%02d-%02d", date_values[,1], date_values[,2], date_values[,3]))
  values <- readBin(con, double(), n=nrow*ncol, size=8L, endian="little")
  if (length(values) != nrow*ncol) stop("truncated binary price data")
  prices <- matrix(values, nrow=nrow, ncol=ncol, dimnames=list(NULL, symbols))
  list(dates=dates, symbols=symbols, prices=prices)
}
