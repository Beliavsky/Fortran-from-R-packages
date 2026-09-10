# SPDX-License-Identifier: MIT
# Independent reference values for the Fortran lubridate tests.

stopifnot(as.integer(as.Date("1970-01-01")) == 0L)
stopifnot(as.integer(as.Date("1969-12-31")) == -1L)
stopifnot(format(as.Date("1970-01-01"), "%w") == "4")
stopifnot(format(as.Date("2021-01-01"), "%G-%V-%u") == "2020-53-5")
stopifnot(as.integer(format(as.Date("2024-03-01"), "%j")) == 61L)

if (requireNamespace("lubridate", quietly = TRUE)) {
  stopifnot(lubridate::`%m+%`(
    lubridate::ymd("2024-01-31"), lubridate::months(1)
  ) == as.Date("2024-02-29"))
  start <- lubridate::ymd_hms("1970-01-01 00:00:00", tz = "UTC")
  stopifnot(lubridate::int_length(lubridate::interval(start, start + 7200)) == 7200)
  rounded <- lubridate::round_date(
    lubridate::ymd_hms("2024-06-19 12:34:40", tz = "UTC"), "minute"
  )
  stopifnot(format(rounded, "%H:%M:%S", tz = "UTC") == "12:35:00")
}

cat("All independent lubridate reference invariants passed.\n")
