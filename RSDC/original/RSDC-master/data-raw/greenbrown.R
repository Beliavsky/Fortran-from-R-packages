# Build the packaged dataset "greenbrown" from the source Excel.
# The workbook lives in data-raw/ (build-ignored, not shipped on CRAN); run this
# script from the package root.

# Only needed for building the data:
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required to rebuild this dataset. Install it with: install.packages(\"readxl\")")

file <- "data-raw/green-brown-ptf.xlsx"

greenbrown <- readxl::read_excel(file)

# Clean/types
greenbrown$DATE <- as.Date(greenbrown$DATE)

# Ensure it's a base data.frame
greenbrown <- as.data.frame(greenbrown)

# Order by DATE ascending
greenbrown <- greenbrown[order(greenbrown$DATE), ]

# Save as compressed .rda in data/
usethis::use_data(greenbrown, overwrite = TRUE, compress = "xz")
