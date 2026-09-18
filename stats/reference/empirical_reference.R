# Regenerate deterministic fixtures in test/test_empirical_distributions.f90.

x <- c(1, 1, 2, 4, 4, 4, 7)
query <- c(0, 1, 1.5, 4, 10)
model <- ecdf(x)
cat("ECDF", format(model(query), digits = 17), "\n")

h <- hist(x, breaks = c(0, 2, 5, 8), plot = FALSE)
cat("HIST_COUNTS", h$counts, "\n")
cat("HIST_DENSITY", format(h$density, digits = 17), "\n")
h <- hist(x, breaks = c(0, 2, 5, 8), right = FALSE, plot = FALSE)
cat("HIST_LEFT_COUNTS", h$counts, "\n")

bandwidth_sample <- c(1, 2, 4, 8, 16, 32)
cat("BW_NRD0", format(bw.nrd0(bandwidth_sample), digits = 17), "\n")
cat("BW_NRD", format(bw.nrd(bandwidth_sample), digits = 17), "\n")

grid <- seq(-1, 9, length.out = 6)
bandwidth <- 1.25
direct_gaussian <- rowMeans(dnorm(outer(grid, x, "-") / bandwidth)) / bandwidth
cat("DENSITY_GRID", format(grid, digits = 17), "\n")
cat("DENSITY_DIRECT", format(direct_gaussian, digits = 17), "\n")
fft_density <- density(x, bw = bandwidth, from = -1, to = 9, n = 6)
cat("DENSITY_R_FFT", format(fft_density$y, digits = 17), "\n")
