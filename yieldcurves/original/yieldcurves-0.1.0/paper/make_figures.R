# Figure generator for the yieldcurves R Journal paper.
#
# Produces six PDF figures and one LaTeX table under paper/figures/ and
# paper/tables/. Uses real US Treasury yields pulled from FRED. Run from
# the package root with:
#   RSTUDIO_PANDOC=/Applications/quarto/bin/tools Rscript paper/make_figures.R

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
  library(showtext)
  library(scales)
})

font_add("HelveticaNeue",
         regular = "/System/Library/Fonts/Helvetica.ttc",
         bold = "/System/Library/Fonts/Helvetica.ttc",
         italic = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto()
showtext_opts(dpi = 300)

fig_dir <- "paper/figures"
tab_dir <- "paper/tables"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

ok_blue   <- "#0072B2"
ok_orange <- "#E69F00"
ok_green  <- "#009E73"
ok_red    <- "#D55E00"
ok_purple <- "#CC79A7"
ok_sky    <- "#56B4E9"
ok_grey   <- "#999999"

fam <- "HelveticaNeue"

theme_wp <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = fam) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey85"),
      axis.line = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks.length = unit(2.5, "pt"),
      axis.text = element_text(size = base_size, colour = "grey20"),
      axis.title = element_text(size = base_size, colour = "grey20"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1, family = fam),
      legend.key.height = unit(10, "pt"),
      legend.key.width = unit(22, "pt"),
      legend.spacing.x = unit(10, "pt"),
      legend.margin = margin(4, 0, 0, 0),
      plot.margin = margin(6, 10, 6, 6)
    )
}

tex_esc <- function(x) gsub("_", "\\\\_", as.character(x))

# -----------------------------------------------------------------------------
# Load FRED Treasury series, stitch into a wide panel (date, maturities).
# -----------------------------------------------------------------------------
series <- list(
  "0.25" = "DGS3MO", "0.5" = "DGS6MO", "1" = "DGS1",
  "2"    = "DGS2",   "3"   = "DGS3",   "5" = "DGS5",
  "7"    = "DGS7",   "10"  = "DGS10",  "20" = "DGS20",
  "30"   = "DGS30"
)

load_series <- function(fname, col_name) {
  df <- read.csv(file.path("paper/data", paste0(fname, ".csv")))
  df$date <- as.Date(df$observation_date)
  yield <- suppressWarnings(as.numeric(df[[fname]]))
  data.frame(date = df$date, yield = yield / 100,
             stringsAsFactors = FALSE)
}

panel_list <- lapply(names(series), function(m) {
  s <- series[[m]]
  df <- load_series(s)
  df$maturity <- as.numeric(m)
  df
})
panel <- do.call(rbind, panel_list)
panel <- panel[!is.na(panel$yield), ]
panel$maturity <- factor(panel$maturity, levels = as.numeric(names(series)))

# Wide: rows = dates, cols = maturities. Require all 10 maturities present.
to_wide <- function(long) {
  mats <- as.numeric(as.character(unique(long$maturity)))
  dates <- sort(unique(long$date))
  mat_cols <- sort(mats)
  wide <- matrix(NA_real_, nrow = length(dates), ncol = length(mat_cols),
                 dimnames = list(as.character(dates),
                                 as.character(mat_cols)))
  for (i in seq_len(nrow(long))) {
    wide[as.character(long$date[i]),
         as.character(as.numeric(as.character(long$maturity[i])))] <-
      long$yield[i]
  }
  keep <- complete.cases(wide)
  wide <- wide[keep, , drop = FALSE]
  wide
}
wide <- to_wide(panel)

# Daily panel is large; reduce to month-end observations for PCA and for
# visual snapshots.
month_end <- function(wide_mat) {
  dates <- as.Date(rownames(wide_mat))
  ym <- format(dates, "%Y-%m")
  last_of <- tapply(seq_along(dates), ym, max)
  wide_mat[as.integer(last_of), , drop = FALSE]
}
wide_m <- month_end(wide)
mats_num <- as.numeric(colnames(wide_m))
dates_m <- as.Date(rownames(wide_m))

# Pick snapshot dates: pre-pandemic, peak, current.
pick <- function(date) {
  d <- as.Date(date)
  idx <- which.min(abs(dates_m - d))
  dates_m[idx]
}
date_pre  <- pick("2019-12-31")
date_peak <- pick("2023-10-31")
date_now  <- dates_m[length(dates_m)]

snap <- function(d) {
  list(date = d,
       maturities = mats_num,
       rates = as.numeric(wide_m[as.character(d), ]))
}
s_pre  <- snap(date_pre)
s_peak <- snap(date_peak)
s_now  <- snap(date_now)

# -----------------------------------------------------------------------------
# Figure 1: Nelson-Siegel vs Svensson fit at the most recent date.
# -----------------------------------------------------------------------------
ns <- yc_nelson_siegel(s_now$maturities, s_now$rates)
sv <- yc_svensson(s_now$maturities, s_now$rates)

grid <- seq(0.25, 30, length.out = 200)
df1 <- rbind(
  data.frame(maturity = grid,
             rate = 100 * yc_predict(ns, maturities = grid)$rate,
             series = "Nelson-Siegel"),
  data.frame(maturity = grid,
             rate = 100 * yc_predict(sv, maturities = grid)$rate,
             series = "Svensson")
)
df1$series <- factor(df1$series, levels = c("Nelson-Siegel", "Svensson"))
pts1 <- data.frame(maturity = s_now$maturities,
                   rate = 100 * s_now$rates)

p1 <- ggplot() +
  geom_point(data = pts1, aes(x = maturity, y = rate),
             colour = "grey20", size = 2.2, alpha = 0.85) +
  geom_line(data = df1, aes(x = maturity, y = rate,
                            colour = series, linetype = series),
            linewidth = 0.8) +
  scale_colour_manual(values = c("Nelson-Siegel" = ok_blue,
                                 "Svensson" = ok_red)) +
  scale_linetype_manual(values = c("Nelson-Siegel" = "solid",
                                   "Svensson" = "longdash")) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 30)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Maturity (years, log scale)", y = "Yield") +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig1_ns_svensson.pdf"),
       p1, width = 5.5, height = 3.2, device = cairo_pdf)

cat(sprintf("fig1: NS RMSE = %.1f bps, Svensson RMSE = %.1f bps\n",
            10000 * sqrt(mean(ns$residuals^2)),
            10000 * sqrt(mean(sv$residuals^2))))

# -----------------------------------------------------------------------------
# Figure 2: forward rate curve from the most-recent Nelson-Siegel fit.
# -----------------------------------------------------------------------------
grid2 <- seq(0.25, 30, length.out = 200)
spot <- 100 * yc_predict(ns, maturities = grid2)$rate
fwd  <- 100 * yc_forward(ns, maturities = grid2)$forward_rate

df2 <- rbind(
  data.frame(maturity = grid2, rate = spot, series = "Spot (zero) rate"),
  data.frame(maturity = grid2, rate = fwd,  series = "Instantaneous forward")
)
df2$series <- factor(df2$series,
                     levels = c("Spot (zero) rate", "Instantaneous forward"))

p2 <- ggplot(df2, aes(x = maturity, y = rate,
                       colour = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Spot (zero) rate" = ok_blue,
                                 "Instantaneous forward" = ok_red)) +
  scale_linetype_manual(values = c("Spot (zero) rate" = "solid",
                                   "Instantaneous forward" = "longdash")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30),
                     expand = c(0, 0.2)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Maturity (years)", y = "Rate") +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig2_forward.pdf"),
       p2, width = 5.5, height = 3.2, device = cairo_pdf)

# -----------------------------------------------------------------------------
# Figure 3: PCA loadings on a five-year panel of US Treasury curves.
# -----------------------------------------------------------------------------
pca_start <- as.Date("2019-01-01")
pca_end   <- dates_m[length(dates_m)]
wide_pca  <- wide_m[dates_m >= pca_start & dates_m <= pca_end, , drop = FALSE]

pc <- yc_pca(wide_pca, n_components = 3, scale = FALSE)

# pc$loadings should be a matrix of shape (maturities, components)
load_df <- NULL
if (!is.null(pc$loadings)) {
  lmat <- as.matrix(pc$loadings)
  # Flip signs so first PC is positive (conventional orientation).
  for (j in seq_len(ncol(lmat))) {
    if (mean(lmat[, j]) < 0) lmat[, j] <- -lmat[, j]
  }
  load_df <- data.frame(
    maturity = as.numeric(rownames(lmat)),
    PC1 = lmat[, 1], PC2 = lmat[, 2], PC3 = lmat[, 3]
  )
}

# Variance explained for legend labels.
ve <- if (!is.null(pc$variance_explained)) pc$variance_explained else pc$prop_var

lab1 <- sprintf("Level (PC1, %.0f%%)", 100 * ve[1])
lab2 <- sprintf("Slope (PC2, %.0f%%)", 100 * ve[2])
lab3 <- sprintf("Curvature (PC3, %.0f%%)", 100 * ve[3])

df3 <- rbind(
  data.frame(maturity = load_df$maturity, load = load_df$PC1, series = lab1),
  data.frame(maturity = load_df$maturity, load = load_df$PC2, series = lab2),
  data.frame(maturity = load_df$maturity, load = load_df$PC3, series = lab3)
)
df3$series <- factor(df3$series, levels = c(lab1, lab2, lab3))

p3 <- ggplot(df3, aes(x = maturity, y = load,
                       colour = series, linetype = series)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_colour_manual(values = setNames(c(ok_blue, ok_red, ok_green),
                                        c(lab1, lab2, lab3))) +
  scale_linetype_manual(values = setNames(c("solid", "longdash", "dotted"),
                                           c(lab1, lab2, lab3))) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 30)) +
  labs(x = "Maturity (years, log scale)",
       y = "Principal component loading") +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig3_pca.pdf"),
       p3, width = 5.5, height = 3.2, device = cairo_pdf)

cat(sprintf("fig3: PC1 = %.1f%%, PC2 = %.1f%%, PC3 = %.1f%%, cum = %.1f%%\n",
            100 * ve[1], 100 * ve[2], 100 * ve[3],
            100 * sum(ve[1:3])))

# -----------------------------------------------------------------------------
# Figure 4: duration vs maturity for zero-coupon and coupon (5%) bonds.
# -----------------------------------------------------------------------------
mats4 <- seq(1, 30, by = 1)

# Zero-coupon: Macaulay duration = maturity; modified duration = m / (1+y).
# Use yc_duration on the fitted NS curve for zero-coupon.
dur_zero <- yc_duration(ns, maturities = mats4,
                         compounding = "continuous")

# Coupon bond at each maturity: duration at the curve's yield at that
# maturity. Use yc_bond_duration scalar-wise.
rates_at <- yc_predict(ns, maturities = mats4)$rate
dur_coupon <- data.frame(
  maturity = mats4,
  duration = vapply(seq_along(mats4), function(i) {
    bd <- yc_bond_duration(coupon_rate = 0.05,
                           maturity = mats4[i],
                           yield = rates_at[i],
                           frequency = 2,
                           compounding = "semi_annual")
    bd$modified_duration
  }, numeric(1))
)

df4 <- rbind(
  data.frame(maturity = dur_zero$maturity,
             duration = dur_zero$modified_duration,
             series = "Zero-coupon"),
  data.frame(maturity = dur_coupon$maturity,
             duration = dur_coupon$duration,
             series = "5% coupon, semi-annual")
)
df4$series <- factor(df4$series,
                     levels = c("Zero-coupon", "5% coupon, semi-annual"))

p4 <- ggplot(df4, aes(x = maturity, y = duration,
                       colour = series, linetype = series)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
              colour = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Zero-coupon" = ok_blue,
                                 "5% coupon, semi-annual" = ok_red)) +
  scale_linetype_manual(values = c("Zero-coupon" = "solid",
                                   "5% coupon, semi-annual" = "longdash")) +
  scale_x_continuous(breaks = seq(0, 30, 5), expand = c(0, 0.3)) +
  scale_y_continuous(expand = c(0, 0.3)) +
  labs(x = "Maturity (years)", y = "Modified duration (years)") +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig4_duration.pdf"),
       p4, width = 5.5, height = 3.2, device = cairo_pdf)

# -----------------------------------------------------------------------------
# Figure 5: carry and roll-down at 3-month horizon across maturities.
# -----------------------------------------------------------------------------
mats5 <- c(2, 3, 5, 7, 10, 20, 30)
carry_df <- yc_carry(ns, maturities = mats5, horizon = 0.25)

# carry_df should include carry, roll_down, total in decimal units.
df5 <- data.frame(
  maturity = factor(carry_df$maturity, levels = rev(mats5)),
  Carry    = 10000 * carry_df$carry,
  Rolldown = 10000 * carry_df$rolldown
)
df5_long <- rbind(
  data.frame(maturity = df5$maturity, value = df5$Carry, component = "Carry"),
  data.frame(maturity = df5$maturity, value = df5$Rolldown, component = "Roll-down")
)
df5_long$component <- factor(df5_long$component, levels = c("Carry", "Roll-down"))

p5 <- ggplot(df5_long, aes(x = maturity, y = value, fill = component)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  scale_fill_manual(values = c("Carry" = ok_blue, "Roll-down" = ok_red)) +
  scale_y_continuous(labels = function(x) paste0(x, " bps"),
                     expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(x = "Maturity (years)",
       y = "Expected return over 3-month horizon") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig5_carry.pdf"),
       p5, width = 5.5, height = 3.2, device = cairo_pdf)

# -----------------------------------------------------------------------------
# Figure 6: the pandemic-era US Treasury curve at three dates.
# -----------------------------------------------------------------------------
date_labels <- c(
  format(date_pre, "%B %Y"),
  format(date_peak, "%B %Y"),
  format(date_now, "%B %Y")
)

df6 <- rbind(
  data.frame(maturity = s_pre$maturities,
             rate = 100 * s_pre$rates, date = date_labels[1]),
  data.frame(maturity = s_peak$maturities,
             rate = 100 * s_peak$rates, date = date_labels[2]),
  data.frame(maturity = s_now$maturities,
             rate = 100 * s_now$rates, date = date_labels[3])
)
df6$date <- factor(df6$date, levels = date_labels)

# Add fitted smooth curves for each date.
fits6 <- list(s_pre, s_peak, s_now)
df6_fit <- do.call(rbind, lapply(seq_along(fits6), function(i) {
  s <- fits6[[i]]
  fit <- yc_nelson_siegel(s$maturities, s$rates)
  grid <- seq(0.25, 30, length.out = 200)
  data.frame(
    maturity = grid,
    rate = 100 * yc_predict(fit, maturities = grid)$rate,
    date = date_labels[i]
  )
}))
df6_fit$date <- factor(df6_fit$date, levels = date_labels)

p6 <- ggplot() +
  geom_line(data = df6_fit,
            aes(x = maturity, y = rate, colour = date, linetype = date),
            linewidth = 0.8) +
  geom_point(data = df6, aes(x = maturity, y = rate, colour = date),
             size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = setNames(c(ok_green, ok_red, ok_blue),
                                        date_labels)) +
  scale_linetype_manual(values = setNames(c("dotted", "longdash", "solid"),
                                           date_labels)) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 30)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Maturity (years, log scale)", y = "Yield") +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 0.8)),
         linetype = guide_legend(nrow = 1)) +
  theme_wp(base_size = 10)

ggsave(file.path(fig_dir, "fig6_curves.pdf"),
       p6, width = 5.5, height = 3.2, device = cairo_pdf)

cat(sprintf("fig6: dates = %s, %s, %s\n",
            date_labels[1], date_labels[2], date_labels[3]))

# -----------------------------------------------------------------------------
# Table: parameter estimates from NS and Svensson at the most recent date.
# -----------------------------------------------------------------------------
tab_lines <- c(
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Parameter & Nelson-Siegel & Svensson \\\\",
  "\\midrule"
)
ns_p <- ns$params; sv_p <- sv$params
param_names <- c("beta0", "beta1", "beta2", "beta3", "tau1", "tau2")
for (pn in param_names) {
  ns_val <- if (!is.null(ns_p[[pn]])) sprintf("%.4f", ns_p[[pn]]) else "---"
  sv_val <- if (!is.null(sv_p[[pn]])) sprintf("%.4f", sv_p[[pn]]) else "---"
  tab_lines <- c(tab_lines,
    sprintf("$\\%s$ & %s & %s \\\\", pn, ns_val, sv_val))
}
tab_lines <- c(tab_lines, "\\midrule")
tab_lines <- c(tab_lines,
  sprintf("RMSE (bps) & %.1f & %.1f \\\\",
          10000 * sqrt(mean(ns$residuals^2)),
          10000 * sqrt(mean(sv$residuals^2))))
tab_lines <- c(tab_lines, "\\bottomrule", "\\end{tabular}")
writeLines(tab_lines, file.path(tab_dir, "params.tex"))

cat("\n--- done ---\n")
