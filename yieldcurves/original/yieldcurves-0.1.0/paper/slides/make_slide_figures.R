# make_slide_figures.R (yieldcurves)
# Prepare slide-ready figures:
#  1. Copy the paper's hero figure (US Treasury curves through the pandemic cycle)
#  2. Copy the four gallery figures (NS/Svensson, PCA, carry, duration)
#  3. Copy the deep-dive figure (PCA loadings, reused for slide 12)
#  4. Generate a QR code pointing to the paper PDF
#
# Usage:  Rscript make_slide_figures.R

suppressPackageStartupMessages({
  if (!requireNamespace("qrcode", quietly = TRUE)) {
    install.packages("qrcode", repos = "https://cloud.r-project.org")
  }
  library(qrcode)
})

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

paper_figs <- file.path("..", "figures")

# ------------------------------------------------------------------
# 1. Hero figure: US Treasury curves at three dates (fig6_curves.pdf)
# ------------------------------------------------------------------
hero_src <- file.path(paper_figs, "fig6_curves.pdf")
hero_dst <- file.path(fig_dir, "hero_figure.pdf")

if (file.exists(hero_src)) {
  file.copy(hero_src, hero_dst, overwrite = TRUE)
  cat("Copied hero figure to", hero_dst, "\n")
} else {
  stop("Source figure not found: ", hero_src,
       ". Run the paper's make_figures.R first.")
}

# ------------------------------------------------------------------
# 2. Gallery figures: NS/Svensson, PCA, carry, duration
# ------------------------------------------------------------------
gallery <- c(
  "fig1_ns_svensson.pdf",
  "fig3_pca.pdf",
  "fig4_duration.pdf",
  "fig5_carry.pdf"
)

for (f in gallery) {
  src <- file.path(paper_figs, f)
  dst <- file.path(fig_dir, f)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("Copied", f, "to", dst, "\n")
  } else {
    warning("Gallery figure missing: ", src)
  }
}

# ------------------------------------------------------------------
# 3. QR code to the paper PDF on the publications page
# ------------------------------------------------------------------
paper_url <- "https://charlescoverdale.github.io/files/coverdale_yieldcurves_2026.pdf"

qr <- qr_code(paper_url, ecl = "M")
png(
  filename = file.path(fig_dir, "qrcode_paper.png"),
  width = 800, height = 800, res = 300, bg = "white"
)
par(mar = rep(0, 4))
plot(qr)
dev.off()

cat("QR code written to", file.path(fig_dir, "qrcode_paper.png"), "\n")
