# SPDX-License-Identifier: GPL-2.0-or-later
# Run from this package directory with Rscript tools/generate_table.R.
# Martin Maechler's qDiptab, unchanged numerical values from diptest 0.77-2.
tab <- readRDS("upstream/inst/extraData/qDiptab.rds")
con <- file("src/diptest_table.f90", "w")
on.exit <- function() close(con)
writeLines(c("! SPDX-License-Identifier: GPL-2.0-or-later",
             "! Martin Maechler's diptest 0.77-2 qDiptab: 1,000,001 uniform samples per row.",
             "! Generated without rounding loss by tools/generate_table.R; see NOTICE.md.",
             "module diptest_table", "   use r_kinds, only : dp", "   implicit none", "   private",
             "   public :: table_sizes, table_probabilities, table_quantiles",
             "   integer, parameter :: table_sizes(21) = [ &",
             paste0("      ", paste(dimnames(tab)$n, collapse=", "), "]"),
             "   real(dp), parameter :: table_probabilities(26) = [ &"), con)
real_lines <- function(x, end) {
    x <- sprintf("%.17e_dp", x)
    for (i in seq.int(1, length(x), by=3)) {
        last <- min(i+2, length(x))
        suffix <- if (last == length(x)) end else ", &"
        writeLines(paste0("      ", paste(x[i:last], collapse=", "), suffix), con)
    }
}
real_lines(as.numeric(dimnames(tab)$Pr), "]")
writeLines("   real(dp), parameter :: table_quantiles(21, 26) = reshape([ &", con)
real_lines(as.vector(tab), "], [21, 26])")
writeLines("end module diptest_table", con)
close(con)
