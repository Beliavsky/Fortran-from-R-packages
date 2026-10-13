# SPDX-License-Identifier: MIT
# Run from the package root with the authentic CRAN credule 0.1.4 installed.
# Fixtures are frozen in test/reference_calls.inc; R is not a build dependency.
stopifnot(as.character(packageVersion("credule")) == "0.1.4")
library(credule)
out <- file("test/reference_calls.inc", "wt")
on.exit(close(out))
writeLines("! Generated with R credule 0.1.4; see test/reference/generate_reference.R.", out)
real_literal <- function(x) paste0(sprintf("%.17e", x), "_dp")
vector_lines <- function(x) {
  s <- real_literal(x)
  groups <- split(s, ceiling(seq_along(s) / 3))
  lines <- vapply(groups, paste, collapse = ", ", FUN.VALUE = "")
  lines[1] <- paste0("[", lines[1])
  lines[length(lines)] <- paste0(lines[length(lines)], "]")
  paste0("      ", lines, ifelse(seq_along(lines) < length(lines), ", &", ", &"))
}
emit <- function(name, arrays, recovery, pf, df, accrued, label) {
  writeLines(paste0("   call ", name, "( &"), out)
  for (x in arrays) writeLines(vector_lines(x), out)
  writeLines(paste0("      ", real_literal(recovery), ", ", pf, ", ", df, ", ",
                    if (accrued) ".true." else ".false.", ", '", label, "')"), out)
}
yt <- c(.25, .5, 1, 2, 5, 10, 20)
yields <- list(rep(0,7), rep(.04,7), c(-.02,-.015,-.005,.005,.015,.022,.03),
               c(.09,.08,.07,.05,.035,.02,.015))
ct <- c(.5, 1.5, 3, 5, 10)
cs <- exp(-cumsum(diff(c(0,ct)) * c(.004,.013,.007,.045,.028)))
price_count <- 0L
for (y in seq_along(yields)) for (pf in c(1L,2L,4L,12L)) {
 for (accrued in c(FALSE,TRUE)) for (recovery in c(0,.4,.8)) {
  df <- switch(as.character(pf), `1`=4L, `2`=12L, `4`=52L, `12`=365L)
  mt <- c(.25,.5,1,1.3,2.25,5,9,10,12)
  mt <- mt[mt >= 1/pf]
  expected <- priceCDS(yt,yields[[y]],ct,cs,mt,recovery,pf,df,accrued)$spread
  price_count <- price_count + 1L
  emit("check_price", list(yt,yields[[y]],ct,cs,mt,expected), recovery,pf,df,accrued,
       paste0("R price ",price_count))
 }
}
# Single credit and yield knots, unsorted requested maturities, sharp default risk.
for (hazard in c(.000001,.5,2)) for (accrued in c(FALSE,TRUE)) {
 mt <- c(5,1,.25,7,2.6)
 expected <- priceCDS(1,-.01,5,exp(-hazard*5),mt,.25,4,12,accrued)$spread
 price_count <- price_count + 1L
 emit("check_price",list(1,-.01,5,exp(-hazard*5),mt,expected),.25,4L,12L,accrued,
      paste0("R price ",price_count))
}
boot_count <- 0L
maturities <- c(1,2,3,5,10)
for (pf in c(1L,2L,4L,12L)) for (accrued in c(FALSE,TRUE)) {
 for (recovery in c(0,.4,.8)) for (curve in 1:3) {
  df <- if (pf == 12L) 52L else 12L
  hazards <- list(c(.005,.01,.02,.03,.04), c(.08,.04,.02,.012,.004), c(.3,.4,.1,.2,.1))[[curve]]
  sp <- exp(-cumsum(diff(c(0,maturities))*hazards))
  yr <- yields[[curve]]
  quotes <- priceCDS(yt,yr,maturities,sp,maturities,recovery,pf,df,accrued)$spread
  ref <- bootstrapCDS(yt,yr,maturities,quotes,recovery,pf,df,accrued)
  boot_count <- boot_count + 1L
  emit("check_bootstrap",list(yt,yr,maturities,quotes,ref$survprob,ref$hazrate),
       recovery,pf,df,accrued,paste0("R bootstrap ",boot_count))
 }
}
# Non-grid knots intentionally exercise the original separate leg anchors.
for (accrued in c(FALSE,TRUE)) for (pf in c(2L,4L,12L)) {
 mt <- c(.6,1.4,2.7,5.3)
 quotes <- c(.005,.011,.015,.02)
 ref <- bootstrapCDS(yt,yields[[4]],mt,quotes,.4,pf,12,accrued)
 boot_count <- boot_count + 1L
 emit("check_bootstrap",list(yt,yields[[4]],mt,quotes,ref$survprob,ref$hazrate),
      .4,pf,12L,accrued,paste0("R irregular bootstrap ",boot_count))
}
close(out)
cat("Frozen",price_count,"pricing cases and",boot_count,"bootstrap cases\n")
writeLines(c(paste("Pricing cases:",price_count),paste("Bootstrap cases:",boot_count),
             capture.output(sessionInfo())), "test/reference/R-session-info.txt")
