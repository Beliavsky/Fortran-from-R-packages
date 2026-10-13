# Run with the official, unmodified CRAN CreditRisk 0.1.7 installed.
# The trace records stats::optim return values without changing the optimization.
library(CreditRisk)
stopifnot(as.character(packageVersion("CreditRisk")) == "0.1.7")
options(digits = 17)
outdir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(outdir, "R-reference-session.txt"))
print(sessionInfo())
sink()
write_fixture <- function(name, value) {
  x <- as.matrix(value)
  write.table(x, file.path(outdir, paste0(name, ".txt")), row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}
t <- c(.5, 1, 2, 4)
r <- c(.015, .017, .022, .028)
sigma <- c(.19, .23, .21, .26)
write_fixture("normal", cum_normal_density(c(-8, -2, 0, .5, 3, 8)))
write_fixture("call", generalized_black_scholes("c", 100, 105, c(.25, 1, 3), .03, .01, .22))
write_fixture("put", generalized_black_scholes("p", 100, 105, c(.25, 1, 3), .03, .01, .22))
write_fixture("merton", Merton(.55, 1, .25, .02, t))
write_fixture("black_cox", BlackCox(.55, .4, 1, .25, .02, .04, t))
write_fixture("at1p", at1p(1, .65, .4, sigma, r, t))
write_fixture("sbtv", sbtv(1, c(.4,.7), c(.8,.2), .15, sigma, r, t))
write_fixture("cds_full", cds(t, c(.01,.03,.025,.045), r))
write_fixture("cds_simple", cds(t, c(.01,.03,.025,.045), r, R=.012, RR=.25, simplified=TRUE))
write_fixture("cds2", cds2(t, c(.25,.75,1.5,3,5), c(.25,2,5), c(.01,.025,.02),
                          c(.1,1.5,4.5), c(.015,.04,.05), R=.009, RR=.35))
# Real R seed stream, followed by its exact innovations for cross-language path equality.
write_fixture("simulation", Merton.sim(20, .05, .2, c(0,.25,1,2), 3, seed=1729))
set.seed(1729)
dt <- diff(c(0,c(0,.25,1,2)))
# R rnorm with zero sd does not consume a standard-normal draw.
innov <- matrix(rnorm(12, mean=(.05-.2^2/2)*dt, sd=.2*sqrt(dt)),4,3)
shocks <- sweep(sweep(innov,1,(.05-.2^2/2)*dt,"-"),1,.2*sqrt(dt),"/")
shocks[1,] <- 0
write_fixture("shocks", shocks)
trace(stats::optim, exit=quote(assign("last_optim", returnValue(), envir=.GlobalEnv)), print=FALSE)
q <- cds(t, BlackCox(.5,.425,1,.24,.02,.07,t)$Default.Intensity,.02)$Rate
write_fixture("quote_black_cox", q)
bc <- calibrate.BlackCox(1,q,.02,t)
write_fixture("calibrate_black_cox_reported", bc)
write_fixture("calibrate_black_cox_parameters", c(last_optim$par[1],last_optim$par[1]-last_optim$par[2],last_optim$par[3:4]))
write_fixture("calibrate_black_cox_error", last_optim$value)
write_fixture("calibrate_black_cox_rates", cds(t,BlackCox(last_optim$par[1],last_optim$par[1]-last_optim$par[2],
                1,last_optim$par[3],.02,last_optim$par[4],t)$Default.Intensity,.02)$Rate)
q <- cds(t,at1p(1,.42,.04,c(.21,.23,.22,.24),r,t)$Default.Intensity,r)$Rate
write_fixture("quote_at1p", q)
a <- calibrate.at1p(1,q,r,t)
write_fixture("calibrate_at1p_parameters", c(a$H0,a$B,a$Sigma))
write_fixture("calibrate_at1p_error", a$Error)
write_fixture("calibrate_at1p_rates", cds(t,at1p(1,a$H0,a$B,a$Sigma,r,t)$Default.Intensity,r)$Rate)
q <- cds(t,sbtv(1,c(.42,.58),c(.8,.2),.04,c(.21,.23,.22,.24),r,t)$Default.Intensity,r)$Rate
write_fixture("quote_sbtv", q)
s <- calibrate.sbtv(1,c(.8,.2),q,r,t)
write_fixture("calibrate_sbtv_parameters", c(s$B,s$Sigma,s$H))
write_fixture("calibrate_sbtv_error", s$Error)
write_fixture("calibrate_sbtv_rates", cds(t,sbtv(1,s$H,c(.8,.2),s$B,s$Sigma,r,t)$Default.Intensity,r)$Rate)
maturity <- c(.75,1.5,3)
c <- calibrate.cds(r,t,maturity,c(.008,.013,.02))
write_fixture("calibrate_cds_parameters", c$min_int)
write_fixture("calibrate_cds_error", c$error)
write_fixture("calibrate_cds_rates", approx(t,cds(t,approx(maturity,c$min_int,t,rule=2)$y,r)$Rate,maturity,rule=2)$y)
untrace(stats::optim)
# Independent-review regression: large trial barriers previously exposed a
# false convergence caused by an excessively large artificial invalid penalty.
t_stress <- c(.5,1.5,3,7)
r_stress <- c(.011,.015,.024,.035)
q_stress <- cds(t_stress,sbtv(1,c(.35,.75),c(.7,.3),.15,c(.22,.27,.31,.26),r_stress,t_stress)$Default.Intensity,r_stress)$Rate
s_stress <- calibrate.sbtv(1,c(.7,.3),q_stress,r_stress,t_stress)
write_fixture("sbtv_stress_quotes", q_stress)
write_fixture("sbtv_stress_parameters", c(s_stress$B,s_stress$Sigma,s_stress$H))
write_fixture("sbtv_stress_error", s_stress$Error)
write_fixture("sbtv_stress_rates", cds(t_stress,sbtv(1,s_stress$H,c(.7,.3),s_stress$B,s_stress$Sigma,r_stress,t_stress)$Default.Intensity,r_stress)$Rate)
