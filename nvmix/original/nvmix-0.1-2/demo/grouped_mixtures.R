## Demo "grouped_mixtures"  (09-12-2020)
## By Erik Hintz, Marius Hofert and Christiane Lemieux 

## Numerical experiments for 'pgnvmix()', 'dgnvmix()', 'corgvnmix()', 
## and 'fitgStudentcopula()' ############

## Load packages
library(nvmix)
library(RColorBrewer) 
library(copula)
library(Matrix)
library(qrmdata)
library(qrmtools)
library(rugarch)
library(quadprog)

doPDF <- FALSE

cat("Note: This demo may take up to 24 hours.")

#### 1. Estimating the df #################################################################

set.seed(271)
runs <- 15 # number of different runs (i.e., different settings of 'scale' and 'upper')
dims <- c(5, 20) # dimensions 
dof <- list(df = 0.5 + 1:5, df = seq(from = 1, to = 5.75, length.out = 20)) # dof-settings
## Algorithm parameters
estims  <- c("pgnvmix", "empirical") 
methods <- c("PRNG", "sobol") # methods
max.iter <- 11
max.fun.evals  <- seq(from = 1920, to = 960000, length.out = max.iter)
## Arrays to store absolute error and the estimate 
abserrs <- 
   array(0, dim = c(length(dims), length(estims), length(methods), 
                    length(max.fun.evals), runs), 
         dimnames = list(d = dims, estims = estims, 
                         method = methods, n = max.fun.evals, run = 1:runs))
ests <- 
   array(0, dim = c(length(dims), length(estims), length(methods), 
                    length(max.fun.evals), runs), 
         dimnames = list(d = dims, estims = estims, 
                         method = methods, n = max.fun.evals, run = 1:runs))

for(dim in seq_along(dims)){ # for each dimension 
   d <- dims[dim] # current dimension
   df <- dof[[dim]] # current dof parameters
   pb. <- txtProgressBar(max = runs, style = 3) # progress bar
   ## Sample 'runs' upper limits 
   up   <- matrix(runif(d * runs) * sqrt(d) * 3, ncol = runs)
   ## Sample 'runs' Wishart matrices 
   sc   <- rWishart(runs, d, diag(d))
   for(run in 1:runs){ # for each run
      ## Grab current 'upper' and 'scale' 
      scale <- as.matrix(nearPD(cov2cor(sc[,,run]))$mat)
      upper <- as.vector(up[, run])
      ## 1. Construct the empirical estimator (PRNG only)
      sample   <- nvmix:::rgStudent(max(max.fun.evals), df = df, scale = scale,
                                    groupings = 1:d)
      inLimits <- sapply(1:max(max.fun.evals), function(i) all(sample[i, ] <= upper))
      ests.emp <- sapply(1:max.iter, function(i) mean(inLimits[1:max.fun.evals[i]]))
      errs.emp  <- 3.5*sqrt(ests.emp*(1-ests.emp)/max.fun.evals) # 3.5 also used in 'pgnvmix()'
      ## Store
      abserrs[dim, "empirical", "PRNG", , run] <- errs.emp
      ests[dim, "empirical", "PRNG", , run] <- ests.emp
      ## 2. Use 'pgnvmix()'
      for(nmax in seq_along(max.fun.evals)){ # for each sample size
         for(met in 1:2){ # for each method 
            ## Call 'pgnvmix()' (here: pgStudent)
            res <- 
               pgStudent(upper, df = df, scale = scale, verbose = FALSE, 
                         groupings = 1:d,
                         control = list(fun.eval = c(2^7, max.fun.evals[nmax]),
                                        pnvmix.abstol = 0, method = methods[met], 
                                        max.iter.rqmc = 1e8))
            abserrs[dim, "pgnvmix", met, nmax, run] <- attr(res, "abs. error")
            ests[dim, "pgnvmix", met, nmax, run] <- as.numeric(res)
         }
      }
      setTxtProgressBar(pb., run) # update progress bar
   }
   close(pb.)
}

## Get *mean* absolute errors in each setting over all length(reps) repetitions
mean.abs.errors <- 
   array(0, dim = c(length(dims), length(estims), length(methods), 
                    length(max.fun.evals)), 
         dimnames = list(d = dims, estims = estims, method = methods, 
                         n = max.fun.evals))
for(i in seq_along(dims)){ 
   for(j in seq_along(estims)){
      for(k in seq_along(methods)){
         for(l in seq_along(max.fun.evals)){
            ## Mean absolute errors in dimension d[index.dim] and for qmix[index.qmix]
            mean.abs.errors[i, j, k, l] <- mean(abserrs[i, j, k, l, ])
         }
      }
   }
}


## Prepare plot 
pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
cols <- pal(3) 
opar <- par(no.readonly = TRUE)
par(pty = "s")
## Each plot has mean estimated errors as fct of n for 'g (MC)', 'g (Sobol'), 'Crude MC'
for(i in seq_along(dims)){
   mean.abs.errors_<- cbind(mean.abs.errors[i, "empirical", "PRNG", ],
                            mean.abs.errors[i, "pgnvmix", "PRNG", ],
                            mean.abs.errors[i, "pgnvmix", "sobol", ])
   nms <- c("Crude MC", "g (MC)", "g (Sobol')")
   colnames(mean.abs.errors_) <- nms
   ## Compute regression coefficients
   coeff <- apply(mean.abs.errors_, 2, function(y) lm(log(y) ~ log(max.fun.evals))$coeff[2])
   names(coeff) <- nms
   ## Plot
   fname <- paste0("fig_pgnvmix_abserr_d", dims[i], ".pdf")
   if(doPDF) pdf(file = fname, width = 6, height = 6)
   plot(NA, log = "xy", xlim = range(max.fun.evals), ylim = range(mean.abs.errors_),
        xlab = "Number of function evaluations", ylab = "Average error")
   lgnd <- character(3)
   lwds <- c(1, 1.1, 1.7) # lwds[k] for 'lty = k'
   for(k in 1:3) {
      lines(max.fun.evals, mean.abs.errors_[, k], col = cols[k], lty = k, lwd = lwds[k])
      lgnd[k] <- paste0(nms[k]," (",round(coeff[k], 2),")")
   }
   legend("bottomleft", bty = "n", lty = 1:3, col = cols, 
          lwd = lwds, legend = lgnd)
   mtext(paste0("d = ", dims[i]), side = 4)
   if(doPDF) dev.off()
}
par(opar)

#### 2. Estimating the (log-)density function ##################################

### 2.1 Helper functions #######################################################

#' Evaluate integrand of a grouped t distribution
#'
#' @param u vector of evaluation points in (0,1)
#' @param d dimension of the distribution
#' @param x input 'x' at which density is to be estimated
#' @param loc location vector
#' @param scale scale matrix
#' @param df vector of dof parameters
#' 
#' @param groupings d-vector specifying group structure
#' @param log logical if log h(u) is to be returned
#' @return length(u)-vector with h(u) evaluations
integrand_h_t <- function(u, d = 2, x = rep(1, d), loc = rep(0, d), scale = diag(d),
                          df = rep(1, d), groupings = 1:d, log = FALSE){
   ## Some checks and declarations
   stopifnot(all(u >= 0), all(u <= 1), length(x) == d, length(loc) == d, 
             all.equal(dim(scale), c(d, d)))
   
   current.n  <- length(u) 
   lrdet      <- log(det(scale))/2
   scale.inv  <- solve(scale) 
   numgroups  <- length(unique(groupings))
   ## Obtain the quantile functions of the mixing variables
   mix_list <- nvmix:::get_mix_(qmix = "inverse.gamma", groupings = groupings,
                                df = df, callingfun = "pnvmix")
   ## Obtain realizations of the mixing variables 
   mixings <- mix_list$mix_(u)
   if(!is.matrix(mixings)) mixings <- cbind(mixings) # can happen in ungrouped case
   ## Transform to  1/sqrt(W_j), j=1,..,*d* 
   mixings <- mixings[, groupings, drop = FALSE] # (current.n, d)
   rt.mix.i <- t(1/sqrt(mixings)) # (d, current.n)  with1/sqrt(W_j), j=1,..,d  
   Dix <- rt.mix.i * matrix(x, ncol = current.n, nrow = d, byrow = FALSE)
   mahasq <- .colSums(Dix * (scale.inv %*% Dix), m = d, n = current.n)
   lres <- -d/2*log(2*pi) - lrdet - .rowSums(log(mixings), m = current.n, n = d)/2 - 
      mahasq/2 
   ## Return
   if(log) lres else exp(lres) 
}


#' Estimate grouped Student density with QUADPACK
#' 
#' @param x (n, d) matrix of evaluation points 
#' @param groupings d-vector specifying group structure
#' @param df vector of dof parameters
#' @param loc location vector
#' @param scale scale matrix
#' @param log logical if the logarithmic density shall be returned
#' @return n-vector with (log)density values estimated via 'integrate()' 
dgStudent_quadpack <- function(x, groupings = rep(1, d), df, loc = rep(0, d), 
                               scale = diag(d), log = FALSE){
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # dimension 
   n <- nrow(x) # number of input 'x' 
   ## Call 'integrate' for each row in 'x' 
   res <- sapply(1:n, function(i) integrate(
      integrand_h_t, lower = 0, upper = 1, d = d, groupings = groupings, 
      df = df, x = x[i, ])$value)
   ## Return 
   if(log) log(res) else res 
}


### 2.2 Numerical examples in the ungrouped case ###############################

## Evaluate the density function of a bivariate t-distribution with 6 dof:
## - True value via 'dStudent()'
## - Estimated using QUADPACK (i.e., using R's 'integrate()')
## - Estimated using dgStudent() 

## Specify quantile function of IG(df/2, df/2) to force estimation below
qmix <- function(u, df) 1 / qgamma(1-u, shape = df/2, rate = df/2)
d <- 2 # dimension
df <- 6 # dof parameter
x <- matrix(rep(c(0, 5, 10, 25, 50), d), ncol = d) # evaluation points
df <- 6 # dof parameter
res <- matrix(NA, ncol = 6, nrow = nrow(x))
colnames(res) <- c("True density", "True log-density", 
                   "Estimated density (QUADPACK)", 
                   "Estimated log-density (QUADPACK)",
                   "Estimated density (dgnvmix)",
                   "Estimated log-density (dgnvmix)")
## For each row in 'x'
for(i in 1:nrow(x)){
   ## True values
   res[i, 1] <- dStudent(x[i, ], df = df) 
   res[i, 2] <- dStudent(x[i, ], df = df, log = TRUE)
   ## Estimated via quadpack
   res[i, 3] <- dgStudent_quadpack(x[i, ], df = df)
   res[i, 4] <- log(res[i, 3])
   ## Estimated via dgnvmix()
   set.seed(271) # reproducibility 
   res[i, 6] <- dgnvmix(x[i, ], groupings = rep(1, d), qmix = qmix, df = df, 
                        log = TRUE)
   res[i, 5] <- exp(res[i, 6])
}
res

### 2.3 Numerical examples in the grouped case #################################

dims.dgnvmix <- c(2, 10) # dimensions
df <- c(3, 6) # dof parameters
groupings <- list(1:2, c(rep(1, 5), rep(2, 5))) # groupings
opar <- par(no.readonly = TRUE)
par(pty = "s") # square plot 
for(k in seq_along(dims.dgnvmix)){
   d <- dims.dgnvmix[k]
   loc <- rep(0, d)
   scale <- diag(d)
   n <- 2500 # sample size 
   x <- rgStudent(n, loc = loc, groupings = groupings[[k]], df = df/3, 
                  scale = scale)
   m <- mahalanobis(x, center = loc, cov = scale) # for plotting
   ## Sort 'x' and 'm' according to the order in 'm' 
   order.m <- order(m)
   m <- m[order.m]
   x <- x[order.m, ]
   ## Call 'dgnvmix()' and 'dgStudent_quadpack()' 
   t1 <- system.time(dg.1 <- dgnvmix(
      x, qmix = "inverse.gamma", loc = loc, scale = scale, 
      groupings = groupings[[k]], df = df, log = TRUE))
   t2 <- system.time(dg.2 <- dgStudent_quadpack(
      x, loc = loc, scale = scale, groupings = groupings[[k]], 
      df = df, log = TRUE))
   ## Also compute ungrouped t density for the two dof parameters
   dt.1 <- dStudent(x, loc = loc, scale = scale, df = df[1], log = TRUE)
   dt.2 <- dStudent(x, loc = loc, scale = scale, df = df[2], log = TRUE)
   ## Prepare plot
   pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
   cols <- pal(4)
   fname <- paste0("fig_dgnvmix_d", d, ".pdf")
   if(doPDF) pdf(file = fname, width = 7, height = 7)
   plot(sqrt(m), dg.1, type = 'l', col = cols[1], log = 'x',
        xlab = expression(paste("Mahalanobis Distance ", sqrt(x^T~x))), ylab = "log-density")
   lines(sqrt(m), dg.2, col = cols[2], lty = 2)
   lines(sqrt(m), dt.1, col = cols[3], lty = 3)
   lines(sqrt(m), dt.2, col = cols[4], lty = 4)
   legend("bottomleft", c(paste0("Estimated with dgnvmix(): ", round(t1[1] , 0), " sec"), 
                          paste0("Estimated with integrate(): ", round(t2[1] , 0), " sec"),
                          "True density (df = 3)", "True density (df = 6)"),
          lty = 1:4, col = cols, bty = 'n')
   mtext(paste0("d = ", d), side = 4)
   if(doPDF) dev.off()
}
par(opar)



#### 3. Estimating Kendall's tau ##################################

qmix <- "inverse.gamma"
## 6 settings of pairwise dof
df <- matrix( c(4, 8,
                1, 2, 
                4, 20,
                1, 5, 
                4, Inf,
                1, Inf),ncol = 2, byrow = TRUE)

l.df <- nrow(df)
scale <- seq(from = 1e-7, to = 1, length.out = 99) # avoid 0 for rel.error below
set.seed(271) # for reproducibility 
## Compute kendall's tau for each setting of 'df' and for each 'scale'
kendalls <- sapply(seq_len(l.df), function(i) 
   corgnvmix(scale, qmix = qmix, method = "kendall", df = df[i, ]))
## Include the elliptical approximation (exact when df1 = df2) 
kendall_ell <- corgnvmix(scale, method = "kendall", ellip.kendall = TRUE)
## Compute relative difference to elliptical approximation 
reldiff <- abs(kendalls - kendall_ell) / kendall_ell 
## Plot
opar <- par(no.readonly = TRUE)
par(pty = "s")
lgnd <- character(l.df + 1)
lgnd[1] <- "elliptical (equal df)"
pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
cols <- pal(l.df+1)
lwds <- c(1, 1.3, 1.8, 1.6, 1.3, 1.5, 1.3) # lwd for lty = 'l'
if(doPDF) pdf(file = "fig_corgnvmix_kendall.pdf", width = 7, height = 7)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = expression(rho),
     ylab = "Kendall's tau")
lines(scale, kendall_ell, lty = 1)
for(i in 1:l.df){
   lines(scale, kendalls[, i], col = cols[i + 1], lty = i + 1, lwd = lwds[i+1])
   lgnd[i+1] <- paste0("df1 = ", df[i, 1], ", df2 = ", df[i, 2])
} 
legend("topleft", lgnd, col = cols[1:(l.df + 1)], lty = 1:(l.df + 1), 
       lwd = lwds, bty = 'n')
if(doPDF) dev.off()
if(doPDF) pdf(file = "fig_corgnvmix_kendall_reldiff.pdf", width = 7, height = 7)
plot(NA, xlim = c(0, 1), ylim = range(reldiff), xlab = expression(rho),
     ylab = "Relative difference wrt elliptical case")
for(i in 1:l.df) lines(scale, reldiff[, i], col = cols[i + 1], lty = i + 1, lwd = lwds[i+1])
legend("topleft", rev(lgnd[-1]), col = rev(cols[2:(l.df + 1)]), lty = rev(2:(l.df + 1)),
       lwd = rev(lwds[2:(l.df+1)]), bty = 'n')
if(doPDF) dev.off()
par(opar)


#### 4. Numerical experiment for 'fitgStudentcopula()' #########################

n_ <- seq(from = 250, to = 2500, by = 250)
num.rep <- 15
d <- 6 # dimension 
df <- c(1, 4, 7) # dof for each group 
groupings <- rep(1:3, each = 2) # 2 components in each group 
control <- list(control.optim = list(maxit = 1000))
set.seed(1)
scale <- cov2cor(rWishart(1, d, diag(d))[,,1]) # same scale for all experiments
num.rep <- 15 # number of repetitions for each sample size
fit.res <- array(, dim = c(length(n_), num.rep, length(df), 2),
                 dimnames = list(n = n_, rep = 1:num.rep, df = c("df1", "df2", "df3"),
                                 est = c("init", "MLE")))
pb <- txtProgressBar(min = 0, max = num.rep, style = 3)
## Perform experiment
for(j in 1:num.rep){
   set.seed(j) # to make results easily reproducible
   sample <- rgStudentcopula(max(n_), groupings = groupings, scale = scale, df = df)
   for(i in seq_along(n_)){
      fit <- fitgStudentcopula(u = sample[1:n_[i], ], groupings = groupings, 
                            scale = scale, control = control)
      fit.res[i, j, , "init"] <- fit$df.init
      fit.res[i, j, , "MLE"] <- fit$df
   }
   setTxtProgressBar(pb, j)
}
close(pb)


## Create plot (initial estimates)
pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
cols <- pal(3)
if(doPDF) pdf(file = "fig_fitgstudentcop_init.pdf", height = 7, width = 7)
plot(NA, xlim = range(n_), ylim = range(fit.res[,,, ]), xlab = "n", ylab = "Estimate")
for(i in 1:num.rep){
   for(j in 1:3){
      points(n_, fit.res[, i, j, "init"], col = cols[j], pch = j)
   } 
}
abline(h = df, col = cols) # horizontal line with true df 
mtext(side = 4,  "Initial estimates for a grouped t copula (d=6), 15 repetitions for each n")
legend("topright", c("df1", "df2", "df3"), col = cols, pch = 1:3, bty = 'n')
if(doPDF) dev.off()
## Create plot (MLEs)
if(doPDF) pdf(file = "fig_fitgstudentcop_MLE.pdf", height = 7, width = 7)
plot(NA, xlim = range(n_), ylim = range(fit.res[,,, ]), xlab = "n", ylab = "Estimate")
for(i in 1:num.rep){
   for(j in 1:3){
      points(n_, fit.res[, i, j, "MLE"], col = cols[j], pch = j)
   } 
}
abline(h = df, col = cols) # horizontal line with true df 
mtext(side = 4,  "MLEs for a grouped t copula (d=6), 15 repetitions for each n")
legend("topright", c("df1", "df2", "df3"), col = cols, pch = 1:3, bty = 'n')
if(doPDF) dev.off()

#### 5. Data analysis example (Example 2) ######################################

## Load data, assign groups, clean #############################################

data(DJ_const) # load data 
stocks_ <- dimnames(DJ_const)[[2]] # stocknames 
d <- length(stocks_)
## Assign business sector (https://en.wikipedia.org/wiki/Dow_Jones_Industrial_Average)
## "InfTech"   Information technology 
## "FinServ"   Financial services 
## "Aero"      Aerospace and defence 
## "ConMin"    Construction and mining (CAT = caterpillar only)
## "Petro"     Petroleum industry 
## "Chem"      Chemical industry (DD = DowDuPont only)
## "Broad"     Broadcasting and entertainment (DIS = Disney only)
## "Congl"     Conglomerate
## "Retail"    Retailing
## "Pharma"    Pharmaceutical industry
## "Food"      Food industry 
## "Apparel"   Apparel (NKE = NIKE only)
## "FMcons"    Fast moving consumer goods (PG = Proctor&Gamble only)
## "ManHealth" Managed health care (UNH = United Health group only)
## "Tele"      Telecommunication (VZ = Verizon only)
sector_ <- c("InfTech", "FinServ", "Aero", "ConMin", "InfTech", "Petro", "Chem",
             "Broad", "Congl", "FinServ", "Retail", "InfTech", "InfTech", 
             "Pharma", "FinServ", "Food", "Food", "Congl", "Pharma", "InfTech",
             "Apparel", "Pharma", "FMcons", "FinServ", "ManHealth", "Aero",
             "FinServ", "Tele", "Retail", "Petro")

table(sector_)
## Avoid groups of size 1:
## - NKE (Nike) from "Apparel" to "Retail"
## - DIS (Disney) from "Broad" to "other"
## - DD (DowDuPont) from "Chem" to "other"
## - CAT (caterpillar) from "ConMin" to "other" 
## - PG (Proctor & Gamble) from "FMcons" to "Retail" 
## - UNH (United Health group) from "ManHealth" to "other"
## - VZ (Verizon) from "Tele" to "other"

sector_[stocks_ == "NKE"] <- "Retail"
sector_[stocks_ == "DIS"] <- "Other"
sector_[stocks_ == "DD"]  <- "Other"
sector_[stocks_ == "CAT"] <- "Other"
sector_[stocks_ == "PG"]  <- "Retail"
sector_[stocks_ == "UNH"] <- "Other"
sector_[stocks_ == "VZ"]  <- "Other"
stopifnot(all(table(sector_) > 1)) # all sectors have at least two constituents 
sector_names_ <- unique(sector_)

## Create vector 'groupings' for fitgStudentcopula() 
groupings <- rep(NA, d)
## Convert 
for(i in seq_along(sector_names_))
   groupings[sector_ == sector_names_[i]] <- i
numgroups <- length(unique(groupings))
## Grab data, compute returns 
time <- c("2014-01-01", "2015-12-31") # time period
x <- DJ_const[paste0(time, collapse = "/"),] # data
stopifnot(all(!is.na(x))) # no NAs 
X <- -returns(x) # -log-returns
n <- nrow(X) # 503
d <- ncol(X) # 30
## Fit marginal ARMA(1,1)-GARCH(1,1) models 
uspec <- rep(list(ugarchspec(distribution.model = "std")), ncol(X)) # GARCH specs
fit.ARMA.GARCH <- fit_ARMA_GARCH(X, ugarchspec.list = uspec) # fit ARMA-GARCH
stopifnot(sapply(fit.ARMA.GARCH$error, is.null)) # NULL = no error
fits <- fit.ARMA.GARCH$fit # fitted models
resi <- lapply(fits, residuals, standardize = TRUE) # grab out standardized residuals
Z <- as.matrix(do.call(merge, resi)) # standardized residuals
stopifnot(is.matrix(Z), nrow(Z) == n, ncol(Z) == d) # fail-safe programming
colnames(Z) <- colnames(x)


## Fit grouped and ungrouped t copula to 'Z' ###################################

### Ungrouped 
fit.nogroup <- fitgStudentcopula(x = Z, groupings = rep(1, d))
print(fit.nogroup$df) # 6.257661

names_ <- c("Ungrouped t", "Initial estimates", "Grouped t (maxit = 1500)", 
            "Grouped t (maxit = 3500)")
n <- 32 
u <- 1 - seq(0.95, to = 0.9995, length.out = n) # small levels for tail probabilities 

### Grouped 
maxit <- seq(from = 500, to = 3500, by = 500) # 'maxit' for 'optim()' 
fit.df <- matrix(NA, ncol = numgroups, nrow = length(maxit)+1)
colnames(fit.df) <- sector_names_
rownames(fit.df) <- c(0, maxit)
df.bounds <-  c(0.5, 30)
pb. <- txtProgressBar(max = length(maxit), style = 3) # progress bar
## Note: initial parameter always the same 
for(i in seq_along(maxit)){
   set.seed(4) # reproducibility 
   fit.group <- 
      fitgStudentcopula(x = Z, groupings = groupings, df.bounds = df.bounds, 
                        control = list(control.optim = list(maxit = maxit[i])))
   if(i == 1){
      ## Also store 'df.init'
      fit.df[1, ] <- fit.group$df.init 
      fit.df[i+1, ] <- fit.group$df
   } else {
      fit.df[i+1, ] <- fit.group$df
   }
   setTxtProgressBar(pb., i) # update progress bar
}

## Compute log-likelihood at the reported estimates
u <- copula:::pobs(Z) 
## Same estimate as returned by 'fitgStudentcopula()':
P <- as.matrix(Matrix::nearPD(pcaPP::cor.fk(u) * pi/2)$mat) 
factor.inv <- solve(t(chol(P)))
n_ <- c(0, maxit)
ll <- matrix(NA, nrow = length(n_), ncol = 6)
for(i in 1:6){
   ll[, i] <- sapply(seq_along(n_), function(k){
      set.seed(i)
      sum(dgStudentcopula(u, groupings = groupings, df = fit.df[k, ], 
                          factor.inv = factor.inv, log = TRUE))
   })
}
## Compute tail probabilities 
names_ <- c("Ungrouped t", "Initial estimates", "Grouped t (maxit = 1500)", 
            "Grouped t (maxit = 3500)")
n <- 32 
u <- 1 - seq(0.95, to = 0.9995, length.out = n) # small levels 
tailprobs <- array(NA, dim = c(length(u), 4), dimnames = list(u = u, names_))
df. <- fit.df[c(1, 4, 8), ]
for(j in 1:4){ # approximately 2.5 hours 
   set.seed(3) 
   tailprobs[, j] <- if(j > 1){
      pgStudentcopula(matrix(rep(u, each = d), ncol = d, byrow = TRUE),  # matrix of upper limits
                      groupings = groupings, scale = fit.nogroup$scale, df = df.[j-1, ], control = list(pnvmix.abstol = 1e-7)) # higher accuracy as prob's are small
   } else {
      pStudentcopula(matrix(rep(u, each = d), ncol = d, byrow = TRUE), scale = fit.nogroup$scale,
                     df = fit.nogroup$df, control = list(pnvmix.abstol = 1e-7))
   }
}


## Plots
n_ <- c(0, maxit)
cols <- pal(length(n_)) # colors 
## Plot of the estimates 
if(doPDF) pdf(file = "fig_dj30_estimates.pdf", width = 6.5, height = 6.5)
plot(NA, xlim = range(n_), ylim = range(fit.df), ylab = "Estimate", xlab = "Maxit")
for(i in 1:9)
   lines(n_, fit.df[, i], col = cols[i], type = 'b')
legend(x=2700, y = 18, sector_names_, col = cols, lty = 1, bty = 'n')
if(doPDF) dev.off()
## Plot of the likelihoods
if(doPDF) pdf(file = "fig_dj30_ll.pdf", width = 6.5, height = 6.5)
plot(NA, xlim = range(n_), ylim = range(ll), ylab = "log-likelihood", xlab = "Maxit")
for(i in 1:6)
   lines(n_, ll[,i ], lty = i)
if(doPDF) dev.off()
## Plot of the likelihoods ("zoomed in")
if(doPDF) pdf(file = "fig_dj30_ll_zoomed.pdf", width = 6.5, height = 6.5)
plot(NA, xlim = range(n_), ylim = c(4000, 4180), ylab = "log-likelihood", xlab = "Maxit")
for(i in 1:6)
   lines(n_, ll[,i ], lty = i)
if(doPDF) dev.off()
## Plot of the copula C(u,...,u)
cols <- pal(4) # colors
if(doPDF) pdf(file = "fig_dj30_shortfallprob.pdf", width = 8, height = 8)
plot(NA, xlim = range(u), ylim = range(tailprobs), xlab = "u",
     #ylab = expression(P(X[1]<=q[u],...,X[30]<=q[u])),
     ylab = expression(C(u,...,u)), 
     log = "y")
for(j in 1:4){
   lines(u, tailprobs[, j], type = 'l', col = cols[j])
}
legend("bottomright", names_[c(2, 1, 4, 3)], col = cols[c(2, 1, 4, 3)], lty = rep(1, 4), box.lty = 0)
if(doPDF) dev.off() 



#### 6. Portfolio management example (Example 3) ######################################

set.seed(271)
## Grab data, compute returns 
time <- c("2013-01-01", "2014-12-31") # time period
x <- DJ_const[paste0(time, collapse = "/"),] # data
stopifnot(all(!is.na(x))) # no NAs 
X <- returns(x) # log-returns
n <- nrow(X) # 503
d <- ncol(X) # 30


##' Given 'mu_t', 'sigma_t', 'next_X', compute portfolio return of optimal portfolio
##'
##' @param mu_t d-vector (expected excess return next period)
##' @param sigma_t (d, d) matrix (VCOV matrix of return next period)
##' @param next_X d-vector (actual returns next period)
##' @return numeric, 'constrained' portfolio return of next period
##'         (constrained = no shortselling)
##' @author Erik Hintz 
next_returns <- function(mu_t, sigma_t, next_X){
   sigma_t_inv <- solve(sigma_t)
   ## Unconstrained weights 
   w.unconstr <- sigma_t_inv %*% mu_t
   w.unconstr <- w.unconstr / abs(sum(w.unconstr)) # relative weights
   ret.unconstr <- sum(next_X * w.unconstr)
   ## Shortselling-constrained weights
   ret.constr <- if(all(w.unconstr >= 0)) ret.unconstr else {
      qp.sol <- solve.QP(Dmat = sigma_t, dvec = mu_t, Amat = diag(d), bvec = rep(0, d))
      w.constr <- qp.sol$solution
      w.constr  <- w.constr  / abs(sum(w.constr)) # relative weights
      sum(next_X * w.constr)
   }
   ## Return
   ret.constr
}

M <- 250 # sampling window
T <- nrow(X) # final time 
## Create result object
models <- c("historical", "ungrouped t", "grouped t")
res.pfret <- matrix(NA, ncol = T-M, nrow = 3)
rownames(res.pfret) <- models 
## Specify time-series model (standard t residuals) 
uspec <- rep(list(ugarchspec(distribution.model = "std")), d) # GARCH specs
n.sample <- 1e4
t.error <- c() # to store indices where 'fit.ARMA.GARCH()' returns an error 
## For each point in time 
for(t in (M+1):T){
   X_t <- X[(t-M):(t-1), ] # relevant data 
   next_X <- as.vector(X[t, ]) # return next period 
   
   ### 1. Model-free (historical) ############################################
   res.pfret["historical", t - M] <- next_returns(mu_t = colMeans(X_t), sigma_t = cov(X_t),
                                                  next_X = next_X)
   
   ### 2. Model based ##########################################################
   
   ## Fit marginal ARMA(1,1)-GARCH(1,1) models 
   fit.ARMA.GARCH <- fit_ARMA_GARCH(X_t, ugarchspec.list = uspec) # fit ARMA-GARCH
   if(!all(sapply(fit.ARMA.GARCH$error, is.null))){
      ## At least one error 
      t.error <- c(t.error, t)
      next 
   }
   Z <- sapply(fit.ARMA.GARCH$fit, residuals, standardize = TRUE) # standardized residuals
   U <- pobs(Z) # pseudo-observations 
   
   ## Helper function with input copula sample 'U.'; returns model based estimates of
   ## 'mu_t' and 'sigma_t' 
   next_mu_sigma <- function(U.){
      nu <- sapply(1:d, function(j) fit.ARMA.GARCH$fit[[j]]@fit$coef["shape"]) 
      ## Realizations of standardized residuals
      Z. <- sapply(1:d, function(j) sqrt((nu[j]-2)/nu[j]) * qt(U.[,j], df = nu[j])) # Z
      ## Simulate from fitted model
      X. <- sapply(1:d, function(j)
         fitted(ugarchsim(fit.ARMA.GARCH$fit[[j]], n.sim = n.sample, m.sim = 1, 
                          startMethod = "sample",
                          rseed = 271, custom.dist = list(name = "sample",
                                                          distfit = Z.[,j, drop = FALSE]))))
      ## Return mu_t, sigma_t 
      list(mu_t = colMeans(X.), sigma_t = cov(X.))
   }
   
   ## 3.1 ungrouped t copula ######################################################
   
   ## Estimate parameter matrix
   P <- sin(pcaPP::cor.fk(U) * pi/2)
   P <- as.matrix(Matrix::nearPD(P)$mat) # ensure positive-definiteness
   
   ## Estimate dof (only in first iteration)
   if(t == (M+1)){
      fit.copula <- fitgStudentcopula(u = U)
      df.ungrouped <- fit.copula$df
   }
   
   ## Sample from the fitted copula
   U. <- rgStudentcopula(n.sample, groupings = rep(1, d), df = df.ungrouped,
                         scale = P)
   ## Estimate 'mu_t' and 'sigma_t' from simulated stock returns based on 'U.'
   mu_sig <- next_mu_sigma(U.)
   ## Compute return for next period
   res.pfret["ungrouped t", t - M] <- 
      next_returns(mu_t = mu_sig$mu_t, sigma_t = mu_sig$sigma_t, next_X = next_X)
   
   ## 3.2 Grouped t copula ######################################################
   
   ## Estimate dof (only in first iteration)
   if(t == (M+1)){
      set.seed(271)
      ## Only in the first iteration
      fit.copula <- fitgStudentcopula(u = U, groupings = groupings, 
                                      control = list(control.optim = list(maxit = 1000)))
      df.grouped <- fit.copula$df
   }
   
   ## Sample from the fitted copula
   U. <- rgStudentcopula(n.sample, groupings = groupings, df = df.grouped, scale = P)
   ## Estimate 'mu_t' and 'sigma_t' from simulated stock returns based on 'U.'
   mu_sig <- next_mu_sigma(U.)
   ## Compute return for next period 
   res.pfret["grouped t", t - M] <- 
      next_returns(mu_t = mu_sig$mu_t, sigma_t = mu_sig$sigma_t, next_X = next_X)
   cat("\n", t, "\n")
}
stopifnot(length(t.error) <= 0.075 * (T-M)) # error in at most 7.5% of the runs 

(end <- (Sys.time() - start))

## Estimate mean return, CER, SR from 'res.pfret' for the different models
res.table <- matrix(NA, ncol = 3, nrow = 3)
colnames(res.table) <- c("Mean return", "CER", "SR") 
rownames(res.table) <- c("historical", "ungrouped t", "grouped t")
res.table[, 1] <- sapply(1:3, function(i) mean(res.pfret[i, ], na.rm = TRUE))
res.table[, 2] <- sapply(1:3, function(i) mean(res.pfret[i, ], na.rm = TRUE) - 
                            sd(res.pfret[i, ], na.rm = TRUE)^2/2)
res.table[, 3] <- sapply(1:3, function(i) mean(res.pfret[i, ], na.rm = TRUE) / 
                            sd(res.pfret[i, ], na.rm = TRUE))
## Print table 
print(round(res.table*100, 3))

save.image("~/Documents/work/overleaf_repos/pro/gnvmix/MDPI_risk/demo_grouped.RData")
##### END DEMO #######