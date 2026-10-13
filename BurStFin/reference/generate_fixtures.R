# Regenerate with BurStFin 1.3 installed. Generated values are test data, not dependency code.
library(BurStFin)
stopifnot(as.character(packageVersion('BurStFin')) == '1.3')
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'test/r_reference.inc'
con <- file(out, 'w')
on.exit(close(con))
writeLines(c('! Actual R 4.5.0 / BurStFin 1.3 fixture values.',
             '! Regenerate with reference/generate_fixtures.R; see reference/README.md.'), con)
emit <- function(name, x) {
 d <- dim(x)
 if(is.null(d)) d <- length(x)
 vals <- formatC(as.numeric(x),digits=17,format='e')
 stopifnot(all(is.finite(x)))
 writeLines(paste0('real(dp), parameter :: ',name,'(',paste(d,collapse=','),') = reshape([ &'),con)
 for(i in seq_along(vals)) {
   end <- if(i == length(vals)) ' &' else ', &'
   writeLines(paste0('   ', vals[i], '_dp', end), con)
 }
 writeLines(paste0('   ], [',paste(d,collapse=','),'])'),con)
}
x <- outer(1:9,1:4,function(i,j) sin(i*(j+1)*.41)*.017+cos((i+2)*j*.31)*.008+i*j*.00015)
colnames(x) <- letters[1:4]
emit('x_reference',x)
m <- factor.model.stat(x, output='factor', verbose=0)
emit('factor_default_cov',fitted(m))
emit('factor_default_systematic',fitted(m,output='systematic'))
emit('factor_default_specific',fitted(m,output='specific'))
emit('factor_default_unique',m$uniquenesses)
emit('factor_default_sd',m$sdev)
emit('factor_default_cumulative',m$cumulative.variance.fraction)
m <- factor.model.stat(x, weights=1, center=FALSE, frac.var=.8, iter.max=0, range.factors=c(1,2),
                       specific.floor=.25, floor.type='fraction', output='factor', verbose=0)
emit('factor_noiter_cov',fitted(m))
s <- var.shrink.eqcor(x,verbose=0)
emit('shrink_default_cov',s)
emit('shrink_default_intensity',attr(s,'shrink'))
s <- var.shrink.eqcor(x,weights=c(0,1,2,3,0,2,1,3,2),shrink=.37,center=c(.002,-.003,.001,.005),
                     vol.shrink=.65,compatible=TRUE,tol=.03,verbose=0)
emit('shrink_explicit_cov',s)
heterogeneous_scales <- c(1e90,1e-90,1e40,1e-40)
s <- var.shrink.eqcor(sweep(x,2,heterogeneous_scales,'*'),shrink=.5,tol=0,sd.min=2,verbose=0)
emit('shrink_heterogeneous_scaled',s / outer(heterogeneous_scales,heterogeneous_scales))
y <- outer(1:15,1:6,function(i,j) sin(i*(j+2)*.23)*.02+cos((i+1)*(j+1)*.37)*.007+i*j*.00002)
colnames(y) <- letters[1:6]
emit('missing_base',y)
y[15,] <- NA
y[2,4] <- NA
y[c(1,3),5] <- NA
y[3:14,6] <- NA
m <- factor.model.stat(y, weights=seq(.2,1.6,length.out=15), output='factor', full.min=2, reg.min=3,
                       sd.min=4, nfac.miss=c(0,1,2), iter.max=3, range.factors=c(2,2),verbose=0)
emit('factor_missing_cov',fitted(m))
emit('factor_missing_unique',m$uniquenesses)
emit('factor_missing_sd',m$sdev)
emit('factor_missing_cumulative',m$cumulative.variance.fraction)
y2 <- y
y2[4,5] <- NA
m <- factor.model.stat(y2, weights=1, center=rep(.001,6), zero.load=TRUE, full.min=2, reg.min=3,
                       sd.min=4, quan.sd=.6, nfac.miss=c(0,1,2), iter.max=8, tol=1e-9,
                       range.factors=c(2,2), output='factor', verbose=0)
emit('factor_multireg_cov',fitted(m))
emit('factor_multireg_unique',m$uniquenesses)
emit('factor_multireg_sd',m$sdev)
s <- var.shrink.eqcor(y, weights=1, sd.min=4, vol.shrink=.3, tol=.005, verbose=0)
emit('shrink_missing_cov',s)
emit('shrink_missing_intensity',attr(s,'shrink'))
s <- var.shrink.eqcor(y, weights=1, sd.min=4, shrink=.2, center=FALSE, tol=0, verbose=0)
emit('shrink_missing_fixed_cov',s)
a <- alpha.proxy(weight=c(.8,.2,.5),vol.man=.24,vol.bench=.15,vol.other=.18,
                 cor.man=c(.5,-.3),cor.bench=.2,plot.it=FALSE,transpose=TRUE)
emit('alpha_grid_expected',a$z)
emit('slide_fraction_expected',slideWeight(8,c(.1875,.6875)))
emit('slide_observations_expected',slideWeight(8,observations=c(9,2)))
v <- matrix(c(.04,.01,-.005,.01,.09,.008,-.005,.008,.0225),3)
dimnames(v) <- list(letters[1:3],letters[1:3])
emit('benchmark_input',v)
b <- var.add.benchmark(v,c(c=-2,a=3))
emit('benchmark_added_expected',b)
emit('benchmark_relative_expected',var.relative.benchmark(b,'benchmark'))
close(con)
cat(R.version.string,'\nBurStFin ',as.character(packageVersion('BurStFin')),'\n',sep='')
