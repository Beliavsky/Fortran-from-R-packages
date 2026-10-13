# Reproduce test_metrics.f90 and test_model_weights.f90 references.
# Run from the package root with exact pinned R versions.
stopifnot(packageVersion("loo") == "2.10.1")
library(loo)
options(digits=17)
s<-6L
x<-matrix(seq_len(24)/4-2,s,4)
x2<-matrix(c(seq(4,1,length.out=s),seq(-2,2,length.out=s),seq(0,1,length.out=s),seq(2,-1,length.out=s)),s,4)
y<-c(.3,1.2,-.5,2.1)
ll<-matrix(-log(seq_len(s)),s,4)
set.seed(319)
pp<-replicate(2,sample(seq_len(s)))
cat('PERMUTATIONS\n')
dput(pp)
for(fun in c('crps','scrps','loo_crps','loo_scrps')) {
set.seed(319)
out<-if(startsWith(fun,'loo')) suppressWarnings(do.call(fun,list(x=x,x2=x2,y=y,log_lik=ll,permutations=2))) else do.call(fun,list(x=x,x2=x2,y=y,permutations=2))
cat(fun,'\n')
dput(out)
}
cat('MARGINAL\n')
dput(drop(weights(suppressWarnings(psis(-ll)),normalize=TRUE,log=FALSE)[,1]))
cat('JOINT\n')
for(p in 1:2) dput(drop(weights(suppressWarnings(psis(-ll-ll[pp[,p],])),normalize=TRUE,log=FALSE)[,1]))
for(metric in c('mae','mse','rmse')) {
cat(metric,'\n')
dput(suppressWarnings(loo_predictive_metric(x,y,ll,metric=metric)))
}
xc<-matrix(c(.1,.2,.3,.7,.8,.9,.9,.8,.7,.2,.1,.0,.0,.1,.2,.3,.4,.5,.9,.8,.7,.6,.5,.4),s,4)
for(metric in c('acc','balanced_acc')) {
cat(metric,'\n')
dput(suppressWarnings(loo_predictive_metric(xc,c(0,0,1,1),ll,metric=metric)))
}
lpd<-log(matrix(c(.8,.2,.1,.7,.3,.2,.1,.8,.6,.1,.8,.2,.1,.1,.7,.2,.1,.8),6,3))
cat('STACKING\n')
dput(unclass(stacking_weights(lpd,optim_control=list(reltol=1e-13,maxit=10000))))
cat('PSEUDOBMA\n')
dput(unclass(pseudobma_weights(lpd,BB=FALSE)))
set.seed(904)
bb<-loo:::dirichlet_rng(4,rep(1,6))
cat('BOOTSTRAP\n')
dput(bb)
set.seed(904)
cat('PSEUDOBMA PLUS\n')
dput(unclass(pseudobma_weights(lpd,BB=TRUE,BB_n=4)))
