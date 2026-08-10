set.seed(7) # for reproducible results
test.join <- function() {
   checkEquals(nlsic::join(" ", c("Hello", "World"), "!"), "Hello World !")
}
test.Nulla <- function() {
   n=5
   v=rnorm(n)
   base=nlsic::Nulla(v)
   checkEqualsNumeric(t(v)%*%base, rep(0, n-1))
}
test.lsi <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   xe=rnorm(nc)
   b=a%*%xe
   # plain ls
   checkEqualsNumeric(nlsic::lsi(a, b), xe)
   # add inequality
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=nlsic::uplo2uco(par, upper=c(x1=xe[1L]-1))
   # inequality must be enforced
   x=nlsic::lsi(a, b, uco$u, uco$co)
   checkEqualsNumeric(x[1], xe[1]-1)
   # Kuhn-Tucker condition: t(a)%*%r can be decomposed in t(u) for enforced inequalities
   checkEqualsNumeric((t(a)%*%(a%*%x-b))[-1], rep(0, nc-1))
}
test.lsi_ln <- function() {
#browser()
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   a[,nc]=a[,nc-1] # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # plain nlsic::ls_ln
   x=nlsic::lsi_ln(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
   # add inequality
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=nlsic::uplo2uco(par, upper=c(x1=xe[1L]-1))
   # inequality must be enforced
   xu=nlsic::lsi_ln(a, b, uco$u, uco$co)
   checkEqualsNumeric(xu[1], xe[1]-1)
   # the last two must be equal
   checkEqualsNumeric(xu[nc], xu[nc-1])
   # Kuhn-Tucker condition: t(a)%*%r can be decomposed in t(u) for enforced inequalities
   checkEqualsNumeric((t(a)%*%(a%*%xu-b))[-1], rep(0, nc-1))
   # mnorm used when no inequality
   mnormid=diag(nrow=nc) # must give the same result as mnorm==NULL
   xmid=nlsic::lsi_ln(a, b, mnorm=mnormid)
   checkEqualsNumeric(xmid,x)
   mnorm=t(c(0., 2., 1.))
   xm=nlsic::lsi_ln(a, b, mnorm=mnorm)
   checkEqualsNumeric(mnorm%*%xm, 0.)
   # mnorm used with inequality
   xmuid=nlsic::lsi_ln(a, b, uco$u, uco$co, mnorm=mnormid)
   checkEqualsNumeric(xmuid, xu)
   xmu=nlsic::lsi_ln(a, b, uco$u, uco$co, mnorm=mnorm)
   checkEqualsNumeric(uco$u%*%xmu, uco$co)
   checkEqualsNumeric(mnorm%*%xmu, 0.)
   # hand made example, solution is (0,0,1)
   a=diag(nrow=3L)[1:2,]
   b=double(2L)
   u=t(rep(1., 3L))
   co=1.
   x=nlsic::lsi_ln(a, b, u, co)
   checkEqualsNumeric(x, c(0., 0., 1.))
   # based on previous example, solution is c(0.25, 0.25, 0.5)
   u2=rbind(u, c(0, 0, -1))
   co2=c(co, -0.5)
   x2=nlsic::lsi_ln(a, b, u2, co2)
   checkEqualsNumeric(x2, c(0.25, 0.25, 0.5))
}
test.ls_ln <- function() {
   nr=5L
   nc=3L
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   # full rank when b is a matrix
   xe=matrix(rnorm(nc*2L), nrow=nc, ncol=2L)
   b=a%*%xe
   x=nlsic::ls_ln(a, b)
   checkEqualsNumeric(x, xe)
   # rank deficient case when b is a matrix
   a[,nc]=a[,nc-1L] # make rank deficient
   b=a%*%xe
   # plain nlsic::ls_ln
   x=nlsic::ls_ln(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc,], x[nc-1L,])
   # rank deficient case for one b
   xe=rnorm(nc)
   b=a%*%xe
   # plain nlsic::ls_ln
   x=nlsic::ls_ln(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
}
test.ls_ln_svd <- function() {
   nr=5
   nc=3
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   a[,nc]=a[,nc-1] # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # plain nlsic::ls_ln
   x=nlsic::ls_ln_svd(a, b)
   # the last two must be equal
   checkEqualsNumeric(x[nc], x[nc-1])
   # must be consistent with QR solution
   checkEqualsNumeric(x, nlsic::ls_ln(a, b))
}
test.ldp <- function() {
   # prepare inequality x1 >= 1
   nc=3L
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=nlsic::uplo2uco(par, lower=c(x1=1))
   # solution must be c(1, 0, 0)
   x=nlsic::ldp(uco$u, uco$co)
   checkEqualsNumeric(x, c(1, 0, 0))
   # x1 >= -1
   uco=nlsic::uplo2uco(par, lower=c(x1=-1))
   # solution must be c(0, 0, 0)
   x=nlsic::ldp(uco$u, uco$co)
   checkEqualsNumeric(x, c(0, 0, 0))
   # unfeasible inequalities
   uco=nlsic::uplo2uco(par, upper=c(x1=-1), lower=c(x1=1))
   x=nlsic::ldp(uco$u, uco$co)
   checkTrue(is.null(x))
}
test.nlsic <- function() {
   # solve min_{a,b} ||exp(a*x+b)-meas||, a,b>=1
   a_true=1; b_true=2; x=0:5
   # simulation function
   sim=function(par, x) exp(par[["a"]]*x+par[["b"]])
   # residual function to be passed to nlsic()
   resid=function(par, cjac, ...) {
      dots=list(...)
      s=sim(par, dots$x)
      result=list(res=s-dots$meas)
      if (cjac) {
         result$jacobian=cbind(a=s*dots$x, b=s)
      }
      result
   }
   # simulated measurements for true parameters
   meas=sim(c(a=a_true, b=b_true), x)
   # starting values for par
   par=c(a=0, b=0)
   # prepare constraints
   uco=nlsic::uplo2uco(par, lower=c(a=1, b=1))
   # main call: solve the problem
   fit <- nlsic::nlsic(par, resid, uco$u, uco$co, x=x, con=list(history=TRUE), meas=meas)
   if (fit$error == 1) {
      stop(fit$mes)
   } else {
      checkEqualsNumeric(fit$par, c(a_true, b_true))
      if (! is.null(fit$mes)) {
         warning(fit$mes)
      }
   }
   checkTrue(any(fit$history$nb_reuse > 0L))
}
test.lsie_ln <- function() {
   nr=5
   nc=4
   a=matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
   isame=c(1L, 3L, 4L) # not just the last two, to be sure that permutation is enabled
   a[,isame[-1L]]=a[,isame[1L]]  # make rank deficient
   xe=rnorm(nc)
   b=a%*%xe
   # impose that x2 be xe[2]+0.5
   e=t(double(nc))
   e[1L,2L]=1.
   ce=xe[2]+0.5
   # plain nlsic::ls_ln
   xln=nlsic::lsie_ln(a, b)
   checkEqualsNumeric(xln, nlsic::ls_ln(a, b))
   # the free must be equal
   checkEqualsNumeric(rep(xln[isame[1L]], length(isame)-1L), xln[isame[-1L]])
   # add equality
   xec=nlsic::lsie_ln(a, b, e=e, ce=ce)
   checkEqualsNumeric(e%*%xec, ce) # it must be enforced
   # add inequality
   par=seq_len(nc)
   names(par)=paste("x", par, sep="")
   uco=nlsic::uplo2uco(par, upper=c(x1=xe[1L]-10))
   xu=nlsic::lsi_ln(a, b, uco$u, uco$co)
   checkEqualsNumeric(uco$u%*%xu, uco$co)   # inequality must be enforced
   xue=nlsic::lsie_ln(a, b, uco$u, uco$co, e=e, ce=ce)
   checkEqualsNumeric(uco$u%*%xue, uco$co)
   checkEqualsNumeric(e%*%xue, ce)
   # mnorm
   mnormid=diag(nrow=nc) # must give the same result as mnorm==NULL
   xmid=nlsic::lsie_ln(a, b, uco$u, uco$co, e=e, ce=ce, mnorm=mnormid)
   checkEqualsNumeric(xmid, xue)
   mnorm=t(c(0., 0., 1., 2.))
   xm=nlsic::lsie_ln(a, b, uco$u, uco$co, e=e, ce=ce, mnorm=mnorm)
   checkEqualsNumeric(uco$u%*%xm, uco$co)
   checkEqualsNumeric(e%*%xm, ce)
   checkEqualsNumeric(mnorm%*%xm, 0.)
}
