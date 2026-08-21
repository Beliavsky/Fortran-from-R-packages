! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_distributions
  use mc2d_kinds, only : dp, pi, nan_dp
  use mc2d_utils, only : clamp01, sort_dp, argsort_dp, is_close
  use mc2d_random, only : random_uniform_open, random_normal, random_gamma, random_beta, random_poisson, random_binomial
  use mc2d_special, only : beta_pdf, beta_cdf, beta_quantile, ncbeta_pdf, ncbeta_cdf, ncbeta_quantile, &
                           normal_cdf, normal_quantile
  use mvtnorm_distributions, only : dmvnorm_one, rmvnorm
  implicit none
  private

  public :: dbern, pbern, qbern, rbern
  public :: dbetagen, pbetagen, qbetagen, rbetagen
  public :: dbetasubj, pbetasubj, qbetasubj, rbetasubj
  public :: dlnormb, plnormb, qlnormb, rlnormb
  public :: dpert, ppert, qpert, rpert, dpert_mean, ppert_mean, qpert_mean, rpert_mean
  public :: dtriang, ptriang, qtriang, rtriang
  public :: dtriang_mean, ptriang_mean, qtriang_mean, rtriang_mean
  public :: dmqi, pmqi, qmqi, rmqi
  public :: dempiricald, pempiricald, qempiricald, rempiricald
  public :: dempiricalc, pempiricalc, qempiricalc, rempiricalc
  public :: ddirichlet, rdirichlet, dmultinomial, rmultinomial
  public :: dmultinormal, rmultinormal, dmultinormal_varying, rmultinormal_varying
  public :: rtrunc, lhs

  abstract interface
    function scalar_cdf(q) result(p)
      import dp
      real(dp), intent(in) :: q
      real(dp) :: p
    end function scalar_cdf
    function scalar_quantile(p) result(q)
      import dp
      real(dp), intent(in) :: p
      real(dp) :: q
    end function scalar_quantile
  end interface

contains

  real(dp) function dbern(x,prob,log_density) result(d)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::prob
    logical,intent(in),optional::log_density
    real(dp)::pr
    logical::llog
    pr=0.5_dp; if(present(prob)) pr=prob
    llog=.false.; if(present(log_density)) llog=log_density
    if(pr<0.0_dp .or. pr>1.0_dp) then; d=nan_dp(); return; end if
    if(is_close(x,0.0_dp,0.0_dp)) then; d=1.0_dp-pr
    else if(is_close(x,1.0_dp,0.0_dp)) then; d=pr
    else; d=0.0_dp; end if
    if(llog) then; if(d>0) d=log(d); if(d==0) d=-huge(1.0_dp); end if
  end function dbern

  real(dp) function pbern(q,prob,lower_tail,log_p) result(p)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::prob
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pr
    logical::lt,lp
    pr=0.5_dp; if(present(prob)) pr=prob
    lt=.true.; if(present(lower_tail)) lt=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    if(q<0) then; p=0; else if(q<1) then; p=1-pr; else; p=1; end if
    if(.not.lt) p=1-p
    if(lp) then; if(p>0) p=log(p); if(p==0) p=-huge(1.0_dp); end if
  end function pbern

  real(dp) function qbern(p,prob,lower_tail,log_p) result(q)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::prob
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp,pr
    logical::lt,lp
    pr=0.5_dp; if(present(prob)) pr=prob
    lt=.true.; if(present(lower_tail)) lt=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    pp=p; if(lp) pp=exp(pp); if(.not.lt) pp=1-pp
    if(pp<0 .or. pp>1 .or. pr<0 .or. pr>1) then; q=nan_dp(); else if(pp<=1-pr) then; q=0; else; q=1; end if
  end function qbern

  integer function rbern(prob) result(x)
    real(dp),intent(in),optional::prob
    real(dp)::pr,u
    pr=0.5_dp; if(present(prob)) pr=prob; call random_number(u); x=merge(1,0,u<pr)
  end function rbern

  real(dp) function dbetagen(x,shape1,shape2,xmin,xmax,ncp,log_density) result(d)
    real(dp),intent(in)::x,shape1,shape2
    real(dp),intent(in),optional::xmin,xmax,ncp
    logical,intent(in),optional::log_density
    real(dp)::a,b,nc,z,lo,hi
    logical::llog
    lo=0; hi=1; nc=0; if(present(xmin)) lo=xmin; if(present(xmax)) hi=xmax; if(present(ncp)) nc=ncp
    llog=.false.; if(present(log_density)) llog=log_density
    if(hi<=lo) then; d=nan_dp(); return; end if
    if(x<lo .or. x>hi) then; d=0; if(llog) d=-huge(1.0_dp); return; end if
    z=(x-lo)/(hi-lo); a=shape1; b=shape2
    if(nc==0) then; d=beta_pdf(z,a,b,llog); else; d=ncbeta_pdf(z,a,b,nc,llog); end if
    if(llog) then; d=d-log(hi-lo); else; d=d/(hi-lo); end if
  end function dbetagen

  real(dp) function pbetagen(q,shape1,shape2,xmin,xmax,ncp,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,shape1,shape2
    real(dp),intent(in),optional::xmin,xmax,ncp
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::lo,hi,nc,z
    logical::lt,lp
    lo=0; hi=1; nc=0; if(present(xmin)) lo=xmin; if(present(xmax)) hi=xmax; if(present(ncp)) nc=ncp
    lt=.true.; if(present(lower_tail)) lt=lower_tail; lp=.false.; if(present(log_p)) lp=log_p
    if(hi<lo) then; p=nan_dp(); return; end if
    if(hi==lo .and. is_close(q,lo)) then; p=merge(1.0_dp,0.0_dp,lt)
    else if(hi==lo) then; p=nan_dp(); return
    else
      z=(q-lo)/(hi-lo); if(nc==0) then; p=beta_cdf(z,shape1,shape2); else; p=ncbeta_cdf(z,shape1,shape2,nc); end if
      if(.not.lt) p=1-p
    end if
    if(lp) then; if(p>0) p=log(p); if(p==0) p=-huge(1.0_dp); end if
  end function pbetagen

  real(dp) function qbetagen(p,shape1,shape2,xmin,xmax,ncp,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,shape1,shape2
    real(dp),intent(in),optional::xmin,xmax,ncp
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::lo,hi,nc,pp,z
    logical::lt,lp
    lo=0; hi=1; nc=0; if(present(xmin)) lo=xmin; if(present(xmax)) hi=xmax; if(present(ncp)) nc=ncp
    lt=.true.; if(present(lower_tail)) lt=lower_tail; lp=.false.; if(present(log_p)) lp=log_p
    pp=p; if(lp) pp=exp(pp); if(.not.lt) pp=1-pp
    if(hi<lo) then; q=nan_dp(); return; end if
    if(nc==0) then; z=beta_quantile(pp,shape1,shape2); else; z=ncbeta_quantile(pp,shape1,shape2,nc); end if
    q=lo+z*(hi-lo)
  end function qbetagen

  real(dp) function rbetagen(shape1,shape2,xmin,xmax,ncp) result(r)
    real(dp),intent(in)::shape1,shape2
    real(dp),intent(in),optional::xmin,xmax,ncp
    real(dp)::lo,hi,nc,z
    integer::k
    lo=0; hi=1; nc=0; if(present(xmin)) lo=xmin; if(present(xmax)) hi=xmax; if(present(ncp)) nc=ncp
    if(hi<lo) then; r=nan_dp(); return; end if
    if(nc==0) then; z=random_beta(shape1,shape2)
    else; k=random_poisson(0.5_dp*nc); z=random_beta(shape1+real(k,dp),shape2); end if
    r=lo+(hi-lo)*z
  end function rbetagen

  subroutine betasubj_shapes(xmin,mode,mean,xmax,a1,a2,ok)
    real(dp),intent(in)::xmin,mode,mean,xmax
    real(dp),intent(out)::a1,a2
    logical,intent(out)::ok
    real(dp)::mid
    ok=.false.; a1=nan_dp(); a2=nan_dp()
    if(xmin>mode .or. mode>xmax .or. xmin>mean .or. mean>xmax .or. mean==mode .or. xmax<=xmin) return
    mid=0.5_dp*(xmin+xmax)
    if((mode>mean .and. mode<=mid) .or. (mode<mean .and. mode>=mid)) return
    a1=2*(mean-xmin)*(mid-mode)/((mean-mode)*(xmax-xmin)); a2=a1*(xmax-mean)/(mean-xmin)
    ok=a1>0 .and. a2>0
  end subroutine betasubj_shapes

  real(dp) function dbetasubj(x,xmin,mode,mean,xmax,log_density) result(d)
    real(dp),intent(in)::x,xmin,mode,mean,xmax
    logical,intent(in),optional::log_density
    real(dp)::a1,a2; logical::ok
    call betasubj_shapes(xmin,mode,mean,xmax,a1,a2,ok); if(.not.ok) then; d=nan_dp(); return; end if
    d=dbetagen(x,a1,a2,xmin,xmax,log_density=log_density)
  end function dbetasubj
  real(dp) function pbetasubj(q,xmin,mode,mean,xmax,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,xmin,mode,mean,xmax
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::a1,a2; logical::ok
    call betasubj_shapes(xmin,mode,mean,xmax,a1,a2,ok); if(.not.ok) then; p=nan_dp(); return; end if
    p=pbetagen(q,a1,a2,xmin,xmax,lower_tail=lower_tail,log_p=log_p)
  end function pbetasubj
  real(dp) function qbetasubj(p,xmin,mode,mean,xmax,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,xmin,mode,mean,xmax
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::a1,a2; logical::ok
    call betasubj_shapes(xmin,mode,mean,xmax,a1,a2,ok); if(.not.ok) then; q=nan_dp(); return; end if
    q=qbetagen(p,a1,a2,xmin,xmax,lower_tail=lower_tail,log_p=log_p)
  end function qbetasubj
  real(dp) function rbetasubj(xmin,mode,mean,xmax) result(r)
    real(dp),intent(in)::xmin,mode,mean,xmax
    real(dp)::a1,a2; logical::ok
    call betasubj_shapes(xmin,mode,mean,xmax,a1,a2,ok); if(.not.ok) then; r=nan_dp(); return; end if
    r=rbetagen(a1,a2,xmin,xmax)
  end function rbetasubj

  subroutine lognormal_params(mean,sd,meanlog,sdlog,ok)
    real(dp),intent(in)::mean,sd; real(dp),intent(out)::meanlog,sdlog; logical,intent(out)::ok
    ok=mean>0 .and. sd>=0
    if(.not.ok) then; meanlog=nan_dp(); sdlog=nan_dp(); return; end if
    meanlog=log(mean*mean/sqrt(sd*sd+mean*mean)); sdlog=sqrt(log(1+sd*sd/(mean*mean)))
  end subroutine lognormal_params
  real(dp) function dlnormb(x,mean,sd,log_density) result(d)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mean,sd
    logical,intent(in),optional::log_density
    real(dp)::m,s,ld,mu,sigma
    logical::ok,llog
    mu=exp(0.5_dp); if(present(mean))mu=mean
    sigma=sqrt(exp(2.0_dp)-exp(1.0_dp)); if(present(sd))sigma=sd
    call lognormal_params(mu,sigma,m,s,ok)
    llog=.false.; if(present(log_density))llog=log_density
    if(.not.ok)then;d=nan_dp();return;end if
    if(x<=0)then;d=0;if(llog)d=-huge(1.0_dp);return;end if
    if(s==0)then
      d=merge(huge(1.0_dp),0.0_dp,is_close(x,exp(m)))
      if(llog)d=log(d)
      return
    end if
    ld=-log(x*s*sqrt(2*pi))-0.5_dp*((log(x)-m)/s)**2
    d=merge(ld,exp(ld),llog)
  end function dlnormb

  real(dp) function plnormb(q,mean,sd,lower_tail,log_p) result(p)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::mean,sd
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,mu,sigma
    logical::ok,lt,lp
    mu=exp(0.5_dp); if(present(mean))mu=mean
    sigma=sqrt(exp(2.0_dp)-exp(1.0_dp)); if(present(sd))sigma=sd
    call lognormal_params(mu,sigma,m,s,ok)
    if(.not.ok)then;p=nan_dp();return;end if
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(q<=0)then
      p=0
    else if(s==0)then
      p=merge(1.0_dp,0.0_dp,q>=exp(m))
    else
      p=normal_cdf((log(q)-m)/s)
    end if
    if(.not.lt)p=1-p
    if(lp)then;if(p>0)p=log(p);if(p==0)p=-huge(1.0_dp);end if
  end function plnormb

  real(dp) function qlnormb(p,mean,sd,lower_tail,log_p) result(q)
    real(dp),intent(in)::p
    real(dp),intent(in),optional::mean,sd
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::m,s,pp,mu,sigma
    logical::ok,lt,lp
    mu=exp(0.5_dp); if(present(mean))mu=mean
    sigma=sqrt(exp(2.0_dp)-exp(1.0_dp)); if(present(sd))sigma=sd
    call lognormal_params(mu,sigma,m,s,ok)
    if(.not.ok)then;q=nan_dp();return;end if
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1-pp
    q=exp(m+s*normal_quantile(pp))
  end function qlnormb

  real(dp) function rlnormb(mean,sd) result(r)
    real(dp),intent(in),optional::mean,sd
    real(dp)::m,s,mu,sigma
    logical::ok
    mu=exp(0.5_dp); if(present(mean))mu=mean
    sigma=sqrt(exp(2.0_dp)-exp(1.0_dp)); if(present(sd))sigma=sd
    call lognormal_params(mu,sigma,m,s,ok)
    if(.not.ok)then;r=nan_dp();else;r=exp(m+s*random_normal());end if
  end function rlnormb

  subroutine pert_shapes(xmin,mode,xmax,shape,a1,a2,ok)
    real(dp),intent(in)::xmin,mode,xmax,shape;real(dp),intent(out)::a1,a2;logical,intent(out)::ok
    ok=xmax>xmin .and. mode>=xmin .and. mode<=xmax .and. shape>0
    if(ok)then;a1=1+shape*(mode-xmin)/(xmax-xmin);a2=1+shape*(xmax-mode)/(xmax-xmin);else;a1=nan_dp();a2=nan_dp();end if
  end subroutine pert_shapes
  real(dp) function dpert(x,xmin,mode,xmax,shape,log_density) result(d)
    real(dp),intent(in)::x,xmin,mode,xmax;real(dp),intent(in),optional::shape;logical,intent(in),optional::log_density
    real(dp)::sh,a1,a2;logical::ok
    sh=4;if(present(shape))sh=shape;call pert_shapes(xmin,mode,xmax,sh,a1,a2,ok);if(.not.ok)then;d=nan_dp();return;end if
    d=dbetagen(x,a1,a2,xmin,xmax,log_density=log_density)
  end function dpert
  real(dp) function ppert(q,xmin,mode,xmax,shape,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,xmin,mode,xmax;real(dp),intent(in),optional::shape;logical,intent(in),optional::lower_tail,log_p
    real(dp)::sh,a1,a2;logical::ok
    sh=4;if(present(shape))sh=shape;call pert_shapes(xmin,mode,xmax,sh,a1,a2,ok);if(.not.ok)then;p=nan_dp();return;end if
    p=pbetagen(q,a1,a2,xmin,xmax,lower_tail=lower_tail,log_p=log_p)
  end function ppert
  real(dp) function qpert(p,xmin,mode,xmax,shape,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,xmin,mode,xmax;real(dp),intent(in),optional::shape;logical,intent(in),optional::lower_tail,log_p
    real(dp)::sh,a1,a2;logical::ok
    sh=4;if(present(shape))sh=shape;call pert_shapes(xmin,mode,xmax,sh,a1,a2,ok);if(.not.ok)then;q=nan_dp();return;end if
    q=qbetagen(p,a1,a2,xmin,xmax,lower_tail=lower_tail,log_p=log_p)
  end function qpert
  real(dp) function rpert(xmin,mode,xmax,shape) result(r)
    real(dp),intent(in)::xmin,mode,xmax;real(dp),intent(in),optional::shape
    real(dp)::sh,a1,a2;logical::ok
    sh=4;if(present(shape))sh=shape;call pert_shapes(xmin,mode,xmax,sh,a1,a2,ok);if(.not.ok)then;r=nan_dp();return;end if
    r=rbetagen(a1,a2,xmin,xmax)
  end function rpert

  real(dp) function dpert_mean(x,xmin,mean,xmax,shape,log_density) result(d)
    real(dp),intent(in)::x,xmin,mean,xmax
    real(dp),intent(in),optional::shape
    logical,intent(in),optional::log_density
    real(dp)::sh,mode
    sh=4.0_dp; if(present(shape)) sh=shape
    mode=((sh+2.0_dp)*mean-xmin-xmax)/sh
    d=dpert(x,xmin,mode,xmax,sh,log_density)
  end function dpert_mean
  real(dp) function ppert_mean(x,xmin,mean,xmax,shape,lower_tail,log_p) result(p)
    real(dp),intent(in)::x,xmin,mean,xmax
    real(dp),intent(in),optional::shape
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::sh,mode
    sh=4.0_dp; if(present(shape)) sh=shape
    mode=((sh+2.0_dp)*mean-xmin-xmax)/sh
    p=ppert(x,xmin,mode,xmax,sh,lower_tail,log_p)
  end function ppert_mean
  real(dp) function qpert_mean(p,xmin,mean,xmax,shape,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,xmin,mean,xmax
    real(dp),intent(in),optional::shape
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::sh,mode
    sh=4.0_dp; if(present(shape)) sh=shape
    mode=((sh+2.0_dp)*mean-xmin-xmax)/sh
    q=qpert(p,xmin,mode,xmax,sh,lower_tail,log_p)
  end function qpert_mean
  real(dp) function rpert_mean(xmin,mean,xmax,shape) result(r)
    real(dp),intent(in)::xmin,mean,xmax
    real(dp),intent(in),optional::shape
    real(dp)::sh,mode
    sh=4.0_dp; if(present(shape)) sh=shape
    mode=((sh+2.0_dp)*mean-xmin-xmax)/sh
    r=rpert(xmin,mode,xmax,sh)
  end function rpert_mean

  real(dp) function dtriang(x,xmin,mode,xmax,log_density) result(d)
    real(dp),intent(in)::x,xmin,mode,xmax;logical,intent(in),optional::log_density
    logical::llog
    llog=.false.;if(present(log_density))llog=log_density
    if(xmax<=xmin .or. mode<xmin .or. mode>xmax)then;d=nan_dp();return;end if
    if(x<xmin .or. x>xmax)then;d=0
    else if(x<mode .or. mode>=xmax)then;d=2*(x-xmin)/((mode-xmin)*(xmax-xmin))
    else;d=2*(xmax-x)/((xmax-mode)*(xmax-xmin));end if
    if(llog)then;if(d>0)d=log(d);if(d==0)d=-huge(1.0_dp);end if
  end function dtriang
  real(dp) function ptriang(q,xmin,mode,xmax,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,xmin,mode,xmax;logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p
    if(xmax<xmin .or. mode<xmin .or. mode>xmax)then;p=nan_dp();return;end if
    if(q<=xmin)then;p=0;else if(q>=xmax)then;p=1
    else if(q<mode .or. mode>=xmax)then;p=(q-xmin)**2/((mode-xmin)*(xmax-xmin))
    else;p=1-(xmax-q)**2/((xmax-mode)*(xmax-xmin));end if
    if(.not.lt)p=1-p;if(lp)then;if(p>0)p=log(p);if(p==0)p=-huge(1.0_dp);end if
  end function ptriang
  real(dp) function qtriang(p,xmin,mode,xmax,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,xmin,mode,xmax;logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp,cut;logical::lt,lp
    lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1-pp
    if(pp<0 .or. pp>1 .or. xmax<xmin .or. mode<xmin .or. mode>xmax)then;q=nan_dp();return;end if
    if(xmax==xmin)then;q=xmin;return;end if
    cut=(mode-xmin)/(xmax-xmin)
    if(pp<=cut)then;q=xmin+sqrt(pp*(mode-xmin)*(xmax-xmin));else;q=xmax-sqrt((1-pp)*(xmax-xmin)*(xmax-mode));end if
  end function qtriang
  real(dp) function rtriang(xmin,mode,xmax) result(r)
    real(dp),intent(in)::xmin,mode,xmax;r=qtriang(random_uniform_open(),xmin,mode,xmax)
  end function rtriang

  real(dp) function dtriang_mean(x,xmin,mean,xmax,log_density) result(d)
    real(dp),intent(in)::x,xmin,mean,xmax
    logical,intent(in),optional::log_density
    d=dtriang(x,xmin,3.0_dp*mean-xmin-xmax,xmax,log_density)
  end function dtriang_mean
  real(dp) function ptriang_mean(x,xmin,mean,xmax,lower_tail,log_p) result(p)
    real(dp),intent(in)::x,xmin,mean,xmax
    logical,intent(in),optional::lower_tail,log_p
    p=ptriang(x,xmin,3.0_dp*mean-xmin-xmax,xmax,lower_tail,log_p)
  end function ptriang_mean
  real(dp) function qtriang_mean(p,xmin,mean,xmax,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,xmin,mean,xmax
    logical,intent(in),optional::lower_tail,log_p
    q=qtriang(p,xmin,3.0_dp*mean-xmin-xmax,xmax,lower_tail,log_p)
  end function qtriang_mean
  real(dp) function rtriang_mean(xmin,mean,xmax) result(r)
    real(dp),intent(in)::xmin,mean,xmax
    r=rtriang(xmin,3.0_dp*mean-xmin-xmax,xmax)
  end function rtriang_mean

  subroutine mqi_bounds(mqi,realization,k,intrinsic,vals,probs,ok)
    real(dp),intent(in)::mqi(3)
    real(dp),intent(in),optional::realization,k,intrinsic(2)
    real(dp),intent(out)::vals(5),probs(5);logical,intent(out)::ok
    real(dp)::kk,l,u
    kk=0.1_dp;if(present(k))kk=k;l=mqi(1);u=mqi(3)
    if(present(realization))then;l=min(l,realization);u=max(u,realization);end if
    vals(1)=l-kk*(u-l);vals(5)=u+kk*(u-l);if(present(intrinsic))vals([1,5])=intrinsic
    vals(2:4)=mqi;probs=[0.0_dp,0.05_dp,0.5_dp,0.95_dp,1.0_dp]
    ok=all(mqi(1:2)<=mqi(2:3)) .and. vals(1)<=mqi(1) .and. vals(5)>=mqi(3)
  end subroutine mqi_bounds
  real(dp) function dmqi(x,mqi,mqi_quantile,realization,k,intrinsic,log_density) result(d)
    real(dp),intent(in)::x,mqi(3);real(dp),intent(in),optional::mqi_quantile(3),realization,k,intrinsic(2)
    logical,intent(in),optional::log_density
    real(dp)::v(5),p(5);logical::ok,llog;integer::i
    call mqi_bounds(mqi,realization,k,intrinsic,v,p,ok);if(present(mqi_quantile))p(2:4)=mqi_quantile
    if(.not.ok)then;d=nan_dp();return;end if;d=0
    do i=1,4;if(x>=v(i) .and. x<v(i+1))then;d=(p(i+1)-p(i))/(v(i+1)-v(i));exit;end if;end do
    llog=.false.;if(present(log_density))llog=log_density;if(llog)then;if(d>0)d=log(d);if(d==0)d=-huge(1.0_dp);end if
  end function dmqi
  real(dp) function pmqi(q,mqi,mqi_quantile,realization,k,intrinsic,lower_tail,log_p) result(pout)
    real(dp),intent(in)::q,mqi(3);real(dp),intent(in),optional::mqi_quantile(3),realization,k,intrinsic(2)
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::v(5),p(5);logical::ok,lt,lp;integer::i
    call mqi_bounds(mqi,realization,k,intrinsic,v,p,ok);if(present(mqi_quantile))p(2:4)=mqi_quantile
    if(.not.ok)then;pout=nan_dp();return;end if
    if(q<=v(1))then;pout=0;else if(q>=v(5))then;pout=1;else
      do i=1,4;if(q<v(i+1))then;pout=p(i)+(p(i+1)-p(i))*(q-v(i))/(v(i+1)-v(i));exit;end if;end do
    end if
    lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p
    if(.not.lt)pout=1-pout;if(lp)then;if(pout>0)pout=log(pout);if(pout==0)pout=-huge(1.0_dp);end if
  end function pmqi
  real(dp) function qmqi(pin,mqi,mqi_quantile,realization,k,intrinsic,lower_tail,log_p) result(q)
    real(dp),intent(in)::pin,mqi(3);real(dp),intent(in),optional::mqi_quantile(3),realization,k,intrinsic(2)
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::v(5),p(5),pp;logical::ok,lt,lp;integer::i
    call mqi_bounds(mqi,realization,k,intrinsic,v,p,ok);if(present(mqi_quantile))p(2:4)=mqi_quantile
    q=nan_dp()
    lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p
    pp=pin;if(lp)pp=exp(pp);if(.not.lt)pp=1-pp
    if(.not.ok .or. pp<0 .or. pp>1)then;q=nan_dp();return;end if
    if(pp==1)then;q=v(5);return;end if
    do i=1,4;if(pp<p(i+1) .or. i==4)then;q=v(i)+(pp-p(i))*(v(i+1)-v(i))/(p(i+1)-p(i));return;end if;end do
  end function qmqi
  real(dp) function rmqi(mqi,mqi_quantile,realization,k,intrinsic) result(r)
    real(dp),intent(in)::mqi(3);real(dp),intent(in),optional::mqi_quantile(3),realization,k,intrinsic(2)
    r=qmqi(random_uniform_open(),mqi,mqi_quantile,realization,k,intrinsic)
  end function rmqi

  subroutine prepare_discrete(values,prob,uval,upr,nout,ok)
    real(dp),intent(in)::values(:),prob(:)
    real(dp),allocatable,intent(out)::uval(:),upr(:)
    integer,intent(out)::nout
    logical,intent(out)::ok
    integer,allocatable::idx(:);integer::i,n;real(dp)::s
    n=size(values);ok=n>0 .and. size(prob)==n .and. all(prob>=0) .and. sum(prob)>0
    if(.not.ok)then;nout=0;allocate(uval(0),upr(0));return;end if
    allocate(idx(n));call argsort_dp(values,idx);allocate(uval(n),upr(n));nout=0
    do i=1,n
      if(nout==0)then
        nout=nout+1
        uval(nout)=values(idx(i))
        upr(nout)=prob(idx(i))
      else if(values(idx(i))/=uval(nout))then
        nout=nout+1
        uval(nout)=values(idx(i))
        upr(nout)=prob(idx(i))
      else
        upr(nout)=upr(nout)+prob(idx(i))
      end if
    end do
    s=sum(upr(1:nout));upr(1:nout)=upr(1:nout)/s
  end subroutine prepare_discrete
  real(dp) function dempiricald(x,values,prob,log_density) result(d)
    real(dp),intent(in)::x,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::log_density
    real(dp),allocatable::pr(:),uv(:),up(:);integer::nu,i;logical::ok,llog
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_discrete(values,pr,uv,up,nu,ok);if(.not.ok)then;d=nan_dp();return;end if;d=0
    do i=1,nu;if(x==uv(i))then;d=up(i);exit;end if;end do
    llog=.false.;if(present(log_density))llog=log_density;if(llog)then;if(d>0)d=log(d);if(d==0)d=-huge(1.0_dp);end if
  end function dempiricald
  real(dp) function pempiricald(q,values,prob,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::lower_tail,log_p
    real(dp),allocatable::pr(:),uv(:),up(:);integer::nu,i;logical::ok,lt,lp
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_discrete(values,pr,uv,up,nu,ok);if(.not.ok)then;p=nan_dp();return;end if;p=0
    do i=1,nu;if(uv(i)<=q)p=p+up(i);end do
    lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p;if(.not.lt)p=1-p
    if(lp)then;if(p>0)p=log(p);if(p==0)p=-huge(1.0_dp);end if
  end function pempiricald
  real(dp) function qempiricald(pin,values,prob,lower_tail,log_p) result(q)
    real(dp),intent(in)::pin,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::lower_tail,log_p
    real(dp),allocatable::pr(:),uv(:),up(:);integer::nu,i;logical::ok,lt,lp;real(dp)::p,c
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_discrete(values,pr,uv,up,nu,ok)
    q=nan_dp()
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    lp=.false.
    if(present(log_p))lp=log_p
    p=pin;if(lp)p=exp(p);if(.not.lt)p=1-p;if(.not.ok .or. p<0 .or. p>1)then;q=nan_dp();return;end if
    c=0;do i=1,nu;c=c+up(i);if(p<=c .or. i==nu)then;q=uv(i);return;end if;end do
  end function qempiricald
  real(dp) function rempiricald(values,prob) result(r)
    real(dp),intent(in)::values(:);real(dp),intent(in),optional::prob(:);r=qempiricald(random_uniform_open(),values,prob)
  end function rempiricald

  subroutine prepare_continuous(xmin,xmax,values,prob,v,hgt,nu,integral,ok)
    real(dp),intent(in)::xmin,xmax,values(:),prob(:);real(dp),allocatable,intent(out)::v(:),hgt(:);integer,intent(out)::nu
    real(dp),intent(out)::integral;logical,intent(out)::ok
    real(dp),allocatable::uv(:),up(:);integer::i;logical::dok
    call prepare_discrete(values,prob,uv,up,nu,dok);ok=dok .and. xmin<=minval(values) .and. xmax>=maxval(values) .and. xmin<xmax
    if(.not.ok)then;allocate(v(0),hgt(0));integral=0;return;end if
    allocate(v(nu+2),hgt(nu+2));v(1)=xmin;v(2:nu+1)=uv(1:nu);v(nu+2)=xmax;hgt(1)=0;hgt(2:nu+1)=up(1:nu);hgt(nu+2)=0
    integral=0;do i=1,nu+1;integral=integral+(v(i+1)-v(i))*(hgt(i)+hgt(i+1))/2;end do
    hgt=hgt/integral
  end subroutine prepare_continuous
  real(dp) function dempiricalc(x,xmin,xmax,values,prob,log_density) result(d)
    real(dp),intent(in)::x,xmin,xmax,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::log_density
    real(dp),allocatable::pr(:),v(:),h(:);real(dp)::integ,t;integer::nu,i;logical::ok,llog
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_continuous(xmin,xmax,values,pr,v,h,nu,integ,ok);if(.not.ok)then;d=nan_dp();return;end if
    if(x<=xmin .or. x>=xmax)then
    d=0
    else
    d=0
    do i=1,nu+1
    if(x<v(i+1))then
    t=(x-v(i))/(v(i+1)-v(i))
    d=h(i)+t*(h(i+1)-h(i))
    exit
    end if
    end do
    end if
    llog=.false.;if(present(log_density))llog=log_density;if(llog)then;if(d>0)d=log(d);if(d==0)d=-huge(1.0_dp);end if
  end function dempiricalc
  real(dp) function pempiricalc(q,xmin,xmax,values,prob,lower_tail,log_p) result(p)
    real(dp),intent(in)::q,xmin,xmax,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::lower_tail,log_p
    real(dp),allocatable::pr(:),v(:),h(:);real(dp)::integ,dx,t;integer::nu,i;logical::ok,lt,lp
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_continuous(xmin,xmax,values,pr,v,h,nu,integ,ok);if(.not.ok)then;p=nan_dp();return;end if
    if(q<=xmin)then;p=0;else if(q>=xmax)then;p=1;else
      p=0;do i=1,nu+1
        if(q>=v(i+1))then;p=p+(v(i+1)-v(i))*(h(i)+h(i+1))/2
        else;dx=q-v(i);t=dx/(v(i+1)-v(i));p=p+dx*(h(i)+h(i)+t*(h(i+1)-h(i)))/2;exit;end if
      end do
    end if
    p=clamp01(p);lt=.true.;if(present(lower_tail))lt=lower_tail;lp=.false.;if(present(log_p))lp=log_p;if(.not.lt)p=1-p
    if(lp)then;if(p>0)p=log(p);if(p==0)p=-huge(1.0_dp);end if
  end function pempiricalc
  real(dp) function qempiricalc(pin,xmin,xmax,values,prob,lower_tail,log_p) result(q)
    real(dp),intent(in)::pin,xmin,xmax,values(:);real(dp),intent(in),optional::prob(:);logical,intent(in),optional::lower_tail,log_p
    real(dp),allocatable::pr(:),v(:),h(:);real(dp)::integ,p,c,nextc,dx,a,b,disc;integer::nu,i;logical::ok,lt,lp
    if(present(prob))then;pr=prob;else;allocate(pr(size(values)));pr=1;end if
    call prepare_continuous(xmin,xmax,values,pr,v,h,nu,integ,ok)
    lt=.true.
    if(present(lower_tail))lt=lower_tail
    lp=.false.
    if(present(log_p))lp=log_p
    p=pin;if(lp)p=exp(p);if(.not.lt)p=1-p;if(.not.ok .or. p<0 .or. p>1)then;q=nan_dp();return;end if
    if(p==0)then;q=xmin;return;else if(p==1)then;q=xmax;return;end if
    c=0;do i=1,nu+1
      nextc=c+(v(i+1)-v(i))*(h(i)+h(i+1))/2
      if(p<=nextc)then
        dx=v(i+1)-v(i);a=(h(i+1)-h(i))/(2*dx);b=h(i)
        if(abs(a)<100*epsilon(1.0_dp))then;q=v(i)+(p-c)/b
        else;disc=max(0.0_dp,b*b+4*a*(p-c));q=v(i)+(-b+sqrt(disc))/(2*a);end if
        return
      end if;c=nextc
    end do;q=xmax
  end function qempiricalc
  real(dp) function rempiricalc(xmin,xmax,values,prob) result(r)
    real(dp),intent(in)::xmin,xmax,values(:)
    real(dp),intent(in),optional::prob(:)
    r=qempiricalc(random_uniform_open(),xmin,xmax,values,prob)
  end function rempiricalc

  real(dp) function ddirichlet(x,alpha,log_density) result(d)
    real(dp),intent(in)::x(:),alpha(:);logical,intent(in),optional::log_density
    real(dp)::ld;logical::llog
    llog=.false.;if(present(log_density))llog=log_density
    if(size(x)/=size(alpha) .or. any(alpha<=0) .or. any(x<0) .or. any(x>1) .or. abs(sum(x)-1)>1e-10_dp)then
    d=0
    if(llog)d=-huge(1.0_dp)
    return
    end if
    ld=log_gamma(sum(alpha))-sum(log_gamma(alpha))+sum((alpha-1)*log(max(x,tiny(1.0_dp))))
    if(llog)then;d=ld;else;d=exp(ld);end if
  end function ddirichlet
  subroutine rdirichlet(alpha,x)
    real(dp),intent(in)::alpha(:);real(dp),intent(out)::x(size(alpha));integer::i
    do i=1,size(alpha);x(i)=random_gamma(alpha(i));end do;x=x/sum(x)
  end subroutine rdirichlet

  real(dp) function dmultinomial(x,prob,log_density) result(d)
    integer,intent(in)::x(:);real(dp),intent(in)::prob(:);logical,intent(in),optional::log_density
    real(dp),allocatable::p(:);real(dp)::ld;integer::n,i;logical::llog
    llog=.false.;if(present(log_density))llog=log_density
    if(size(x)/=size(prob) .or. any(x<0) .or. any(prob<0) .or. sum(prob)<=0)then;d=nan_dp();return;end if
    p=prob/sum(prob);n=sum(x);ld=log_gamma(real(n+1,dp))
    do i=1,size(x)
      if(p(i)==0 .and. x(i)>0)then;d=0;if(llog)d=-huge(1.0_dp);return;end if
      if(x(i)>0)ld=ld+real(x(i),dp)*log(p(i));ld=ld-log_gamma(real(x(i)+1,dp))
    end do
    if(llog)then;d=ld;else;d=exp(ld);end if
  end function dmultinomial
  subroutine rmultinomial(nsize,prob,x)
    integer,intent(in)::nsize;real(dp),intent(in)::prob(:);integer,intent(out)::x(size(prob))
    real(dp)::remaining_p,pcond;integer::remaining_n,i,k
    if(nsize<0 .or. any(prob<0) .or. sum(prob)<=0)then;x=-1;return;end if
    remaining_n=nsize;remaining_p=sum(prob);x=0
    do i=1,size(prob)-1
      if(remaining_n<=0)exit
      pcond=prob(i)/remaining_p
      k=random_binomial(remaining_n,pcond)
      x(i)=k
      remaining_n=remaining_n-k
      remaining_p=remaining_p-prob(i)
    end do;x(size(prob))=remaining_n
  end subroutine rmultinomial

  real(dp) function dmultinormal(x,mean,sigma,log_density) result(d)
    real(dp),intent(in)::x(:),mean(:),sigma(:,:)
    logical,intent(in),optional::log_density
    logical::llog
    llog=.false.;if(present(log_density))llog=log_density
    d=dmvnorm_one(x,mean,sigma,llog)
  end function dmultinormal
  subroutine rmultinormal(n,mean,sigma,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer,intent(in),optional::seed
    x=rmvnorm(n,mean,sigma,seed)
  end subroutine rmultinormal
  subroutine dmultinormal_varying(x,mean,sigma,d,log_density)
    real(dp),intent(in)::x(:,:),mean(:,:),sigma(:,:,:)
    real(dp),intent(out)::d(size(x,1))
    logical,intent(in),optional::log_density
    integer::i,nm,ns
    logical::llog
    nm=size(mean,1);ns=size(sigma,3);llog=.false.
    if(present(log_density))llog=log_density
    do i=1,size(x,1)
      d(i)=dmvnorm_one(x(i,:),mean(1+mod(i-1,nm),:), &
        sigma(:,:,1+mod(i-1,ns)),llog)
    end do
  end subroutine dmultinormal_varying
  subroutine rmultinormal_varying(n,mean,sigma,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:,:),sigma(:,:,:)
    real(dp),intent(out)::x(n,size(mean,2))
    integer,intent(in),optional::seed
    real(dp),allocatable::one(:,:);integer::i,nm,ns
    if(present(seed))call seed_rng_local(seed);nm=size(mean,1);ns=size(sigma,3)
    do i=1,n;one=rmvnorm(1,mean(1+mod(i-1,nm),:),sigma(:,:,1+mod(i-1,ns)));x(i,:)=one(1,:);end do
  contains
    subroutine seed_rng_local(s)
      use mc2d_random,only:seed_random
      integer,intent(in)::s;call seed_random(s)
    end subroutine
  end subroutine rmultinormal_varying

  subroutine rtrunc(pfun,qfun,n,linf,lsup,x)
    procedure(scalar_cdf)::pfun;procedure(scalar_quantile)::qfun
    integer,intent(in)::n;real(dp),intent(in)::linf,lsup;real(dp),intent(out)::x(n)
    real(dp)::plo,phi,u;integer::i
    if(linf>=lsup)then;x=nan_dp();return;end if
    plo=pfun(linf);phi=pfun(lsup)
    do i=1,n;call random_number(u);u=plo+(phi-plo)*u;x(i)=qfun(u);if(x(i)<=linf .or. x(i)>lsup)x(i)=nan_dp();end do
  end subroutine rtrunc

  subroutine lhs(qfun,nsv,nsu,nvariates,x)
    procedure(scalar_quantile)::qfun
    integer,intent(in)::nsv,nsu,nvariates;real(dp),intent(out)::x(nsv,nsu,nvariates)
    integer,allocatable::perm(:);integer::i,j,k,m,t;real(dp)::u
    allocate(perm(nsv))
    do k=1,nvariates;do j=1,nsu
      do i=1,nsv;perm(i)=i;end do
      do i=nsv,2,-1;call random_number(u);m=1+int(u*real(i,dp));if(m>i)m=i;t=perm(i);perm(i)=perm(m);perm(m)=t;end do
      do i=1,nsv;call random_number(u);x(i,j,k)=qfun((real(perm(i)-1,dp)+u)/real(nsv,dp));end do
    end do;end do
  end subroutine lhs
end module mc2d_distributions
