! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_distributions
  use degreenet_kinds, only : dp, pi
  use degreenet_math, only : zeta_r, logdiffexp, gauss_hermite, &
    nbinom_logpmf, nbinom_pmf
  use degreenet_rng, only : runif01, rpoisson_basic
  implicit none
  private
  public :: dyule, ldyule, ddp, lddp, dtp, ldtp
  public :: dwar, ldwar, ddqe, lddqe, dghdi, ldghdi
  public :: dnb_degree, ldnb_degree, ddpe_weight, dpe_pmf
  public :: cmp_z, dcmp_natural, cmp_moments, cmp_mutonatural
  public :: dcmp_mu, rcmp_mu
  public :: dpln, ldpln, simpln_one
  public :: rdiscrete_pmf

contains
  real(dp) function yule_base(rho,x) result(p)
    real(dp),intent(in)::rho
    integer,intent(in)::x
    if(rho<=1.0_dp.or.x<1) then
      p=0.0_dp
    else
      p=exp(log(rho-1.0_dp)+log_gamma(real(x,dp))+log_gamma(rho)- &
        log_gamma(real(x,dp)+rho))
    end if
  end function yule_base

  real(dp) function dyule(rho,x,cutoff) result(p)
    real(dp),intent(in)::rho
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j
    real(dp)::z
    c=1; if(present(cutoff)) c=cutoff
    if(x<c) then; p=0.0_dp; return; end if
    p=yule_base(rho,x)
    if(c>1) then
      z=1.0_dp
      do j=1,c-1; z=z-yule_base(rho,j); end do
      if(z>0.0_dp) p=p/z
    end if
  end function dyule

  real(dp) function ldyule(rho,x,cutoff) result(lp)
    real(dp),intent(in)::rho
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::p
    p=dyule(rho,x,cutoff)
    lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldyule

  real(dp) function dp_base(alpha,x) result(p)
    real(dp),intent(in)::alpha
    integer,intent(in)::x
    if(alpha<=1.0_dp.or.x<1) then; p=0.0_dp
    else; p=real(x,dp)**(-alpha)/zeta_r(alpha); end if
  end function dp_base

  real(dp) function ddp(alpha,x,cutoff) result(p)
    real(dp),intent(in)::alpha
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j
    real(dp)::z
    c=1; if(present(cutoff)) c=cutoff
    if(x<c) then; p=0.0_dp; return; end if
    p=dp_base(alpha,x)
    if(c>1) then
      z=1.0_dp
      do j=1,c-1; z=z-dp_base(alpha,j); end do
      if(z>0.0_dp) p=p/z
    end if
  end function ddp

  real(dp) function lddp(alpha,x,cutoff) result(lp)
    real(dp),intent(in)::alpha
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::p
    p=ddp(alpha,x,cutoff); lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function lddp

  real(dp) function dtp(alpha,x,cutoff) result(p)
    real(dp),intent(in)::alpha
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j
    real(dp)::c0,z,raw
    c=1; if(present(cutoff)) c=cutoff
    if(x<1.or.alpha<=1.0_dp) then; p=0.0_dp; return; end if
    if(c<=1) then
      p=real(x,dp)**(-alpha)/zeta_r(alpha); return
    end if
    c0=0.0_dp
    do j=1,c-1
      c0=c0+exp(-alpha*(real(j,dp)/real(c,dp)-1.0_dp))
    end do
    c0=c0 + real(c,dp)**alpha/zeta_r(alpha)
    if(x<c) then
      raw=exp(-alpha*(real(x,dp)/real(c,dp)-1.0_dp))/c0
    else
      raw=(real(x,dp)/real(c,dp))**(-alpha)/c0
    end if
    if(x<c) then; p=0.0_dp; return; end if
    z=0.0_dp
    do j=c,c+200000
      if(j<c) cycle
      if(j<c) then
        z=z+exp(-alpha*(real(j,dp)/real(c,dp)-1.0_dp))/c0
      else
        z=z+(real(j,dp)/real(c,dp))**(-alpha)/c0
      end if
      if(j>c+1000.and.(real(j,dp)/real(c,dp))**(-alpha)/c0<1.0e-14_dp*z) exit
    end do
    p=raw/max(z,tiny(1.0_dp))
  end function dtp

  real(dp) function ldtp(alpha,x,cutoff) result(lp)
    real(dp),intent(in)::alpha
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::p
    p=dtp(alpha,x,cutoff); lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldtp

  real(dp) function war_base(rho,a,x) result(p)
    real(dp),intent(in)::rho,a
    integer,intent(in)::x
    if(rho<=1.0_dp.or.a<=-1.0_dp.or.x<1) then; p=0.0_dp; return; end if
    p=exp(log(rho-1.0_dp)+log_gamma(real(x,dp)+a)+log_gamma(rho+a)- &
      log_gamma(1.0_dp+a)-log_gamma(real(x,dp)+rho+a))
  end function war_base

  real(dp) function dwar(v,x,cutoff) result(p)
    real(dp),intent(in)::v(2)
    integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j
    real(dp)::z
    c=1;if(present(cutoff))c=cutoff
    if(x<c) then;p=0.0_dp;return;end if
    p=war_base(v(1),v(2),x)
    if(c>1) then
      z=1.0_dp; do j=1,c-1; z=z-war_base(v(1),v(2),j); end do
      if(z>0.0_dp)p=p/z
    end if
  end function dwar

  real(dp) function ldwar(v,x,cutoff) result(lp)
    real(dp),intent(in)::v(2); integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::p
    p=dwar(v,x,cutoff); lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldwar

  real(dp) function dqe_base(alpha,beta,x) result(p)
    real(dp),intent(in)::alpha,beta; integer,intent(in)::x
    real(dp)::a,b
    if(alpha<=0.0_dp.or.beta<=0.0_dp.or.x<1) then;p=0.0_dp;return;end if
    a=-alpha*log(1.0_dp+real(x-1,dp)/beta)
    b=-alpha*log(1.0_dp+real(x,dp)/beta)
    p=exp(logdiffexp(a,b))
  end function dqe_base

  real(dp) function ddqe(v,x,cutoff) result(p)
    real(dp),intent(in)::v(2); integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j; real(dp)::z
    c=1;if(present(cutoff))c=cutoff
    if(x<c) then;p=0.0_dp;return;end if
    p=dqe_base(v(1),v(2),x)
    if(c>1) then
      z=1.0_dp;do j=1,c-1;z=z-dqe_base(v(1),v(2),j);end do
      if(z>0.0_dp)p=p/z
    end if
  end function ddqe

  real(dp) function lddqe(v,x,cutoff) result(lp)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::p
    p=ddqe(v,x,cutoff);lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function lddqe

  real(dp) function ghdi_lograw(v,x) result(lq)
    real(dp),intent(in)::v(4);integer,intent(in)::x
    real(dp)::m10,m20
    if(x<0.or.v(4)<=0.0_dp.or.v(2)<=0.0_dp.or.v(3)<=v(2).or.v(1)<=0.0_dp) then
      lq=-huge(1.0_dp);return
    end if
    m10=v(3)-v(2);m20=v(2)
    lq=real(x,dp)*log(v(4))+log_gamma(real(x,dp)+v(1))+log_gamma(m10)+ &
      log_gamma(m20+real(x,dp))-log_gamma(real(x+1,dp))- &
      log_gamma(m10+m20+real(x,dp))
  end function ghdi_lograw

  real(dp) function dghdi(v,x,cutoff,maxk) result(p)
    real(dp),intent(in)::v(4);integer,intent(in)::x
    integer,intent(in),optional::cutoff,maxk
    integer::c,k,mk;real(dp)::z,mx,lr
    c=1;if(present(cutoff))c=cutoff;mk=10000;if(present(maxk))mk=maxk
    if(x<c) then;p=0.0_dp;return;end if
    mx=-huge(1.0_dp)
    do k=c,mk;mx=max(mx,ghdi_lograw(v,k));end do
    if(.not.(mx>-huge(1.0_dp)))then;p=0.0_dp;return;end if
    z=0.0_dp
    do k=c,mk
      lr=ghdi_lograw(v,k)
      z=z+exp(lr-mx)
    end do
    p=exp(ghdi_lograw(v,x)-mx)/max(z,tiny(1.0_dp))
  end function dghdi

  real(dp) function ldghdi(v,x,cutoff,maxk) result(lp)
    real(dp),intent(in)::v(4);integer,intent(in)::x
    integer,intent(in),optional::cutoff,maxk;real(dp)::p
    p=dghdi(v,x,cutoff,maxk);lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldghdi

  real(dp) function dnb_degree(v,x,cutoff,shifted) result(pmf)
    real(dp),intent(in)::v(2); integer,intent(in)::x
    integer,intent(in),optional::cutoff
    logical,intent(in),optional::shifted
    integer::c,k;logical::sh
    real(dp)::size,z
    c=0;if(present(cutoff))c=cutoff; sh=.false.;if(present(shifted))sh=shifted
    if(v(1)<=0.0_dp.or.v(2)<=0.0_dp.or.v(2)>1.0_dp) then;pmf=0.0_dp;return;end if
    size=v(1)*v(2); k=x; if(sh)k=x-c
    pmf=nbinom_pmf(k,size,v(2))
    if(.not.sh.and.c>0) then
      z=0.0_dp;do k=c,c+100000;z=z+nbinom_pmf(k,size,v(2));if(k>c+100.and.nbinom_pmf(k,size,v(2))<1e-14_dp*z)exit;end do
      pmf=pmf/max(z,tiny(1.0_dp))
    end if
  end function dnb_degree

  real(dp) function ldnb_degree(v,x,cutoff,shifted) result(lp)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff;logical,intent(in),optional::shifted
    real(dp)::p
    p=dnb_degree(v,x,cutoff,shifted);lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldnb_degree

  pure real(dp) function ddpe_weight(v,x) result(w)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    if(x<1.or.v(2)<=0.0_dp)then;w=0.0_dp;else;w=real(x,dp)**(-v(1))*exp(-real(x,dp)/v(2));end if
  end function ddpe_weight

  real(dp) function dpe_pmf(v,x,cutoff,maxk) result(p)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff,maxk
    integer::c,mk,k;real(dp)::z,t
    c=1;if(present(cutoff))c=cutoff;mk=100000;if(present(maxk))mk=maxk
    if(x<c)then;p=0.0_dp;return;end if
    z=0.0_dp
    do k=c,mk
      t=ddpe_weight(v,k);z=z+t
      if(k>c+100.and.t<1e-14_dp*z)exit
    end do
    p=ddpe_weight(v,x)/max(z,tiny(1.0_dp))
  end function dpe_pmf

  real(dp) function cmp_z(lambda,nu,err,kmin,log_value) result(z)
    real(dp),intent(in)::lambda,nu
    real(dp),intent(in),optional::err
    integer,intent(in),optional::kmin
    logical,intent(in),optional::log_value
    real(dp)::e,ratio,lsum,lterm
    integer::j,km
    logical::lv
    e=1e-10_dp;if(present(err))e=err;km=200;if(present(kmin))km=max(2,kmin)
    lv=.false.;if(present(log_value))lv=log_value
    if(lambda<0.0_dp.or.nu<0.0_dp)then;z=huge(1.0_dp);return;end if
    if(lambda<=0.0_dp)then;z=merge(0.0_dp,1.0_dp,lv);return;end if
    ! Log-sum-exp recurrence is stable in the high-mode regime.
    lterm=0.0_dp;lsum=0.0_dp
    do j=1,200000
      lterm=lterm+log(lambda)-nu*log(real(j,dp))
      if(lsum>lterm)then
        lsum=lsum+log(1.0_dp+exp(lterm-lsum))
      else
        lsum=lterm+log(1.0_dp+exp(lsum-lterm))
      end if
      if(j>=km)then
        ratio=lambda/real(j+1,dp)**nu
        if(ratio<1.0_dp.and.exp(lterm-lsum)<e*max(1.0_dp,1.0_dp-ratio))exit
      end if
    end do
    if(lv)then;z=lsum;else;z=exp(lsum);end if
  end function cmp_z

  recursive real(dp) function dcmp_natural(lambda,nu,x,cutoff,err) result(p)
    real(dp),intent(in)::lambda,nu;integer,intent(in)::x
    integer,intent(in),optional::cutoff;real(dp),intent(in),optional::err
    integer::c,j;real(dp)::lz,zcut
    c=0;if(present(cutoff))c=cutoff
    if(x<c.or.x<0.or.lambda<0.0_dp.or.nu<0.0_dp)then;p=0.0_dp;return;end if
    lz=cmp_z(lambda,nu,err,max(200,2*x),.true.)
    p=exp(real(x,dp)*log(max(lambda,tiny(1.0_dp)))-nu*log_gamma(real(x+1,dp))-lz)
    if(lambda<=0.0_dp)p=merge(1.0_dp,0.0_dp,x==0)
    if(c>0)then
      zcut=1.0_dp
      do j=0,c-1;zcut=zcut-dcmp_natural(lambda,nu,j,0,err);end do
      p=p/max(zcut,tiny(1.0_dp))
    end if
  end function dcmp_natural

  subroutine cmp_moments(lambda,nu,mu,sd,kmax)
    real(dp),intent(in)::lambda,nu;real(dp),intent(out)::mu,sd
    integer,intent(in),optional::kmax
    integer::k,km;real(dp)::p,m2,s
    km=2000;if(present(kmax))km=kmax
    mu=0.0_dp;m2=0.0_dp;s=0.0_dp
    do k=0,km
      p=dcmp_natural(lambda,nu,k,0,1e-12_dp);s=s+p
      mu=mu+real(k,dp)*p;m2=m2+real(k*k,dp)*p
      if(k>100.and.p<1e-14_dp*max(s,1e-30_dp))exit
    end do
    sd=sqrt(max(0.0_dp,m2-mu*mu))
  end subroutine cmp_moments

  subroutine cmp_mutonatural(mu_target,sd_target,lambda,nu,ok)
    real(dp),intent(in)::mu_target,sd_target
    real(dp),intent(out)::lambda,nu
    logical,intent(out),optional::ok
    real(dp)::a,b,mu,sd,f1,f2,ma,sa,mb,sb,h,det,da,db
    integer::it
    logical::good
    a=log(max(mu_target,1e-3_dp));b=0.0_dp;good=.false.
    do it=1,60
      call cmp_moments(exp(a),exp(b),mu,sd,1500)
      f1=mu-mu_target;f2=sd-sd_target
      if(sqrt(f1*f1+f2*f2)<1e-7_dp*max(1.0_dp,mu_target+sd_target))then;good=.true.;exit;end if
      h=1e-4_dp
      call cmp_moments(exp(a+h),exp(b),ma,sa,1200)
      call cmp_moments(exp(a),exp(b+h),mb,sb,1200)
      det=((ma-mu)/h)*((sb-sd)/h)-((mb-mu)/h)*((sa-sd)/h)
      if(abs(det)<1e-12_dp)exit
      da=(-f1*((sb-sd)/h)+f2*((mb-mu)/h))/det
      db=(-(ma-mu)/h*f2+(sa-sd)/h*f1)/det
      da=max(-1.0_dp,min(1.0_dp,da));db=max(-1.0_dp,min(1.0_dp,db))
      a=a+da;b=b+db
    end do
    lambda=exp(a);nu=exp(b);if(present(ok))ok=good
  end subroutine cmp_mutonatural

  real(dp) function dcmp_mu(v,x,cutoff) result(p)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff
    real(dp)::lambda,nu;logical::ok
    call cmp_mutonatural(v(1),v(2),lambda,nu,ok)
    p=dcmp_natural(lambda,nu,x,cutoff,1e-10_dp)
  end function dcmp_mu

  integer function rcmp_mu(v,kmax) result(x)
    real(dp),intent(in)::v(2);integer,intent(in),optional::kmax
    integer::km;real(dp)::lambda,nu,u,c;logical::ok
    km=1000;if(present(kmax))km=kmax
    call cmp_mutonatural(v(1),v(2),lambda,nu,ok)
    u=runif01();c=0.0_dp
    do x=0,km
      c=c+dcmp_natural(lambda,nu,x,0,1e-10_dp)
      if(u<=c)return
    end do
    x=km
  end function rcmp_mu

  real(dp) function dpln(v,x,cutoff,logn,points) result(p)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff,points;logical,intent(in),optional::logn
    integer::c,nq,j,k;logical::ln,ok
    real(dp)::mu0,sig0,meanv,sdv,lambda,raw,norm
    real(dp),allocatable::gx(:),gw(:)
    c=0;if(present(cutoff))c=cutoff;nq=32;if(present(points))nq=points
    ln=.true.;if(present(logn))ln=logn
    if(x<c.or.x<0)then;p=0.0_dp;return;end if
    if(ln)then;mu0=v(1);sig0=v(2)
    else
      meanv=v(1);sdv=v(2)
      if(meanv<=0.0_dp.or.sdv<=0.0_dp)then;p=0.0_dp;return;end if
      if(sdv*sdv>meanv)then
        sig0=sqrt(log(1.0_dp+(sdv*sdv-meanv)/(meanv*meanv)))
      else
        sig0=1e-7_dp
      end if
      mu0=log(meanv)-0.5_dp*sig0*sig0
    end if
    if(sig0<=0.0_dp)then;p=0.0_dp;return;end if
    allocate(gx(nq),gw(nq));call gauss_hermite(nq,gx,gw,ok)
    raw=0.0_dp
    do j=1,nq
      lambda=exp(mu0+sig0*gx(j))
      raw=raw+gw(j)*exp(real(x,dp)*log(lambda)-lambda-log_gamma(real(x+1,dp)))
    end do
    p=raw
    if(c>0)then
      norm=1.0_dp
      do k=0,c-1
        raw=0.0_dp
        do j=1,nq
          lambda=exp(mu0+sig0*gx(j))
          raw=raw+gw(j)*exp(real(k,dp)*log(lambda)-lambda-log_gamma(real(k+1,dp)))
        end do
        norm=norm-raw
      end do
      p=p/max(norm,tiny(1.0_dp))
    end if
    deallocate(gx,gw)
  end function dpln

  real(dp) function ldpln(v,x,cutoff,logn,points) result(lp)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff,points;logical,intent(in),optional::logn
    real(dp)::p
    p=dpln(v,x,cutoff,logn,points);lp=merge(log(p),-huge(1.0_dp),p>0.0_dp)
  end function ldpln

  integer function simpln_one(v) result(x)
    real(dp),intent(in)::v(2)
    real(dp)::lambda
    lambda=exp(v(1)+v(2)*sqrt(2.0_dp)*0.0_dp)
    ! Draw lognormal mixing variable explicitly.
    block
      use degreenet_rng, only : rnorm01
      lambda=exp(v(1)+v(2)*rnorm01())
    end block
    x=rpoisson_basic(lambda)
  end function simpln_one

  integer function rdiscrete_pmf(p) result(ix)
    real(dp),intent(in)::p(:)
    real(dp)::u,c,s
    integer::i
    s=sum(max(p,0.0_dp));if(s<=0.0_dp)then;ix=1;return;end if
    u=runif01()*s;c=0.0_dp
    do i=1,size(p);c=c+max(p(i),0.0_dp);if(u<=c)then;ix=i;return;end if;end do
    ix=size(p)
  end function rdiscrete_pmf
end module degreenet_distributions
