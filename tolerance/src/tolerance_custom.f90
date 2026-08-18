! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_custom
  use tolerance_kinds, only : dp
  use tolerance_types, only : discrete_tolerance_interval
  use tolerance_math, only : normal_quantile, invert_matrix
  use tolerance_distributions, only : ddpareto, qdpareto, dpoislind, qpoislind, &
       dzipfman, qzipfman
  use tolerance_optimize, only : golden_minimize, nelder_mead, numerical_hessian
  implicit none
  private
  public :: dpareto_mle, dparetotol_int, poislind_mle, poislindtol_int
  public :: zipf_mle, zipfman_mle, zeta_mle, zipftol_int
contains

  real(dp) function dpareto_mle(x,se) result(theta)
    integer,intent(in)::x(:)
    real(dp),intent(out),optional::se
    real(dp)::h,eps,f0,fp,fm
    theta=golden_minimize(nll,1.0e-12_dp,1.0_dp-1.0e-12_dp,1.0e-11_dp)
    if(present(se))then
      h=max(1.0e-5_dp,1.0e-4_dp*theta)
      h=min(h,0.49_dp*min(theta,1.0_dp-theta))
      f0=nll(theta);fp=nll(theta+h);fm=nll(theta-h)
      eps=(fp-2.0_dp*f0+fm)/(h*h)
      se=sqrt(1.0_dp/max(eps,1.0e-30_dp))
    end if
  contains
    real(dp) function nll(t) result(v)
      real(dp),intent(in)::t
      integer::i
      v=0.0_dp
      do i=1,size(x);v=v-ddpareto(x(i),t,.true.);end do
    end function nll
  end function dpareto_mle

  function dparetotol_int(x,m,alpha,p,side,theta_hat) result(out)
    integer,intent(in)::x(:)
    integer,intent(in),optional::m
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    real(dp),intent(out),optional::theta_hat
    type(discrete_tolerance_interval)::out
    real(dp)::a,pp,theta,se,lo_t,hi_t,z,ql,qu
    integer::n,mm,sd
    n=size(x);mm=n;if(present(m))mm=m;a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    theta=dpareto_mle(x,se);z=normal_quantile(1.0_dp-a)*se*sqrt(real(n,dp)/real(mm,dp))
    lo_t=max(theta-z,1.0e-14_dp);hi_t=min(theta+z,1.0_dp)
    ql=qdpareto(1.0_dp-pp,lo_t);qu=qdpareto(pp,min(hi_t,1.0_dp-1.0e-14_dp))
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2);out%estimate=theta
    out%lower=max(0,int(ql));out%upper_infinite=hi_t>=1.0_dp-1.0e-14_dp
    if(out%upper_infinite)then;out%upper=huge(out%upper);else;out%upper=max(0,int(qu));end if
    if(present(theta_hat))theta_hat=theta
  end function dparetotol_int

  real(dp) function poislind_mle(x,se) result(theta)
    integer,intent(in)::x(:)
    real(dp),intent(out),optional::se
    real(dp)::xb,start,hi,h,f0,fp,fm,cur
    integer::i
    xb=sum(real(x,dp))/real(size(x),dp)
    if(xb<=1.0e-15_dp)then;start=10.0_dp;else
      start=(-(xb-1.0_dp)+sqrt((xb-1.0_dp)**2+8.0_dp*xb))/(2.0_dp*xb)
    end if
    hi=max(10.0_dp,10.0_dp*start)
    cur=nll(hi)
    do i=1,20
      if(nll(2.0_dp*hi)>=cur)exit
      hi=2.0_dp*hi;cur=nll(hi)
    end do
    theta=golden_minimize(nll,1.0e-12_dp,hi,1.0e-10_dp)
    if(present(se))then
      h=max(1.0e-5_dp,1.0e-4_dp*theta);h=min(h,0.49_dp*theta)
      f0=nll(theta);fp=nll(theta+h);fm=nll(theta-h)
      se=sqrt(1.0_dp/max((fp-2.0_dp*f0+fm)/(h*h),1.0e-30_dp))
    end if
  contains
    real(dp) function nll(t) result(v)
      real(dp),intent(in)::t
      integer::j
      if(t<=0.0_dp)then;v=huge(1.0_dp);return;end if
      v=0.0_dp;do j=1,size(x);v=v-dpoislind(x(j),t,.true.);end do
    end function nll
  end function poislind_mle

  function poislindtol_int(x,m,alpha,p,side,theta_hat) result(out)
    integer,intent(in)::x(:)
    integer,intent(in),optional::m
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    real(dp),intent(out),optional::theta_hat
    type(discrete_tolerance_interval)::out
    real(dp)::a,pp,theta,se,z,lo_t,hi_t,ql,qu
    integer::n,mm,sd
    n=size(x);mm=n;if(present(m))mm=m;a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    theta=poislind_mle(x,se);z=normal_quantile(1.0_dp-a)*se*sqrt(real(n,dp)/real(mm,dp))
    lo_t=max(theta-z,1.0e-12_dp);hi_t=max(theta+z,1.0e-12_dp)
    ql=qpoislind(1.0_dp-pp,hi_t);qu=qpoislind(pp,lo_t)
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2);out%estimate=theta
    out%lower=max(0,int(ql));out%upper=max(0,int(qu));out%upper_infinite=.false.
    if(present(theta_hat))theta_hat=theta
  end function poislindtol_int

  real(dp) function zipf_mle(x,nmax,se) result(s)
    integer,intent(in)::x(:),nmax
    real(dp),intent(out),optional::se
    real(dp)::hi,h,f0,fp,fm
    hi=2.0_dp
    do while(nll(2.0_dp*hi)<nll(hi) .and. hi<100.0_dp);hi=2.0_dp*hi;end do
    s=golden_minimize(nll,1.0e-10_dp,hi,1.0e-10_dp)
    if(present(se))then
      h=max(1.0e-5_dp,1.0e-4_dp*s);h=min(h,0.49_dp*s)
      f0=nll(s);fp=nll(s+h);fm=nll(s-h);se=sqrt(1.0_dp/max((fp-2*f0+fm)/(h*h),1.0e-30_dp))
    end if
  contains
    real(dp) function nll(t) result(v)
      real(dp),intent(in)::t
      integer::i
      if(t<=0.0_dp)then;v=huge(1.0_dp);return;end if
      v=0.0_dp;do i=1,size(x);v=v-dzipfman(x(i),t,0.0_dp,nmax,.true.);end do
    end function nll
  end function zipf_mle

  subroutine zipfman_mle(x,nmax,s,b,cov)
    integer,intent(in)::x(:),nmax
    real(dp),intent(out)::s,b
    real(dp),intent(out),optional::cov(2,2)
    real(dp)::q(2),hess(2,2),invh(2,2)
    integer::info
    q=[log(1.0_dp),log(1.0_dp)]
    call nelder_mead(nll,q,step=0.2_dp,tol=1.0e-10_dp,max_iter=4000)
    s=exp(q(1));b=exp(q(2))
    if(present(cov))then
      call numerical_hessian(nll,q,hess);call invert_matrix(hess,invh,info)
      if(info==0)then
        cov(1,1)=invh(1,1)*s*s;cov(2,2)=invh(2,2)*b*b
        cov(1,2)=invh(1,2)*s*b;cov(2,1)=cov(1,2)
      else;cov=0.0_dp;end if
    end if
  contains
    real(dp) function nll(qv) result(v)
      real(dp),intent(in)::qv(:)
      real(dp)::ss,bb
      integer::i
      ss=exp(qv(1));bb=exp(qv(2));v=0.0_dp
      do i=1,size(x);v=v-dzipfman(x(i),ss,bb,nmax,.true.);end do
    end function nll
  end subroutine zipfman_mle

  real(dp) function zeta_mle(x,se) result(s)
    integer,intent(in)::x(:)
    real(dp),intent(out),optional::se
    real(dp)::hi,h,f0,fp,fm
    hi=3.0_dp;do while(nll(2.0_dp*hi)<nll(hi) .and. hi<100.0_dp);hi=2.0_dp*hi;end do
    s=golden_minimize(nll,1.0_dp+1.0e-10_dp,hi,1.0e-10_dp)
    if(present(se))then
      h=max(1.0e-5_dp,1.0e-4_dp*s);h=min(h,0.49_dp*(s-1.0_dp))
      f0=nll(s);fp=nll(s+h);fm=nll(s-h);se=sqrt(1.0_dp/max((fp-2*f0+fm)/(h*h),1.0e-30_dp))
    end if
  contains
    real(dp) function nll(t) result(v)
      real(dp),intent(in)::t
      integer::i
      if(t<=1.0_dp)then;v=huge(1.0_dp);return;end if
      v=0.0_dp;do i=1,size(x);v=v-dzipfman(x(i),t,0.0_dp,-1,.true.);end do
    end function nll
  end function zeta_mle

  function zipftol_int(x,nmax,m,alpha,p,side,dist,s_hat,b_hat) result(out)
    integer,intent(in)::x(:)
    integer,intent(in),optional::nmax,m
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::dist
    real(dp),intent(out),optional::s_hat,b_hat
    type(discrete_tolerance_interval)::out
    real(dp)::a,pp,s,b,se,cov(2,2),z,sl,su,bl,bu,ql,qu
    integer::nm,mm,n,sd
    character(len=12)::di
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side;di='Zipf';if(present(dist))di=adjustl(dist)
    n=size(x);nm=maxval(x);if(present(nmax))nm=nmax;mm=n;if(present(m))mm=m
    if(sd==2)then;a=a/2.0_dp;pp=(pp+1.0_dp)/2.0_dp;end if
    z=normal_quantile(1.0_dp-a)*sqrt(real(n,dp)/real(mm,dp));b=0.0_dp
    select case(trim(di))
    case('Zipf')
      s=zipf_mle(x,nm,se);sl=max(s-z*se,1.0e-14_dp);su=s+z*se
      ql=qzipfman(1.0_dp-pp,su,0.0_dp,nm);qu=qzipfman(pp,sl,0.0_dp,nm)
    case('Zipf-Man','ZipfMan')
      call zipfman_mle(x,nm,s,b,cov);sl=max(s-z*sqrt(abs(cov(1,1))),1.0e-14_dp)
      su=s+z*sqrt(abs(cov(1,1)));bl=max(b-z*sqrt(abs(cov(2,2))),0.0_dp);bu=b+z*sqrt(abs(cov(2,2)))
      ql=qzipfman(1.0_dp-pp,su,bl,nm);qu=qzipfman(pp,sl,bu,nm)
    case('Zeta')
      s=zeta_mle(x,se);sl=max(s-z*se,1.0_dp+1.0e-10_dp);su=s+z*se
      ql=qzipfman(1.0_dp-pp,su,0.0_dp,-1);qu=qzipfman(pp,sl,0.0_dp,-1)
    case default
      s=0.0_dp;ql=1.0_dp;qu=real(nm,dp)
    end select
    out%alpha=merge(2.0_dp*a,a,sd==2);out%p=merge(2.0_dp*pp-1.0_dp,pp,sd==2);out%estimate=s
    out%lower=max(1,int(ql));out%upper_infinite=qu>real(huge(out%upper)/2,dp)
    if(out%upper_infinite)then;out%upper=huge(out%upper);else;out%upper=max(out%lower,int(qu));end if
    if(present(s_hat))s_hat=s;if(present(b_hat))b_hat=b
  end function zipftol_int

end module tolerance_custom
