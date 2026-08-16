! SPDX-License-Identifier: GPL-2.0-only
module mbbefd_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mbbefd_kinds, only : dp
  use mbbefd_math, only : nan_dp, log_beta, beta_pdf, beta_cdf, beta_quantile, beta_random, &
    ncbeta_pdf, ncbeta_cdf, ncbeta_quantile, ncbeta_random, probability_to_lower, &
    gauss_legendre_integrate, digamma_dp
  use actuar, only : mbeta, levbeta, munif, levunif
  implicit none
  private

  public :: g2a, swiss_re
  public :: dmbbefd1, dmbbefd2, dmbbefd_gb1, dmbbefd_gb2, dgbeta1
  public :: doifun, poifun, qoifun, roifun, ecoifun, moifun, tloifun
  public :: dmbbefd, pmbbefd, qmbbefd, rmbbefd, ecmbbefd, mmbbefd, tlmbbefd
  public :: dmbbefd_gb, pmbbefd_gb, qmbbefd_gb, rmbbefd_gb, ecmbbefd_gb, mmbbefd_gb, tlmbbefd_gb
  public :: dstpareto, pstpareto, qstpareto, rstpareto, ecstpareto, mstpareto
  public :: dgbeta, pgbeta, qgbeta, rgbeta, ecgbeta, mgbeta
  public :: ecunif, ecbeta, theil_theoretical, theil_shape0
  public :: doiunif, poiunif, qoiunif, roiunif, ecoiunif, moiunif, tloiunif
  public :: doibeta, poibeta, qoibeta, roibeta, ecoibeta, moibeta, tloibeta
  public :: doistpareto, poistpareto, qoistpareto, roistpareto, ecoistpareto, moistpareto, tloistpareto
  public :: doigbeta, poigbeta, qoigbeta, roigbeta, ecoigbeta, moigbeta, tloigbeta
  public :: trans_m10, invt_m10, trans_01, invt_01, trans_0inf, invt_0inf, trans_1inf, invt_1inf
  public :: etl, theil_empirical, eecf_t, make_eecf
  public :: loglik_mbbefd, gradloglik_mbbefd, hessloglik_mbbefd
  public :: loglik_mbbefd_gb, gradloglik_mbbefd_gb

  abstract interface
    pure function oi_density_callback(x,par) result(v)
      import dp
      real(dp),intent(in)::x,par(:)
      real(dp)::v
    end function oi_density_callback
    pure function oi_cdf_callback(x,par) result(v)
      import dp
      real(dp),intent(in)::x,par(:)
      real(dp)::v
    end function oi_cdf_callback
    pure function oi_quantile_callback(p,par) result(v)
      import dp
      real(dp),intent(in)::p,par(:)
      real(dp)::v
    end function oi_quantile_callback
    subroutine oi_random_callback(x,par)
      import dp
      real(dp),intent(out)::x(:)
      real(dp),intent(in)::par(:)
    end subroutine oi_random_callback
    pure function oi_exposure_callback(x,par) result(v)
      import dp
      real(dp),intent(in)::x,par(:)
      real(dp)::v
    end function oi_exposure_callback
    function oi_moment_callback(order,par) result(v)
      import dp
      real(dp),intent(in)::order,par(:)
      real(dp)::v
    end function oi_moment_callback
  end interface

  type :: mb_moment_context
    integer :: family = 0
    real(dp) :: p1 = 0.0_dp
    real(dp) :: p2 = 0.0_dp
    real(dp) :: order = 1.0_dp
  end type mb_moment_context

  type, public :: eecf_t
    real(dp), allocatable :: x(:)
    real(dp) :: mean_x = 0.0_dp
  contains
    procedure :: evaluate => eecf_evaluate
  end type eecf_t

contains

  pure function dmbbefd1(x,a,b,log_value) result(v)
    real(dp),intent(in)::x,a,b;logical,intent(in),optional::log_value;real(dp)::v
    v=dmbbefd(x,trans_m10(a),trans_1inf(b),log_value)
  end function dmbbefd1

  pure function dmbbefd2(x,a,b,log_value) result(v)
    real(dp),intent(in)::x,a,b;logical,intent(in),optional::log_value;real(dp)::v
    v=dmbbefd(x,trans_0inf(a),trans_01(b),log_value)
  end function dmbbefd2

  pure function dmbbefd_gb1(x,g,b,log_value) result(v)
    real(dp),intent(in)::x,g,b;logical,intent(in),optional::log_value;real(dp)::v
    v=dmbbefd_gb(x,trans_1inf(g),trans_1inf(b),log_value)
  end function dmbbefd_gb1

  pure function dmbbefd_gb2(x,g,b,log_value) result(v)
    real(dp),intent(in)::x,g,b;logical,intent(in),optional::log_value;real(dp)::v
    v=dmbbefd_gb(x,trans_1inf(g),trans_01(b),log_value)
  end function dmbbefd_gb2

  pure function dgbeta1(x,shape0,shape1,shape2,log_value) result(v)
    real(dp),intent(in)::x,shape0,shape1,shape2;logical,intent(in),optional::log_value;real(dp)::v
    v=dgbeta(x,trans_0inf(shape0),trans_0inf(shape1),trans_0inf(shape2),log_value)
  end function dgbeta1

  pure function doifun(x,dfun,par,p1,log_value) result(v)
    real(dp),intent(in)::x,par(:),p1;procedure(oi_density_callback)::dfun
    logical,intent(in),optional::log_value;real(dp)::v;logical::l
    l=.false.;if(present(log_value))l=log_value
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    if(x==1.0_dp)then;v=p1;else;v=(1.0_dp-p1)*dfun(x,par);end if
    if(l)then;if(v>0.0_dp)then;v=log(v);else if(v==0.0_dp)then;v=-huge(1.0_dp);end if;end if
  end function doifun

  pure function poifun(q,pfun,par,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,par(:),p1;procedure(oi_cdf_callback)::pfun
    logical,intent(in),optional::lower_tail,log_p;real(dp)::v;logical::low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    v=(1.0_dp-p1)*pfun(q,par)+merge(p1,0.0_dp,q>=1.0_dp)
    if(.not.low)v=1.0_dp-v
    if(lp)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function poifun

  pure function qoifun(p,qfun,par,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::p,par(:),p1;procedure(oi_quantile_callback)::qfun
    logical,intent(in),optional::lower_tail,log_p;real(dp)::v,pp
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    pp=probability_to_lower(p,lower_tail,log_p);if(pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp();return;end if
    if(p1>=1.0_dp.or.pp>=1.0_dp-p1)then;v=1.0_dp;else;v=qfun(pp/(1.0_dp-p1),par);end if
  end function qoifun

  subroutine roifun(x,rfun,par,p1)
    real(dp),intent(out)::x(:);real(dp),intent(in)::par(:),p1;procedure(oi_random_callback)::rfun
    real(dp)::u;integer::i
    if(p1<0.0_dp.or.p1>1.0_dp)then;x=nan_dp();return;end if
    call rfun(x,par);do i=1,size(x);call random_number(u);if(u<p1)x(i)=1.0_dp;end do
  end subroutine roifun

  function ecoifun(x,ecfun,mfun,par,p1) result(v)
    real(dp),intent(in)::x,par(:),p1;procedure(oi_exposure_callback)::ecfun
    procedure(oi_moment_callback)::mfun;real(dp)::v,g0,e0
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    g0=ecfun(x,par);e0=mfun(1.0_dp,par);v=((1.0_dp-p1)*g0+p1*x/e0)/(1.0_dp-p1+p1/e0)
  end function ecoifun

  function moifun(order,mfun,par,p1) result(v)
    real(dp),intent(in)::order,par(:),p1;procedure(oi_moment_callback)::mfun;real(dp)::v
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();else;v=p1+(1.0_dp-p1)*mfun(order,par);end if
  end function moifun

  pure function tloifun(p1) result(v)
    real(dp),intent(in)::p1;real(dp)::v;v=merge(p1,nan_dp(),p1>=0.0_dp.and.p1<=1.0_dp)
  end function tloifun

  pure function g2a(g,b) result(a)
    real(dp), intent(in) :: g,b
    real(dp) :: a
    a=((g-1.0_dp)*b)/(1.0_dp-g*b)
  end function g2a

  pure function swiss_re(c) result(bg)
    real(dp), intent(in) :: c
    real(dp) :: bg(2)
    bg(1)=exp(3.1_dp-0.15_dp*c*(1.0_dp+c))
    bg(2)=exp(c*(0.78_dp+0.12_dp*c))
  end function swiss_re

  pure function dmbbefd(x,a,b,log_value) result(v)
    real(dp), intent(in) :: x,a,b
    logical, intent(in), optional :: log_value
    real(dp) :: v
    logical :: l
    l=.false.;if(present(log_value))l=log_value
    v=nan_dp()
    if ((b==0.0_dp .and. a>0.0_dp) .or. (a==-1.0_dp .and. b>1.0_dp)) then
      v=merge(1.0_dp,0.0_dp,x==1.0_dp)
    else if ((b==1.0_dp .and. a>-1.0_dp) .or. &
             ((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) .or. (a==0.0_dp.and.b>0.0_dp)) then
      v=merge(1.0_dp,0.0_dp,x==1.0_dp)
    else if ((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.ieee_is_finite(b)) then
      if(x>0.0_dp.and.x<1.0_dp) then; v=-log(b)*b**x
      else if(x==1.0_dp) then; v=1.0_dp-b
      else; v=0.0_dp; end if
    else if (((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
             .and.ieee_is_finite(a).and.ieee_is_finite(b)) then
      if(x>0.0_dp.and.x<1.0_dp) then
        v=-a*(a+1.0_dp)*b**x*log(b)/(a+b**x)**2
      else if(x==1.0_dp) then
        v=(a+1.0_dp)*b/(a+b)
      else
        v=0.0_dp
      end if
    end if
    if(l) then
      if(v>0.0_dp) then; v=log(v); else if(v==0.0_dp) then; v=-huge(1.0_dp); end if
    end if
  end function dmbbefd

  pure function pmbbefd(q,a,b,lower_tail,log_p) result(v)
    real(dp), intent(in) :: q,a,b
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: v
    logical :: low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    v=nan_dp()
    if ((b==0.0_dp.and.a>0.0_dp).or.(a==-1.0_dp.and.b>1.0_dp)) then
      v=merge(1.0_dp,0.0_dp,q>=1.0_dp)
    else if ((b==1.0_dp.and.a>-1.0_dp).or.((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) &
             .or.(a==0.0_dp.and.b>0.0_dp)) then
      v=merge(1.0_dp,0.0_dp,q>=1.0_dp)
    else if ((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.ieee_is_finite(b)) then
      if(q<=0.0_dp) then; v=0.0_dp; else if(q<1.0_dp) then; v=1.0_dp-b**q; else; v=1.0_dp; end if
    else if (((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
             .and.ieee_is_finite(a).and.ieee_is_finite(b)) then
      if(q<=0.0_dp) then; v=0.0_dp
      else if(q<1.0_dp) then; v=1.0_dp-(a+1.0_dp)*b**q/(a+b**q)
      else; v=1.0_dp; end if
    end if
    if(.not.low .and. ieee_is_finite(v)) v=1.0_dp-v
    if(lp .and. ieee_is_finite(v)) then
      if(v>0.0_dp) then; v=log(v); else; v=-huge(1.0_dp); end if
    end if
  end function pmbbefd

  pure function qmbbefd(p,a,b,lower_tail,log_p) result(x)
    real(dp), intent(in) :: p,a,b
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: x,pp,cut
    pp=probability_to_lower(p,lower_tail,log_p);x=nan_dp()
    if(pp<0.0_dp.or.pp>1.0_dp)return
    if((b==0.0_dp.and.a>0.0_dp).or.(a==-1.0_dp.and.b>1.0_dp)) then; x=1.0_dp
    else if((b==1.0_dp.and.a>-1.0_dp).or.((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) &
            .or.(a==0.0_dp.and.b>0.0_dp)) then; x=1.0_dp
    else if((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.ieee_is_finite(b)) then
      if(pp<1.0_dp-b) then; x=log(1.0_dp-pp)/log(b); else; x=1.0_dp; end if
    else if(((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
            .and.ieee_is_finite(a).and.ieee_is_finite(b)) then
      cut=1.0_dp-(a+1.0_dp)*b/(a+b)
      if(pp==0.0_dp) then; x=0.0_dp
      else if(pp<cut) then; x=log((1.0_dp-pp)*a/(a+pp))/log(b)
      else; x=1.0_dp; end if
    end if
  end function qmbbefd

  subroutine rmbbefd(x,a,b)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: a,b
    integer :: i
    real(dp) :: u
    do i=1,size(x);call random_number(u);x(i)=qmbbefd(u,a,b);end do
  end subroutine rmbbefd

  pure function ecmbbefd(x,a,b) result(v)
    real(dp), intent(in)::x,a,b
    real(dp)::v
    v=nan_dp();if(x<0.0_dp.or.x>1.0_dp)return
    if((b==0.0_dp.and.a>0.0_dp).or.(a==-1.0_dp.and.b>1.0_dp))then;v=merge(1.0_dp,0.0_dp,x==1.0_dp)
    else if((b==1.0_dp.and.a>-1.0_dp).or.((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) &
      .or.(a==0.0_dp.and.b>0.0_dp))then;v=x
    else if((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.ieee_is_finite(b))then
      v=(1.0_dp-b**x)/(1.0_dp-b)
    else if(((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
      .and.ieee_is_finite(a).and.ieee_is_finite(b))then
      v=log((a+b**x)/(a+1.0_dp))/log((a+b)/(a+1.0_dp))
    end if
  end function ecmbbefd

  function mmbbefd(order,a,b) result(v)
    real(dp),intent(in)::order,a,b
    real(dp)::v
    type(mb_moment_context)::ctx
    if(order<=0.0_dp)then;v=merge(1.0_dp,nan_dp(),order==0.0_dp);return;end if
    if((b==0.0_dp.and.a>0.0_dp).or.(a==-1.0_dp.and.b>1.0_dp).or. &
       (b==1.0_dp.and.a>-1.0_dp).or.((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) &
       .or.(a==0.0_dp.and.b>0.0_dp))then;v=1.0_dp;return;end if
    if(order==1.0_dp)then
      if((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)then;v=(b-1.0_dp)/log(b);return;end if
      if(((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
        .and.ieee_is_finite(a).and.ieee_is_finite(b))then
        v=log((a+b)/(a+1.0_dp))/log(b)*(a+1.0_dp);return
      end if
    end if
    if((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)then
      ctx=mb_moment_context(1,0.0_dp,b,order);v=gauss_legendre_integrate(moment_integrand,ctx,0.0_dp,1.0_dp,128)
    else if(((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
      .and.ieee_is_finite(a).and.ieee_is_finite(b))then
      ctx=mb_moment_context(2,a,b,order);v=gauss_legendre_integrate(moment_integrand,ctx,0.0_dp,1.0_dp,128)
    else;v=nan_dp();end if
  end function mmbbefd

  pure function tlmbbefd(a,b) result(v)
    real(dp),intent(in)::a,b;real(dp)::v
    v=nan_dp()
    if((b==0.0_dp.and.a>0.0_dp).or.(a==-1.0_dp.and.b>1.0_dp).or. &
       (b==1.0_dp.and.a>-1.0_dp).or.((.not.ieee_is_finite(b)).and.a>-1.0_dp.and.a<0.0_dp) &
       .or.(a==0.0_dp.and.b>0.0_dp))then;v=1.0_dp
    else if((.not.ieee_is_finite(a)).and.a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)then;v=b
    else if(((a>-1.0_dp.and.a<0.0_dp.and.b>1.0_dp).or.(a>0.0_dp.and.b>0.0_dp.and.b<1.0_dp)) &
      .and.ieee_is_finite(a).and.ieee_is_finite(b))then;v=(a+1.0_dp)*b/(a+b)
    end if
  end function tlmbbefd

  pure function dmbbefd_gb(x,g,b,log_value) result(v)
    real(dp),intent(in)::x,g,b;logical,intent(in),optional::log_value
    real(dp)::v,den;logical::l
    l=.false.;if(present(log_value))l=log_value;v=nan_dp()
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;v=merge(1.0_dp,0.0_dp,x==1.0_dp)
    else if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or.(g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then
      v=merge(1.0_dp,0.0_dp,x==1.0_dp)
    else if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp.and.ieee_is_finite(b))then
      if(x>0.0_dp.and.x<1.0_dp)then;v=-log(b)*b**x;else if(x==1.0_dp)then;v=1.0_dp-b;else;v=0.0_dp;end if
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then
      if(x>0.0_dp.and.x<1.0_dp)then;v=(g-1.0_dp)/(1.0_dp+(g-1.0_dp)*x)**2
      else if(x==1.0_dp)then;v=1.0_dp/g;else;v=0.0_dp;end if
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then
      if(x>0.0_dp.and.x<1.0_dp)then
        den=(g-1.0_dp)*b**(1.0_dp-x)+1.0_dp-g*b
        v=-(1.0_dp-b)*(g-1.0_dp)*log(b)*b**(1.0_dp-x)/(den*den)
      else if(x==1.0_dp)then;v=1.0_dp/g;else;v=0.0_dp;end if
    end if
    if(l)then;if(v>0.0_dp)then;v=log(v);else if(v==0.0_dp)then;v=-huge(1.0_dp);end if;end if
  end function dmbbefd_gb

  pure function pmbbefd_gb(q,g,b,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,g,b;logical,intent(in),optional::lower_tail,log_p
    real(dp)::v;logical::low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p;v=nan_dp()
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;v=merge(1.0_dp,0.0_dp,q>=1.0_dp)
    else if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or.(g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then
      v=merge(1.0_dp,0.0_dp,q>=1.0_dp)
    else if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp.and.ieee_is_finite(b))then
      if(q<=0.0_dp)then;v=0.0_dp;else if(q<1.0_dp)then;v=1.0_dp-b**q;else;v=1.0_dp;end if
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then
      if(q<=0.0_dp)then;v=0.0_dp;else if(q<1.0_dp)then;v=1.0_dp-1.0_dp/(1.0_dp+(g-1.0_dp)*q);else;v=1.0_dp;end if
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then
      if(q<=0.0_dp)then;v=0.0_dp
      else if(q<1.0_dp)then;v=1.0_dp-(1.0_dp-b)/((g-1.0_dp)*b**(1.0_dp-q)+1.0_dp-g*b)
      else;v=1.0_dp;end if
    end if
    if(.not.low.and.ieee_is_finite(v))v=1.0_dp-v
    if(lp.and.ieee_is_finite(v))then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function pmbbefd_gb

  pure function qmbbefd_gb(p,g,b,lower_tail,log_p) result(x)
    real(dp),intent(in)::p,g,b;logical,intent(in),optional::lower_tail,log_p
    real(dp)::x,pp,t
    pp=probability_to_lower(p,lower_tail,log_p);x=nan_dp();if(pp<0.0_dp.or.pp>1.0_dp)return
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;x=1.0_dp
    else if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or. &
      (g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then
      x=1.0_dp
    else if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp)then
      if(pp<1.0_dp-b)then;x=log(1.0_dp-pp)/log(b);else;x=1.0_dp;end if
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then
      if(pp<1.0_dp-1.0_dp/g)then;x=pp/((1.0_dp-pp)*(g-1.0_dp));else;x=1.0_dp;end if
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then
      if(pp==0.0_dp)then;x=0.0_dp
      else if(pp<1.0_dp-1.0_dp/g)then
        t=(g*b-1.0_dp)/(g-1.0_dp)+(1.0_dp-b)/((1.0_dp-pp)*(g-1.0_dp));x=1.0_dp-log(t)/log(b)
      else;x=1.0_dp;end if
    end if
  end function qmbbefd_gb

  subroutine rmbbefd_gb(x,g,b)
    real(dp),intent(out)::x(:);real(dp),intent(in)::g,b;integer::i;real(dp)::u
    do i=1,size(x);call random_number(u);x(i)=qmbbefd_gb(u,g,b);end do
  end subroutine rmbbefd_gb

  pure function ecmbbefd_gb(x,g,b) result(v)
    real(dp),intent(in)::x,g,b;real(dp)::v,t
    v=nan_dp();if(x<0.0_dp.or.x>1.0_dp)return
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;v=merge(1.0_dp,0.0_dp,x>0.0_dp)
    else if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or.(g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then;v=x
    else if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp)then;v=(1.0_dp-b**x)/(1.0_dp-b)
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then;v=log(1.0_dp+(g-1.0_dp)*x)/log(g)
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then
      t=(g-1.0_dp)*b+(1.0_dp-g*b)*b**x;v=log(t/(1.0_dp-b))/log(g*b)
    end if
  end function ecmbbefd_gb

  function mmbbefd_gb(order,g,b) result(v)
    real(dp),intent(in)::order,g,b;real(dp)::v;type(mb_moment_context)::ctx
    if(order<=0.0_dp)then;v=merge(1.0_dp,nan_dp(),order==0.0_dp);return;end if
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;v=1.0_dp;return;end if
    if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or. &
      (g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then
      v=1.0_dp;return
    end if
    if(order==1.0_dp)then
      if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp)then;v=(b-1.0_dp)/log(b);return;end if
      if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then;v=log(g)/(g-1.0_dp);return;end if
      if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then
        v=(1.0_dp-b)*log(g*b)/((1.0_dp-g*b)*log(b));return
      end if
    end if
    if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp)then;ctx=mb_moment_context(1,0.0_dp,b,order)
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then;ctx=mb_moment_context(3,g,b,order)
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and. &
      ieee_is_finite(g).and.ieee_is_finite(b))then
      ctx=mb_moment_context(4,g,b,order)
    else;v=nan_dp();return;end if
    v=gauss_legendre_integrate(moment_integrand,ctx,0.0_dp,1.0_dp,128)
  end function mmbbefd_gb

  pure function tlmbbefd_gb(g,b) result(v)
    real(dp),intent(in)::g,b;real(dp)::v
    v=nan_dp()
    if((.not.ieee_is_finite(g)).and.g>0.0_dp.and.b>0.0_dp.and.b/=1.0_dp)then;v=1.0_dp
    else if((g>1.0_dp.and.b==0.0_dp).or.(g>1.0_dp.and.(.not.ieee_is_finite(b))).or. &
      (g==1.0_dp.and.b>0.0_dp.and.b/=1.0_dp))then
      v=1.0_dp
    else if(g>1.0_dp.and.b>0.0_dp.and.b<1.0_dp.and.b*g==1.0_dp)then;v=b
    else if(g>1.0_dp.and.b==1.0_dp.and.ieee_is_finite(g))then;v=1.0_dp/g
    else if(g>1.0_dp.and.b>0.0_dp.and.b/=1.0_dp.and.b*g/=1.0_dp.and.ieee_is_finite(g).and.ieee_is_finite(b))then;v=1.0_dp/g
    end if
  end function tlmbbefd_gb

  function moment_integrand(x,context) result(v)
    real(dp),intent(in)::x;class(*),intent(in)::context;real(dp)::v,t
    select type(c=>context);type is(mb_moment_context)
      if(x==0.0_dp)then;t=0.0_dp;else;t=x**(1.0_dp/c%order);end if
      select case(c%family)
      case(1);v=c%p2**t
      case(2);v=(c%p1+1.0_dp)*c%p2**t/(c%p1+c%p2**t)
      case(3);v=1.0_dp/(1.0_dp+(c%p1-1.0_dp)*t)
      case(4);v=(1.0_dp-c%p2)/((c%p1-1.0_dp)*c%p2**(1.0_dp-t)+1.0_dp-c%p1*c%p2)
      case default;v=0.0_dp
      end select
    class default;v=0.0_dp;end select
  end function moment_integrand

  pure function dstpareto(x,a,log_value) result(v)
    real(dp),intent(in)::x,a;logical,intent(in),optional::log_value;real(dp)::v;logical::l
    l=.false.;if(present(log_value))l=log_value
    if(a<=0.0_dp)then;v=nan_dp();return;end if
    if(x<0.0_dp.or.x>1.0_dp)then;v=0.0_dp;else;v=a*(x+1.0_dp)**(-a-1.0_dp)/(1.0_dp-2.0_dp**(-a));end if
    if(l)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function dstpareto

  pure function pstpareto(q,a,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,a;logical,intent(in),optional::lower_tail,log_p;real(dp)::v;logical::low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    if(a<=0.0_dp)then;v=nan_dp();return;end if
    if(q<=0.0_dp)then;v=0.0_dp;else if(q>=1.0_dp)then;v=1.0_dp;else;v=(1.0_dp-(q+1.0_dp)**(-a))/(1.0_dp-2.0_dp**(-a));end if
    if(.not.low)v=1.0_dp-v;if(lp)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function pstpareto

  pure function qstpareto(p,a,lower_tail,log_p) result(x)
    real(dp),intent(in)::p,a;logical,intent(in),optional::lower_tail,log_p;real(dp)::x,pp
    pp=probability_to_lower(p,lower_tail,log_p)
    if(a<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;x=nan_dp();else;x=(1.0_dp-pp*(1.0_dp-2.0_dp**(-a)))**(-1.0_dp/a)-1.0_dp;end if
  end function qstpareto

  subroutine rstpareto(x,a)
    real(dp),intent(out)::x(:);real(dp),intent(in)::a;integer::i;real(dp)::u
    do i=1,size(x);call random_number(u);x(i)=qstpareto(u,a);end do
  end subroutine rstpareto

  pure function ecstpareto(x,a) result(v)
    real(dp),intent(in)::x,a;real(dp)::v
    if(a<=0.0_dp)then;v=nan_dp();return;end if
    if(x<=0.0_dp)then;v=0.0_dp;return;else if(x>=1.0_dp)then;v=1.0_dp;return;end if
    if(a==1.0_dp)then;v=(2.0_dp*log(x+1.0_dp)-x)/(2.0_dp*log(2.0_dp)-1.0_dp)
    else;v=((x+1.0_dp)**(1.0_dp-a)-2.0_dp**(-a)*x*(1.0_dp-a)-1.0_dp)/ &
      (2.0_dp**(1.0_dp-a)-2.0_dp**(-a)*(1.0_dp-a)-1.0_dp);end if
  end function ecstpareto

  function mstpareto(order,a) result(v)
    real(dp),intent(in)::order,a;real(dp)::v
    type(mb_moment_context)::ctx
    if(a<=0.0_dp.or.order<0.0_dp)then;v=nan_dp();return;end if
    if(order==0.0_dp)then;v=1.0_dp;return;end if
    if(order==1.0_dp)then
      if(a==1.0_dp)then;v=2.0_dp*log(2.0_dp)-1.0_dp
      else;v=(2.0_dp**(1.0_dp-a)-2.0_dp**(-a)*(1.0_dp-a)-1.0_dp)/((1.0_dp-a)*(1.0_dp-2.0_dp**(-a)));end if
      return
    end if
    ctx=mb_moment_context(5,a,0.0_dp,order);v=gauss_legendre_integrate(stpareto_moment_integrand,ctx,0.0_dp,1.0_dp,128)
  end function mstpareto

  function stpareto_moment_integrand(x,context) result(v)
    real(dp),intent(in)::x;class(*),intent(in)::context;real(dp)::v
    select type(c=>context);type is(mb_moment_context);v=x**c%order*dstpareto(x,c%p1);class default;v=0.0_dp;end select
  end function stpareto_moment_integrand

  recursive pure function dgbeta(x,shape0,shape1,shape2,log_value) result(v)
    real(dp),intent(in)::x,shape0,shape1,shape2;logical,intent(in),optional::log_value
    real(dp)::v,lv;logical::l
    l=.false.;if(present(log_value))l=log_value
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp)then;v=nan_dp();return;end if
    if(x<0.0_dp.or.x>1.0_dp)then;v=merge(-huge(1.0_dp),0.0_dp,l);return;end if
    if(l)then
      if(x==0.0_dp)then
        if(shape0*shape1<1.0_dp)then
          v=huge(1.0_dp)
        else if(shape0*shape1==1.0_dp)then
          v=log(shape0)-log_beta(shape1,shape2)
        else
          v=-huge(1.0_dp)
        end if
      else;lv=log(shape0)+(shape0-1.0_dp)*log(x)+beta_pdf(x**shape0,shape1,shape2,.true.);v=lv;end if
    else
      if(x==0.0_dp.and.shape0<1.0_dp)then
        ! evaluate through log form to avoid 0**negative before beta cancellation
        lv=dgbeta(x,shape0,shape1,shape2,.true.);if(lv>=log(huge(1.0_dp)))then;v=huge(1.0_dp);else;v=exp(lv);end if
      else;v=shape0*x**(shape0-1.0_dp)*beta_pdf(x**shape0,shape1,shape2);end if
    end if
  end function dgbeta

  pure function pgbeta(q,shape0,shape1,shape2,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,shape0,shape1,shape2;logical,intent(in),optional::lower_tail,log_p
    real(dp)::v;logical::low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp)then;v=nan_dp();return;end if
    if(q<=0.0_dp)then;v=0.0_dp;else if(q>=1.0_dp)then;v=1.0_dp;else;v=beta_cdf(q**shape0,shape1,shape2);end if
    if(.not.low)v=1.0_dp-v;if(lp)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function pgbeta

  pure function qgbeta(p,shape0,shape1,shape2,lower_tail,log_p) result(x)
    real(dp),intent(in)::p,shape0,shape1,shape2;logical,intent(in),optional::lower_tail,log_p
    real(dp)::x,pp
    pp=probability_to_lower(p,lower_tail,log_p)
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;x=nan_dp();return;end if
    x=beta_quantile(pp,shape1,shape2)**(1.0_dp/shape0)
  end function qgbeta

  subroutine rgbeta(x,shape0,shape1,shape2)
    real(dp),intent(out)::x(:);real(dp),intent(in)::shape0,shape1,shape2;integer::i
    do i=1,size(x);x(i)=beta_random(shape1,shape2)**(1.0_dp/shape0);end do
  end subroutine rgbeta

  pure function mgbeta(order,shape0,shape1,shape2) result(v)
    real(dp),intent(in)::order,shape0,shape1,shape2;real(dp)::v,k
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp.or.order<0.0_dp)then;v=nan_dp();return;end if
    if(order==0.0_dp)then;v=1.0_dp;return;end if
    k=order/shape0;v=exp(log_beta(shape1+k,shape2)-log_beta(shape1,shape2))
  end function mgbeta

  pure function ecgbeta(x,shape0,shape1,shape2) result(v)
    real(dp),intent(in)::x,shape0,shape1,shape2;real(dp)::v,cst2,xx
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp)then;v=nan_dp();return;end if
    if(x<=0.0_dp)then;v=0.0_dp;return;else if(x>=1.0_dp)then;v=1.0_dp;return;end if
    xx=x**shape0;cst2=exp(log_beta(shape1,1.0_dp/shape0)-log_beta(shape1+shape2,1.0_dp/shape0))
    v=beta_cdf(xx,shape1+1.0_dp/shape0,shape2)+x*(1.0_dp-beta_cdf(xx,shape1,shape2))*cst2
  end function ecgbeta

  pure function theil_theoretical(shape0,shape1,shape2) result(v)
    real(dp),intent(in)::shape0,shape1,shape2;real(dp)::v,ex
    if(shape0<=0.0_dp.or.shape1<=0.0_dp.or.shape2<=0.0_dp)then;v=nan_dp();return;end if
    ex=mgbeta(1.0_dp,shape0,shape1,shape2)
    v=(digamma_dp(shape1+1.0_dp/shape0)-digamma_dp(shape1+shape2+1.0_dp/shape0))/shape0-log(ex)
  end function theil_theoretical

  function theil_shape0(shape0,obs) result(v)
    real(dp),intent(in)::shape0,obs(:);real(dp)::v,m,var,aux,s1,s2;real(dp),allocatable::z(:)
    integer::n
    if(shape0<=0.0_dp.or.size(obs)<2)then;v=nan_dp();return;end if
    allocate(z(size(obs)));z=obs**shape0;n=size(z);m=sum(z)/real(n,dp);var=sum((z-m)**2)/real(n,dp)
    if(var<=0.0_dp)then;v=nan_dp();return;end if
    aux=m*(1.0_dp-m)/var-1.0_dp;s1=m*aux;s2=(1.0_dp-m)*aux;v=theil_theoretical(shape0,s1,s2)
  end function theil_shape0

  pure function ecunif(x,xmin,xmax) result(v)
    real(dp),intent(in)::x;real(dp),intent(in),optional::xmin,xmax;real(dp)::v,lo,hi
    lo=0.0_dp;hi=1.0_dp;if(present(xmin))lo=xmin;if(present(xmax))hi=xmax
    v=levunif(x,lo,hi,1.0_dp)/munif(1.0_dp,lo,hi)
  end function ecunif

  pure function ecbeta(x,shape1,shape2) result(v)
    real(dp),intent(in)::x,shape1,shape2;real(dp)::v
    v=levbeta(x,shape1,shape2,1.0_dp)/mbeta(1.0_dp,shape1,shape2)
  end function ecbeta

  pure function oi_density(base_density,x,p1,log_value) result(v)
    real(dp),intent(in)::base_density,x,p1;logical,intent(in),optional::log_value;real(dp)::v;logical::l
    l=.false.;if(present(log_value))l=log_value
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    if(x==1.0_dp)then;v=p1;else;v=(1.0_dp-p1)*base_density;end if
    if(l)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function oi_density

  pure function oi_cdf(base_cdf,q,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::base_cdf,q,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v;logical::low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    if(p1<0.0_dp.or.p1>1.0_dp)then;v=nan_dp();return;end if
    v=(1.0_dp-p1)*base_cdf+merge(p1,0.0_dp,q>=1.0_dp)
    if(.not.low)v=1.0_dp-v;if(lp)then;if(v>0.0_dp)then;v=log(v);else;v=-huge(1.0_dp);end if;end if
  end function oi_cdf

  pure function doiunif(x,p1,log_value) result(v)
    real(dp),intent(in)::x,p1;logical,intent(in),optional::log_value;real(dp)::v,d
    d=merge(1.0_dp,0.0_dp,x>=0.0_dp.and.x<=1.0_dp);v=oi_density(d,x,p1,log_value)
  end function doiunif
  pure function poiunif(q,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v,f
    f=min(1.0_dp,max(0.0_dp,q));v=oi_cdf(f,q,p1,lower_tail,log_p)
  end function poiunif
  pure function qoiunif(p,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::p,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v,pp
    pp=probability_to_lower(p,lower_tail,log_p);if(p1<0.0_dp.or.p1>1.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp()
    else if(pp>=1.0_dp-p1.or.p1==1.0_dp)then;v=1.0_dp;else;v=pp/(1.0_dp-p1);end if
  end function qoiunif
  subroutine roiunif(x,p1)
    real(dp),intent(out)::x(:);real(dp),intent(in)::p1;integer::i;real(dp)::u
    do i=1,size(x);call random_number(u);x(i)=qoiunif(u,p1);end do
  end subroutine roiunif
  pure function moiunif(order,p1) result(v)
    real(dp),intent(in)::order,p1;real(dp)::v;v=p1+(1.0_dp-p1)*munif(order,0.0_dp,1.0_dp)
  end function moiunif
  pure function ecoiunif(x,p1) result(v)
    real(dp),intent(in)::x,p1;real(dp)::v,e0,g0
    e0=munif(1.0_dp,0.0_dp,1.0_dp);g0=ecunif(x)
    v=((1.0_dp-p1)*g0+p1*x/e0)/(1.0_dp-p1+p1/e0)
  end function ecoiunif
  pure function tloiunif(p1) result(v)
    real(dp),intent(in)::p1;real(dp)::v
    v=merge(p1,nan_dp(),p1>=0.0_dp.and.p1<=1.0_dp)
  end function tloiunif

  function doibeta(x,shape1,shape2,p1,ncp,log_value) result(v)
    real(dp),intent(in)::x,shape1,shape2,p1;real(dp),intent(in),optional::ncp;logical,intent(in),optional::log_value
    real(dp)::v,nc,d;nc=0.0_dp;if(present(ncp))nc=ncp;d=ncbeta_pdf(x,shape1,shape2,nc);v=oi_density(d,x,p1,log_value)
  end function doibeta
  function poibeta(q,shape1,shape2,p1,ncp,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,shape1,shape2,p1;real(dp),intent(in),optional::ncp;logical,intent(in),optional::lower_tail,log_p
    real(dp)::v,nc,f;nc=0.0_dp;if(present(ncp))nc=ncp;f=ncbeta_cdf(q,shape1,shape2,nc);v=oi_cdf(f,q,p1,lower_tail,log_p)
  end function poibeta
  function qoibeta(p,shape1,shape2,p1,ncp,lower_tail,log_p) result(v)
    real(dp),intent(in)::p,shape1,shape2,p1;real(dp),intent(in),optional::ncp;logical,intent(in),optional::lower_tail,log_p
    real(dp)::v,nc,pp;nc=0.0_dp;if(present(ncp))nc=ncp;pp=probability_to_lower(p,lower_tail,log_p)
    if(p1<0.0_dp.or.p1>1.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp();else if(pp>=1.0_dp-p1.or.p1==1.0_dp)then;v=1.0_dp
    else;v=ncbeta_quantile(pp/(1.0_dp-p1),shape1,shape2,nc);end if
  end function qoibeta
  subroutine roibeta(x,shape1,shape2,p1,ncp)
    real(dp),intent(out)::x(:);real(dp),intent(in)::shape1,shape2,p1;real(dp),intent(in),optional::ncp;integer::i;real(dp)::u,nc
    nc=0.0_dp;if(present(ncp))nc=ncp
    do i=1,size(x);call random_number(u);if(u<p1)then;x(i)=1.0_dp;else;x(i)=ncbeta_random(shape1,shape2,nc);end if;end do
  end subroutine roibeta
  pure function moibeta(order,shape1,shape2,p1,ncp) result(v)
    real(dp),intent(in)::order,shape1,shape2,p1;real(dp),intent(in),optional::ncp;real(dp)::v,nc
    nc=0.0_dp;if(present(ncp))nc=ncp;if(nc/=0.0_dp)then;v=nan_dp();else;v=p1+(1.0_dp-p1)*mbeta(order,shape1,shape2);end if
  end function moibeta
  pure function ecoibeta(x,shape1,shape2,p1,ncp) result(v)
    real(dp),intent(in)::x,shape1,shape2,p1;real(dp),intent(in),optional::ncp;real(dp)::v,nc,e0,g0
    nc=0.0_dp;if(present(ncp))nc=ncp;if(nc/=0.0_dp)then;v=nan_dp();return;end if
    e0=mbeta(1.0_dp,shape1,shape2);g0=ecbeta(x,shape1,shape2);v=((1.0_dp-p1)*g0+p1*x/e0)/(1.0_dp-p1+p1/e0)
  end function ecoibeta
  pure function tloibeta(shape1,shape2,p1,ncp) result(v)
    real(dp),intent(in)::shape1,shape2,p1;real(dp),intent(in),optional::ncp;real(dp)::v
    v=merge(p1,nan_dp(),p1>=0.0_dp.and.p1<=1.0_dp)
  end function tloibeta

  pure function doistpareto(x,a,p1,log_value) result(v)
    real(dp),intent(in)::x,a,p1;logical,intent(in),optional::log_value;real(dp)::v;v=oi_density(dstpareto(x,a),x,p1,log_value)
  end function doistpareto
  pure function poistpareto(q,a,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,a,p1
    logical,intent(in),optional::lower_tail,log_p;real(dp)::v
    v=oi_cdf(pstpareto(q,a),q,p1,lower_tail,log_p)
  end function poistpareto
  pure function qoistpareto(p,a,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::p,a,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v,pp
    pp=probability_to_lower(p,lower_tail,log_p);if(p1<0.0_dp.or.p1>1.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp()
    else if(pp>=1.0_dp-p1.or.p1==1.0_dp)then;v=1.0_dp;else;v=qstpareto(pp/(1.0_dp-p1),a);end if
  end function qoistpareto
  subroutine roistpareto(x,a,p1)
    real(dp),intent(out)::x(:);real(dp),intent(in)::a,p1;integer::i;real(dp)::u
    do i=1,size(x);call random_number(u);x(i)=qoistpareto(u,a,p1);end do
  end subroutine roistpareto
  function moistpareto(order,a,p1) result(v)
    real(dp),intent(in)::order,a,p1;real(dp)::v
    v=p1+(1.0_dp-p1)*mstpareto(order,a)
  end function moistpareto
  function ecoistpareto(x,a,p1) result(v)
    real(dp),intent(in)::x,a,p1;real(dp)::v,e0,g0
    e0=mstpareto(1.0_dp,a);g0=ecstpareto(x,a)
    v=((1.0_dp-p1)*g0+p1*x/e0)/(1.0_dp-p1+p1/e0)
  end function ecoistpareto
  pure function tloistpareto(a,p1) result(v)
    real(dp),intent(in)::a,p1;real(dp)::v
    v=merge(p1,nan_dp(),p1>=0.0_dp.and.p1<=1.0_dp)
  end function tloistpareto

  pure function doigbeta(x,shape0,shape1,shape2,p1,log_value) result(v)
    real(dp),intent(in)::x,shape0,shape1,shape2,p1;logical,intent(in),optional::log_value;real(dp)::v
    v=oi_density(dgbeta(x,shape0,shape1,shape2),x,p1,log_value)
  end function doigbeta
  pure function poigbeta(q,shape0,shape1,shape2,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,shape0,shape1,shape2,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v
    v=oi_cdf(pgbeta(q,shape0,shape1,shape2),q,p1,lower_tail,log_p)
  end function poigbeta
  pure function qoigbeta(p,shape0,shape1,shape2,p1,lower_tail,log_p) result(v)
    real(dp),intent(in)::p,shape0,shape1,shape2,p1;logical,intent(in),optional::lower_tail,log_p;real(dp)::v,pp
    pp=probability_to_lower(p,lower_tail,log_p);if(p1<0.0_dp.or.p1>1.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then;v=nan_dp()
    else if(pp>=1.0_dp-p1.or.p1==1.0_dp)then;v=1.0_dp;else;v=qgbeta(pp/(1.0_dp-p1),shape0,shape1,shape2);end if
  end function qoigbeta
  subroutine roigbeta(x,shape0,shape1,shape2,p1)
    real(dp),intent(out)::x(:);real(dp),intent(in)::shape0,shape1,shape2,p1;integer::i;real(dp)::u
    do i=1,size(x);call random_number(u);x(i)=qoigbeta(u,shape0,shape1,shape2,p1);end do
  end subroutine roigbeta
  pure function moigbeta(order,shape0,shape1,shape2,p1) result(v)
    real(dp),intent(in)::order,shape0,shape1,shape2,p1;real(dp)::v;v=p1+(1.0_dp-p1)*mgbeta(order,shape0,shape1,shape2)
  end function moigbeta
  pure function ecoigbeta(x,shape0,shape1,shape2,p1) result(v)
    real(dp),intent(in)::x,shape0,shape1,shape2,p1;real(dp)::v,e0,g0
    e0=mgbeta(1.0_dp,shape0,shape1,shape2)
    g0=ecgbeta(x,shape0,shape1,shape2)
    v=((1.0_dp-p1)*g0+p1*x/e0)/(1.0_dp-p1+p1/e0)
  end function ecoigbeta
  pure function tloigbeta(shape0,shape1,shape2,p1) result(v)
    real(dp),intent(in)::shape0,shape1,shape2,p1;real(dp)::v;v=merge(p1,nan_dp(),p1>=0.0_dp.and.p1<=1.0_dp)
  end function tloigbeta

  pure function trans_m10(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=-1.0_dp/(1.0_dp+exp(-x));end function trans_m10
  pure function invt_m10(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=log(-x/(1.0_dp+x));end function invt_m10
  pure function trans_01(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=1.0_dp/(1.0_dp+exp(-x));end function trans_01
  pure function invt_01(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=log(x/(1.0_dp-x));end function invt_01
  pure function trans_0inf(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=exp(x);end function trans_0inf
  pure function invt_0inf(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=log(x);end function invt_0inf
  pure function trans_1inf(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=1.0_dp+exp(x);end function trans_1inf
  pure function invt_1inf(x) result(v);real(dp),intent(in)::x;real(dp)::v;v=log(x-1.0_dp);end function invt_1inf

  pure function etl(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::v
    if(size(x)==0)then;v=nan_dp();else;v=real(count(x==1.0_dp),dp)/real(size(x),dp);end if
  end function etl

  pure function theil_empirical(x) result(v)
    real(dp),intent(in)::x(:);real(dp)::v,m,r;integer::i
    if(size(x)==0)then;v=nan_dp();return;end if;m=sum(x)/real(size(x),dp);if(m<=0.0_dp)then;v=nan_dp();return;end if
    r=0.0_dp;do i=1,size(x);if(x(i)>0.0_dp)r=r+(x(i)/m)*log(x(i)/m);end do;v=r/real(size(x),dp)
  end function theil_empirical

  function make_eecf(data) result(obj)
    real(dp),intent(in)::data(:);type(eecf_t)::obj;integer::i,j;real(dp)::tmp
    allocate(obj%x(size(data)));obj%x=data
    do i=2,size(obj%x)
      tmp=obj%x(i);j=i-1
      do while(j>=1)
        if(obj%x(j)<=tmp)exit
        obj%x(j+1)=obj%x(j);j=j-1
      end do
      obj%x(j+1)=tmp
    end do
    if(size(data)>0)obj%mean_x=sum(data)/real(size(data),dp)
  end function make_eecf

  pure function eecf_evaluate(self,d) result(v)
    class(eecf_t),intent(in)::self;real(dp),intent(in)::d;real(dp)::v
    if(size(self%x)==0.or.self%mean_x<=0.0_dp)then;v=nan_dp();else;v=sum(min(self%x,d))/(real(size(self%x),dp)*self%mean_x);end if
  end function eecf_evaluate

  function loglik_mbbefd(obs,a,b) result(v)
    real(dp),intent(in)::obs(:),a,b;real(dp)::v,d;integer::i
    v=0.0_dp
    do i=1,size(obs)
      d=dmbbefd(obs(i),a,b)
      if(d<=0.0_dp.or..not.ieee_is_finite(d))then;v=-huge(1.0_dp);return;end if
      v=v+log(d)
    end do
  end function loglik_mbbefd

  function gradloglik_mbbefd(obs,a,b) result(gv)
    real(dp),intent(in)::obs(:),a,b;real(dp)::gv(2),x;integer::i
    gv=0.0_dp;do i=1,size(obs);x=obs(i)
      if(x==1.0_dp)then
        gv(1)=gv(1)+(b-1.0_dp)/((a+1.0_dp)*(a+b));gv(2)=gv(2)+a/(b*(a+b))
      else
        gv(1)=gv(1)+(2.0_dp*a+1.0_dp)/(a*(a+1.0_dp))-2.0_dp/(a+b**x)
        gv(2)=gv(2)+x/b+1.0_dp/(b*log(b))-2.0_dp*b**x*x/(b*(a+b**x))
      end if
    end do
  end function gradloglik_mbbefd

  function hessloglik_mbbefd(obs,a,b) result(h)
    real(dp),intent(in)::obs(:),a,b;real(dp)::h(2,2),x,h11,h21,h22;integer::i
    h=0.0_dp;do i=1,size(obs);x=obs(i)
      if(x==1.0_dp)then
        h11=1.0_dp/(a+b)**2-1.0_dp/(a+1.0_dp)**2;h21=1.0_dp/(a+b)**2;h22=1.0_dp/(a+b)**2-1.0_dp/b**2
      else
        h11=2.0_dp/(a+b**x)**2-1.0_dp/a**2-1.0_dp/(a+1.0_dp)**2
        h21=2.0_dp*x*b**(x-1.0_dp)/(a+b**x)**2
        h22=x/b**2-(log(b)+1.0_dp)/(b**2*log(b)**2)-2.0_dp*a*x/(b**2*(a+b**x))- &
          2.0_dp*a*x*x*b**x/(b**2*(a+b**x)**2)
      end if
      h(1,1)=h(1,1)+h11;h(1,2)=h(1,2)+h21;h(2,1)=h(2,1)+h21;h(2,2)=h(2,2)+h22
    end do
  end function hessloglik_mbbefd

  function loglik_mbbefd_gb(obs,g,b) result(v)
    real(dp),intent(in)::obs(:),g,b;real(dp)::v,d;integer::i
    v=0.0_dp
    do i=1,size(obs)
      d=dmbbefd_gb(obs(i),g,b)
      if(d<=0.0_dp.or..not.ieee_is_finite(d))then;v=-huge(1.0_dp);return;end if
      v=v+log(d)
    end do
  end function loglik_mbbefd_gb

  function gradloglik_mbbefd_gb(obs,g,b) result(gv)
    real(dp),intent(in)::obs(:),g,b;real(dp)::gv(2),x,den;integer::i
    gv=0.0_dp;do i=1,size(obs);x=obs(i)
      if(x==1.0_dp)then;gv(1)=gv(1)-1.0_dp/g
      else
        den=(g-1.0_dp)*b**(1.0_dp-x)+1.0_dp-g*b
        gv(1)=gv(1)+1.0_dp/(g-1.0_dp)-2.0_dp*(b**(1.0_dp-x)-b)/den
        gv(2)=gv(2)+1.0_dp/(b-1.0_dp)+1.0_dp/(b*log(b))+(1.0_dp-x)/b+ &
          ((g-1.0_dp)*(1.0_dp-x)*b**(-x)-g)/den
      end if
    end do
  end function gradloglik_mbbefd_gb

end module mbbefd_distributions
