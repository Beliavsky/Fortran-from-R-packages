! SPDX-License-Identifier: GPL-3.0-only
module tscopula_margins
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use tscopula_kinds, only : dp
  use tscopula_status, only : tsc_error, tsc_success, tsc_invalid_input, clear_error, set_error
  use tscopula_math, only : normal_pdf, normal_cdf, normal_quantile, student_pdf, &
    student_cdf, student_quantile, uniform_random, normal_random, student_random, &
    optimizer_result, minimize_nelder_mead, finite_hessian, safe_standard_errors
  implicit none
  private

  integer, parameter, public :: margin_uniform = 0
  integer, parameter, public :: margin_gauss = 1
  integer, parameter, public :: margin_gauss0 = 2
  integer, parameter, public :: margin_laplace = 3
  integer, parameter, public :: margin_laplace0 = 4
  integer, parameter, public :: margin_slaplace = 5
  integer, parameter, public :: margin_doubleweibull = 6
  integer, parameter, public :: margin_sdoubleweibull = 7
  integer, parameter, public :: margin_student = 8
  integer, parameter, public :: margin_student0 = 9
  integer, parameter, public :: margin_sstudent = 10
  integer, parameter, public :: margin_empirical = 11

  type, public :: margin_spec
    integer :: family = margin_uniform
    real(dp) :: mu = 0.0_dp
    real(dp) :: sigma = 1.0_dp
    real(dp) :: shape = 1.0_dp
    real(dp) :: gamma = 1.0_dp
    real(dp) :: df = 10.0_dp
  end type margin_spec

  type, public :: margin_fit_result
    type(margin_spec) :: margin
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: standard_errors(:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: convergence = 1
    integer :: iterations = 0
  end type margin_fit_result

  public :: margin, edf, dmarg, pmarg, qmarg, rmarg, fit_margin
  public :: dgauss, pgauss, qgauss, rgauss
  public :: dgauss0, pgauss0, qgauss0, rgauss0
  public :: dlaplace, plaplace, qlaplace, rlaplace
  public :: dlaplace0, plaplace0, qlaplace0, rlaplace0
  public :: dslaplace, pslaplace, qslaplace, rslaplace
  public :: ddoubleweibull, pdoubleweibull, qdoubleweibull, rdoubleweibull
  public :: dsdoubleweibull, psdoubleweibull, qsdoubleweibull, rsdoubleweibull
  public :: dst, pst, qst, rst, dst0, pst0, qst0, rst0
  public :: dsst, psst, qsst, rsst

  type :: margin_context
    type(margin_spec) :: template
    real(dp), allocatable :: data(:)
  end type margin_context

contains

  function margin(name, pars) result(spec)
    character(len=*), intent(in) :: name
    real(dp), intent(in), optional :: pars(:)
    type(margin_spec) :: spec
    character(len=:), allocatable :: lname
    lname = lower_string(trim(name))
    select case(lname)
    case('unif','uniform'); spec%family = margin_uniform
    case('gauss','norm','normal'); spec%family = margin_gauss
    case('gauss0','norm0'); spec%family = margin_gauss0
    case('laplace'); spec%family = margin_laplace
    case('laplace0'); spec%family = margin_laplace0
    case('slaplace'); spec%family = margin_slaplace
    case('doubleweibull'); spec%family = margin_doubleweibull
    case('sdoubleweibull'); spec%family = margin_sdoubleweibull
    case('st','student','t'); spec%family = margin_student
    case('st0','student0','t0'); spec%family = margin_student0
    case('sst'); spec%family = margin_sstudent
    case('edf'); spec%family = margin_empirical
    case default; spec%family = -1
    end select
    if (present(pars)) call assign_parameters(spec, pars)
  end function margin

  function edf() result(spec)
    type(margin_spec) :: spec
    spec%family = margin_empirical
  end function edf

  subroutine assign_parameters(spec, pars)
    type(margin_spec), intent(inout) :: spec
    real(dp), intent(in) :: pars(:)
    select case(spec%family)
    case(margin_gauss)
      if(size(pars)>=1)spec%mu=pars(1);if(size(pars)>=2)spec%sigma=pars(2)
    case(margin_gauss0)
      if(size(pars)>=1)spec%sigma=pars(1)
    case(margin_laplace)
      if(size(pars)>=1)spec%mu=pars(1);if(size(pars)>=2)spec%sigma=pars(2)
    case(margin_laplace0)
      if(size(pars)>=1)spec%sigma=pars(1)
    case(margin_slaplace)
      if(size(pars)>=1)spec%mu=pars(1);if(size(pars)>=2)spec%sigma=pars(2);if(size(pars)>=3)spec%gamma=pars(3)
    case(margin_doubleweibull)
      if(size(pars)>=1)spec%mu=pars(1);if(size(pars)>=2)spec%shape=pars(2);if(size(pars)>=3)spec%sigma=pars(3)
    case(margin_sdoubleweibull)
      if(size(pars)>=1)spec%mu=pars(1);if(size(pars)>=2)spec%shape=pars(2);if(size(pars)>=3)spec%sigma=pars(3);if(size(pars)>=4)spec%gamma=pars(4)
    case(margin_student)
      if(size(pars)>=1)spec%df=pars(1);if(size(pars)>=2)spec%mu=pars(2);if(size(pars)>=3)spec%sigma=pars(3)
    case(margin_student0)
      if(size(pars)>=1)spec%df=pars(1);if(size(pars)>=2)spec%sigma=pars(2)
    case(margin_sstudent)
      if(size(pars)>=1)spec%df=pars(1);if(size(pars)>=2)spec%gamma=pars(2);if(size(pars)>=3)spec%mu=pars(3);if(size(pars)>=4)spec%sigma=pars(4)
    end select
  end subroutine assign_parameters

  pure function lower_string(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i,c
    do i=1,len(s)
      c=iachar(s(i:i));if(c>=iachar('A').and.c<=iachar('Z'))then;out(i:i)=achar(c+32);else;out(i:i)=s(i:i);end if
    end do
  end function lower_string

  elemental real(dp) function dgauss(x,mu,sigma,log_density) result(value)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::mu,sigma
    logical,intent(in),optional::log_density
    real(dp)::m,s;logical::ld
    m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(sigma))s=sigma;ld=.false.;if(present(log_density))ld=log_density
    if(s<0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(s<=tiny(1.0_dp))then;value=merge(huge(1.0_dp),-huge(1.0_dp),.not.ld);return;end if
    value=-0.5_dp*((x-m)/s)**2-log(s)-0.5_dp*log(2.0_dp*acos(-1.0_dp))
    if(.not.ld)value=exp(value)
  end function dgauss

  elemental real(dp) function pgauss(q,mu,sigma) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,sigma
    real(dp)::m,s;m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(sigma))s=sigma
    if(s<0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;value=normal_cdf((q-m)/s);end if
  end function pgauss

  elemental real(dp) function qgauss(p,mu,sigma) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::mu,sigma
    real(dp)::m,s;m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(sigma))s=sigma
    if(s<0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;value=m+s*normal_quantile(p);end if
  end function qgauss

  function rgauss(n,mu,sigma) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,sigma
    real(dp),allocatable::values(:);real(dp)::m,s;integer::i
    m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(sigma))s=sigma;allocate(values(n))
    do i=1,n;values(i)=m+s*normal_random();end do
  end function rgauss

  elemental real(dp) function dgauss0(x,sigma,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::sigma;logical,intent(in),optional::log_density
    value=dgauss(x,0.0_dp,sigma,log_density)
  end function dgauss0
  elemental real(dp) function pgauss0(q,sigma) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::sigma;value=pgauss(q,0.0_dp,sigma)
  end function pgauss0
  elemental real(dp) function qgauss0(p,sigma) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::sigma;value=qgauss(p,0.0_dp,sigma)
  end function qgauss0
  function rgauss0(n,sigma) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::sigma;real(dp),allocatable::values(:);values=rgauss(n,0.0_dp,sigma)
  end function rgauss0

  elemental real(dp) function dlaplace(x,mu,scale,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,scale;logical,intent(in),optional::log_density
    real(dp)::m,s;logical::ld;m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(scale))s=scale;ld=.false.;if(present(log_density))ld=log_density
    if(s<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    value=-abs((x-m)/s)-log(2.0_dp*s);if(.not.ld)value=exp(value)
  end function dlaplace
  elemental real(dp) function plaplace(q,mu,scale) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,scale;real(dp)::m,s,y
    m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(scale))s=scale
    if(s<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    y=(q-m)/s;if(y<0.0_dp)then;value=0.5_dp*exp(y);else;value=1.0_dp-0.5_dp*exp(-y);end if
  end function plaplace
  elemental real(dp) function qlaplace(p,mu,scale) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::mu,scale;real(dp)::m,s
    m=0.0_dp;if(present(mu))m=mu;s=1.0_dp;if(present(scale))s=scale
    if(s<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else if(p<0.5_dp)then;value=m+s*log(2.0_dp*p);else;value=m-s*log(2.0_dp-2.0_dp*p);end if
  end function qlaplace
  function rlaplace(n,mu,scale) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,scale;real(dp),allocatable::values(:);integer::i;allocate(values(n));do i=1,n;values(i)=qlaplace(uniform_random(),mu,scale);end do
  end function rlaplace

  elemental real(dp) function dlaplace0(x,scale,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::scale;logical,intent(in),optional::log_density;value=dlaplace(x,0.0_dp,scale,log_density)
  end function dlaplace0
  elemental real(dp) function plaplace0(q,scale) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::scale;value=plaplace(q,0.0_dp,scale)
  end function plaplace0
  elemental real(dp) function qlaplace0(p,scale) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::scale;value=qlaplace(p,0.0_dp,scale)
  end function qlaplace0
  function rlaplace0(n,scale) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::scale;real(dp),allocatable::values(:);values=rlaplace(n,0.0_dp,scale)
  end function rlaplace0

  elemental real(dp) function dsdoubleweibull(x,mu,shape,scale,gamma,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,shape,scale,gamma;logical,intent(in),optional::log_density
    real(dp)::m,a,s,g,y,arg;logical::ld
    m=0.05_dp;if(present(mu))m=mu;a=1.0_dp;if(present(shape))a=shape;s=1.0_dp;if(present(scale))s=scale;g=1.0_dp;if(present(gamma))g=gamma;ld=.false.;if(present(log_density))ld=log_density
    if(s<=0.0_dp.or.a<=0.0_dp.or.g<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    y=(x-m)/s;if(y<0.0_dp)then;arg=g*abs(y);else;arg=y/g;end if
    if(arg<=tiny(1.0_dp))then
      if(a<1.0_dp)then;value=merge(huge(1.0_dp),huge(1.0_dp),ld);return
      else if(a>1.0_dp)then;value=merge(0.0_dp,-huge(1.0_dp),.not.ld);return
      else;value=log(a)-log(s)-log(g+1.0_dp/g);if(.not.ld)value=exp(value);return;end if
    end if
    value=(a-1.0_dp)*log(arg)-arg**a+log(a)-log(s)-log(g+1.0_dp/g);if(.not.ld)value=exp(value)
  end function dsdoubleweibull
  elemental real(dp) function psdoubleweibull(q,mu,shape,scale,gamma) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,shape,scale,gamma;real(dp)::m,a,s,g,y,arg
    m=0.05_dp;if(present(mu))m=mu;a=1.0_dp;if(present(shape))a=shape;s=1.0_dp;if(present(scale))s=scale;g=1.0_dp;if(present(gamma))g=gamma
    if(s<=0.0_dp.or.a<=0.0_dp.or.g<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    y=(q-m)/s;if(y<0.0_dp)then;arg=g*abs(y);value=exp(-arg**a)/(1.0_dp+g*g);else;arg=y/g;value=1.0_dp-g*g*exp(-arg**a)/(1.0_dp+g*g);end if
  end function psdoubleweibull
  elemental real(dp) function qsdoubleweibull(p,mu,shape,scale,gamma) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::mu,shape,scale,gamma;real(dp)::m,a,s,g,w
    m=0.05_dp;if(present(mu))m=mu;a=1.0_dp;if(present(shape))a=shape;s=1.0_dp;if(present(scale))s=scale;g=1.0_dp;if(present(gamma))g=gamma;w=1.0_dp+g*g
    if(p<=1.0_dp/w)then;value=m-s*(-log(w*p))**(1.0_dp/a)/g;else;value=m+s*(-log(w*(1.0_dp-p)/(w-1.0_dp)))**(1.0_dp/a)*g;end if
  end function qsdoubleweibull
  function rsdoubleweibull(n,mu,shape,scale,gamma) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,shape,scale,gamma;real(dp),allocatable::values(:);integer::i;allocate(values(n));do i=1,n;values(i)=qsdoubleweibull(uniform_random(),mu,shape,scale,gamma);end do
  end function rsdoubleweibull

  elemental real(dp) function ddoubleweibull(x,mu,shape,scale,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,shape,scale;logical,intent(in),optional::log_density;value=dsdoubleweibull(x,mu,shape,scale,1.0_dp,log_density)
  end function ddoubleweibull
  elemental real(dp) function pdoubleweibull(q,mu,shape,scale) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,shape,scale;value=psdoubleweibull(q,mu,shape,scale,1.0_dp)
  end function pdoubleweibull
  elemental real(dp) function qdoubleweibull(p,mu,shape,scale) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::mu,shape,scale;value=qsdoubleweibull(p,mu,shape,scale,1.0_dp)
  end function qdoubleweibull
  function rdoubleweibull(n,mu,shape,scale) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,shape,scale;real(dp),allocatable::values(:);values=rsdoubleweibull(n,mu,shape,scale,1.0_dp)
  end function rdoubleweibull

  elemental real(dp) function dslaplace(x,mu,scale,gamma,log_density) result(value)
    real(dp),intent(in)::x;real(dp),intent(in),optional::mu,scale,gamma;logical,intent(in),optional::log_density;value=dsdoubleweibull(x,mu,1.0_dp,scale,gamma,log_density)
  end function dslaplace
  elemental real(dp) function pslaplace(q,mu,scale,gamma) result(value)
    real(dp),intent(in)::q;real(dp),intent(in),optional::mu,scale,gamma;value=psdoubleweibull(q,mu,1.0_dp,scale,gamma)
  end function pslaplace
  elemental real(dp) function qslaplace(p,mu,scale,gamma) result(value)
    real(dp),intent(in)::p;real(dp),intent(in),optional::mu,scale,gamma;value=qsdoubleweibull(p,mu,1.0_dp,scale,gamma)
  end function qslaplace
  function rslaplace(n,mu,scale,gamma) result(values)
    integer,intent(in)::n;real(dp),intent(in),optional::mu,scale,gamma;real(dp),allocatable::values(:);values=rsdoubleweibull(n,mu,1.0_dp,scale,gamma)
  end function rslaplace

  elemental real(dp) function dst(x,df,mu,sigma,log_density) result(value)
    real(dp),intent(in)::x,df,mu,sigma;logical,intent(in),optional::log_density;logical::ld
    ld=.false.;if(present(log_density))ld=log_density
    if(sigma<0.0_dp.or.df<0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    value=log(student_pdf((x-mu)/sigma,df))-log(sigma);if(.not.ld)value=exp(value)
  end function dst
  elemental real(dp) function pst(q,df,mu,sigma) result(value)
    real(dp),intent(in)::q,df,mu,sigma;value=student_cdf((q-mu)/sigma,df)
  end function pst
  elemental real(dp) function qst(p,df,mu,sigma) result(value)
    real(dp),intent(in)::p,df,mu,sigma;value=student_quantile(p,df)*sigma+mu
  end function qst
  function rst(n,df,mu,sigma) result(values)
    integer,intent(in)::n;real(dp),intent(in)::df,mu,sigma;real(dp),allocatable::values(:);integer::i;allocate(values(n));do i=1,n;values(i)=mu+sigma*student_random(df);end do
  end function rst

  elemental real(dp) function dst0(x,df,sigma,log_density) result(value)
    real(dp),intent(in)::x,df,sigma;logical,intent(in),optional::log_density;value=dst(x,df,0.0_dp,sigma,log_density)
  end function dst0
  elemental real(dp) function pst0(q,df,sigma) result(value)
    real(dp),intent(in)::q,df,sigma;value=pst(q,df,0.0_dp,sigma)
  end function pst0
  elemental real(dp) function qst0(p,df,sigma) result(value)
    real(dp),intent(in)::p,df,sigma;value=qst(p,df,0.0_dp,sigma)
  end function qst0
  function rst0(n,df,sigma) result(values)
    integer,intent(in)::n;real(dp),intent(in)::df,sigma;real(dp),allocatable::values(:);values=rst(n,df,0.0_dp,sigma)
  end function rst0

  elemental real(dp) function psst(q,df,gamma,mu,sigma) result(value)
    real(dp),intent(in)::q,df,gamma,mu,sigma;real(dp)::x
    x=(q-mu)/sigma
    if(x<0.0_dp)then;value=2.0_dp/(gamma*gamma+1.0_dp)*student_cdf(gamma*x,df)
    else;value=1.0_dp/(gamma*gamma+1.0_dp)+2.0_dp/(1.0_dp+1.0_dp/(gamma*gamma))*(student_cdf(x/gamma,df)-0.5_dp);end if
  end function psst
  elemental real(dp) function qsst(p,df,gamma,mu,sigma) result(value)
    real(dp),intent(in)::p,df,gamma,mu,sigma;real(dp)::p0,z
    p0=1.0_dp/(gamma*gamma+1.0_dp)
    if(p<p0)then;z=student_quantile((gamma*gamma+1.0_dp)*p/2.0_dp,df)/gamma
    else;z=gamma*student_quantile((1.0_dp+1.0_dp/(gamma*gamma))*(p-p0)/2.0_dp+0.5_dp,df);end if
    value=mu+sigma*z
  end function qsst
  elemental real(dp) function dsst(x,df,gamma,mu,sigma,log_density) result(value)
    real(dp),intent(in)::x,df,gamma,mu,sigma;logical,intent(in),optional::log_density;logical::ld;real(dp)::z
    ld=.false.;if(present(log_density))ld=log_density
    if(sigma<0.0_dp.or.df<0.0_dp.or.gamma<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    z=(x-mu)/sigma;if(z<0.0_dp)then;value=log(student_pdf(gamma*z,df));else;value=log(student_pdf(z/gamma,df));end if
    value=value+log(2.0_dp/(gamma+1.0_dp/gamma))-log(sigma);if(.not.ld)value=exp(value)
  end function dsst
  function rsst(n,df,gamma,mu,sigma) result(values)
    integer,intent(in)::n;real(dp),intent(in)::df,gamma,mu,sigma;real(dp),allocatable::values(:);integer::i;allocate(values(n));do i=1,n;values(i)=qsst(uniform_random(),df,gamma,mu,sigma);end do
  end function rsst

  elemental real(dp) function dmarg(spec,x,log_density) result(value)
    type(margin_spec),intent(in)::spec;real(dp),intent(in)::x;logical,intent(in),optional::log_density
    logical::ld;ld=.false.;if(present(log_density))ld=log_density
    select case(spec%family)
    case(margin_uniform);if(x>=0.0_dp.and.x<=1.0_dp)then;value=merge(0.0_dp,1.0_dp,ld);else;value=merge(-huge(1.0_dp),0.0_dp,ld);end if
    case(margin_gauss);value=dgauss(x,spec%mu,spec%sigma,ld)
    case(margin_gauss0);value=dgauss0(x,spec%sigma,ld)
    case(margin_laplace);value=dlaplace(x,spec%mu,spec%sigma,ld)
    case(margin_laplace0);value=dlaplace0(x,spec%sigma,ld)
    case(margin_slaplace);value=dslaplace(x,spec%mu,spec%sigma,spec%gamma,ld)
    case(margin_doubleweibull);value=ddoubleweibull(x,spec%mu,spec%shape,spec%sigma,ld)
    case(margin_sdoubleweibull);value=dsdoubleweibull(x,spec%mu,spec%shape,spec%sigma,spec%gamma,ld)
    case(margin_student);value=dst(x,spec%df,spec%mu,spec%sigma,ld)
    case(margin_student0);value=dst0(x,spec%df,spec%sigma,ld)
    case(margin_sstudent);value=dsst(x,spec%df,spec%gamma,spec%mu,spec%sigma,ld)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function dmarg

  elemental real(dp) function pmarg(spec,q) result(value)
    type(margin_spec),intent(in)::spec;real(dp),intent(in)::q
    select case(spec%family)
    case(margin_uniform);value=min(max(q,0.0_dp),1.0_dp)
    case(margin_gauss);value=pgauss(q,spec%mu,spec%sigma)
    case(margin_gauss0);value=pgauss0(q,spec%sigma)
    case(margin_laplace);value=plaplace(q,spec%mu,spec%sigma)
    case(margin_laplace0);value=plaplace0(q,spec%sigma)
    case(margin_slaplace);value=pslaplace(q,spec%mu,spec%sigma,spec%gamma)
    case(margin_doubleweibull);value=pdoubleweibull(q,spec%mu,spec%shape,spec%sigma)
    case(margin_sdoubleweibull);value=psdoubleweibull(q,spec%mu,spec%shape,spec%sigma,spec%gamma)
    case(margin_student);value=pst(q,spec%df,spec%mu,spec%sigma)
    case(margin_student0);value=pst0(q,spec%df,spec%sigma)
    case(margin_sstudent);value=psst(q,spec%df,spec%gamma,spec%mu,spec%sigma)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function pmarg

  elemental real(dp) function qmarg(spec,p) result(value)
    type(margin_spec),intent(in)::spec;real(dp),intent(in)::p
    select case(spec%family)
    case(margin_uniform);value=p
    case(margin_gauss);value=qgauss(p,spec%mu,spec%sigma)
    case(margin_gauss0);value=qgauss0(p,spec%sigma)
    case(margin_laplace);value=qlaplace(p,spec%mu,spec%sigma)
    case(margin_laplace0);value=qlaplace0(p,spec%sigma)
    case(margin_slaplace);value=qslaplace(p,spec%mu,spec%sigma,spec%gamma)
    case(margin_doubleweibull);value=qdoubleweibull(p,spec%mu,spec%shape,spec%sigma)
    case(margin_sdoubleweibull);value=qsdoubleweibull(p,spec%mu,spec%shape,spec%sigma,spec%gamma)
    case(margin_student);value=qst(p,spec%df,spec%mu,spec%sigma)
    case(margin_student0);value=qst0(p,spec%df,spec%sigma)
    case(margin_sstudent);value=qsst(p,spec%df,spec%gamma,spec%mu,spec%sigma)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function qmarg

  function rmarg(spec,n) result(values)
    type(margin_spec),intent(in)::spec;integer,intent(in)::n;real(dp),allocatable::values(:);integer::i
    allocate(values(n));do i=1,n;values(i)=qmarg(spec,uniform_random());end do
  end function rmarg

  subroutine fit_margin(spec,data,result,error,compute_hessian,max_iter)
    type(margin_spec),intent(in)::spec;real(dp),intent(in)::data(:)
    type(margin_fit_result),intent(out)::result;type(tsc_error),intent(out)::error
    logical,intent(in),optional::compute_hessian;integer,intent(in),optional::max_iter
    type(margin_context)::context;type(optimizer_result)::opt
    real(dp),allocatable::start(:),lower(:),upper(:);logical::do_hess
    call clear_error(error);do_hess=.false.;if(present(compute_hessian))do_hess=compute_hessian
    if(size(data)<2.or.any(.not.ieee_is_finite(data)))then;call set_error(error,tsc_invalid_input,'margin data must contain at least two finite values');return;end if
    context%template=spec;context%data=data
    call parameter_vectors(spec,data,start,lower,upper,error);if(.not.error%ok())return
    call minimize_nelder_mead(margin_objective,start,lower,upper,context,opt,max_iter=max_iter)
    result%margin=spec;call assign_parameters(result%margin,opt%par);result%parameters=opt%par
    result%log_likelihood=-opt%value;result%convergence=opt%convergence;result%iterations=opt%iterations
    if(do_hess)then;call finite_hessian(margin_objective,opt%par,context,result%hessian);call safe_standard_errors(result%hessian,result%standard_errors);end if
  end subroutine fit_margin

  subroutine parameter_vectors(spec,data,start,lower,upper,error)
    type(margin_spec),intent(in)::spec;real(dp),intent(in)::data(:);real(dp),allocatable,intent(out)::start(:),lower(:),upper(:);type(tsc_error),intent(out)::error
    real(dp)::meanx,sdx;call clear_error(error);meanx=sum(data)/real(size(data),dp);sdx=sqrt(sum((data-meanx)**2)/real(max(1,size(data)-1),dp));sdx=max(sdx,1.0e-4_dp)
    select case(spec%family)
    case(margin_gauss,margin_laplace);start=[meanx,sdx];lower=[minval(data)-10.0_dp*sdx,1.0e-8_dp];upper=[maxval(data)+10.0_dp*sdx,100.0_dp*sdx]
    case(margin_gauss0,margin_laplace0);start=[sdx];lower=[1.0e-8_dp];upper=[100.0_dp*sdx]
    case(margin_slaplace);start=[meanx,sdx,max(spec%gamma,0.2_dp)];lower=[minval(data)-10.0_dp*sdx,1.0e-8_dp,0.05_dp];upper=[maxval(data)+10.0_dp*sdx,100.0_dp*sdx,20.0_dp]
    case(margin_doubleweibull);start=[meanx,max(spec%shape,1.0_dp),sdx];lower=[minval(data)-10.0_dp*sdx,0.05_dp,1.0e-8_dp];upper=[maxval(data)+10.0_dp*sdx,20.0_dp,100.0_dp*sdx]
    case(margin_sdoubleweibull);start=[meanx,max(spec%shape,1.0_dp),sdx,max(spec%gamma,0.2_dp)];lower=[minval(data)-10.0_dp*sdx,0.05_dp,1.0e-8_dp,0.05_dp];upper=[maxval(data)+10.0_dp*sdx,20.0_dp,100.0_dp*sdx,20.0_dp]
    case(margin_student);start=[max(spec%df,4.0_dp),meanx,sdx];lower=[1.01_dp,minval(data)-10.0_dp*sdx,1.0e-8_dp];upper=[200.0_dp,maxval(data)+10.0_dp*sdx,100.0_dp*sdx]
    case(margin_student0);start=[max(spec%df,4.0_dp),sdx];lower=[1.01_dp,1.0e-8_dp];upper=[200.0_dp,100.0_dp*sdx]
    case(margin_sstudent);start=[max(spec%df,4.0_dp),max(spec%gamma,0.2_dp),meanx,sdx];lower=[1.01_dp,0.05_dp,minval(data)-10.0_dp*sdx,1.0e-8_dp];upper=[200.0_dp,20.0_dp,maxval(data)+10.0_dp*sdx,100.0_dp*sdx]
    case default;allocate(start(0),lower(0),upper(0));call set_error(error,tsc_invalid_input,'margin family is not parametric')
    end select
  end subroutine parameter_vectors

  real(dp) function margin_objective(theta,context_any) result(value)
    real(dp),intent(in)::theta(:);class(*),intent(inout)::context_any
    type(margin_spec)::spec;integer::i
    select type(context=>context_any);type is(margin_context);spec=context%template;call assign_parameters(spec,theta);value=0.0_dp;do i=1,size(context%data);value=value-dmarg(spec,context%data(i),.true.);if(.not.ieee_is_finite(value))then;value=huge(1.0_dp)/100.0_dp;return;end if;end do;class default;value=huge(1.0_dp);end select
  end function margin_objective
end module tscopula_margins
