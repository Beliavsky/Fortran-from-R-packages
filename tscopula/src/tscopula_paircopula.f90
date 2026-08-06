! SPDX-License-Identifier: GPL-3.0-only
module tscopula_paircopula
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tscopula_kinds, only : dp, pi
  use tscopula_math, only : normal_cdf, normal_quantile, student_cdf, &
    student_quantile, clamp_probability, uniform_random, integrate_simpson
  implicit none
  private

  integer, parameter, public :: cop_indep = 0
  integer, parameter, public :: cop_gauss = 1
  integer, parameter, public :: cop_student = 2
  integer, parameter, public :: cop_clayton = 3
  integer, parameter, public :: cop_gumbel = 4
  integer, parameter, public :: cop_frank = 5
  integer, parameter, public :: cop_joe = 6
  integer, parameter, public :: cop_bb1 = 7

  type, public :: pair_copula
    integer :: family = cop_indep
    integer :: rotation = 0
    real(dp) :: par1 = 0.0_dp
    real(dp) :: par2 = 0.0_dp
  end type pair_copula

  public :: paircop, pair_cdf, pair_density, pair_log_density
  public :: pair_h1, pair_h2, pair_hinv1, pair_hinv2, pair_simulate
  public :: pair_kendall, kendall_to_parameter, family_from_name, family_name

  type :: cdf_context
    type(pair_copula) :: copula
    real(dp) :: fixed = 0.5_dp
  end type cdf_context

  type :: frank_context
    real(dp) :: theta = 1.0_dp
  end type frank_context

  type :: joe_context
    real(dp) :: theta = 1.0_dp
  end type joe_context
contains
  function paircop(family,par1,par2,rotation) result(cop)
    character(len=*),intent(in)::family
    real(dp),intent(in),optional::par1,par2
    integer,intent(in),optional::rotation
    type(pair_copula)::cop
    cop%family=family_from_name(family)
    if(present(par1))cop%par1=par1
    if(present(par2))cop%par2=par2
    if(present(rotation))cop%rotation=modulo(rotation,360)
  end function paircop

  integer function family_from_name(name) result(family)
    character(len=*),intent(in)::name
    character(len=:),allocatable::s
    s=lower_string(trim(name))
    select case(s)
    case('indep','independence');family=cop_indep
    case('gauss','gaussian','normal');family=cop_gauss
    case('t','student','student-t');family=cop_student
    case('clayton');family=cop_clayton
    case('gumbel');family=cop_gumbel
    case('frank');family=cop_frank
    case('joe');family=cop_joe
    case('bb1');family=cop_bb1
    case default;family=-1
    end select
  end function family_from_name

  function family_name(family) result(name)
    integer,intent(in)::family
    character(len=:),allocatable::name
    select case(family)
    case(cop_indep);name='indep'
    case(cop_gauss);name='gauss'
    case(cop_student);name='t'
    case(cop_clayton);name='clayton'
    case(cop_gumbel);name='gumbel'
    case(cop_frank);name='frank'
    case(cop_joe);name='joe'
    case(cop_bb1);name='bb1'
    case default;name='unknown'
    end select
  end function family_name

  pure function lower_string(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,c
    do i=1,len(s);c=iachar(s(i:i));if(c>=65.and.c<=90)then;out(i:i)=achar(c+32);else;out(i:i)=s(i:i);end if;end do
  end function lower_string

  real(dp) function pair_cdf(cop,u,v) result(value)
    type(pair_copula),intent(in)::cop
    real(dp),intent(in)::u,v
    real(dp)::uu,vv,base
    uu=min(max(u,0.0_dp),1.0_dp);vv=min(max(v,0.0_dp),1.0_dp)
    select case(modulo(cop%rotation,360))
    case(0);value=base_cdf(cop,uu,vv)
    case(90);base=base_cdf(cop,1.0_dp-uu,vv);value=vv-base
    case(180);base=base_cdf(cop,1.0_dp-uu,1.0_dp-vv);value=uu+vv-1.0_dp+base
    case(270);base=base_cdf(cop,uu,1.0_dp-vv);value=uu-base
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    value=min(max(value,0.0_dp),min(uu,vv))
  end function pair_cdf

  real(dp) function base_cdf(cop,u,v) result(value)
    type(pair_copula),intent(in)::cop
    real(dp),intent(in)::u,v
    real(dp)::a,b,t,theta,delta,d,eu,ev,e1
    type(cdf_context)::context
    if(u<=0.0_dp.or.v<=0.0_dp)then;value=0.0_dp;return
    else if(u>=1.0_dp)then;value=v;return
    else if(v>=1.0_dp)then;value=u;return
    end if
    select case(cop%family)
    case(cop_indep);value=u*v
    case(cop_gauss,cop_student)
      context%copula=cop;context%fixed=v
      value=integrate_simpson(base_h1_integrand,0.0_dp,u,context,2.0e-8_dp,16)
    case(cop_clayton)
      theta=cop%par1;if(theta<=0.0_dp)then;value=u*v;else;a=u**(-theta)+v**(-theta)-1.0_dp;value=max(a,0.0_dp)**(-1.0_dp/theta);end if
    case(cop_gumbel)
      theta=cop%par1;if(theta<1.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;a=(-log(u))**theta+(-log(v))**theta;value=exp(-a**(1.0_dp/theta));end if
    case(cop_frank)
      theta=cop%par1
      if(abs(theta)<=1.0e-10_dp)then;value=u*v
      else;eu=exp(-theta*u)-1.0_dp;ev=exp(-theta*v)-1.0_dp;e1=exp(-theta)-1.0_dp;value=-log(1.0_dp+eu*ev/e1)/theta;end if
    case(cop_joe)
      theta=cop%par1;if(theta<1.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;a=(1.0_dp-u)**theta;b=(1.0_dp-v)**theta;t=a+b-a*b;value=1.0_dp-t**(1.0_dp/theta);end if
    case(cop_bb1)
      theta=cop%par1;delta=cop%par2
      if(theta<=0.0_dp.or.delta<1.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan)
      else;a=(u**(-theta)-1.0_dp)**delta+(v**(-theta)-1.0_dp)**delta;d=1.0_dp+a**(1.0_dp/delta);value=d**(-1.0_dp/theta);end if
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function base_cdf

  real(dp) function base_h1_integrand(x,context_any) result(value)
    real(dp),intent(in)::x;class(*),intent(inout)::context_any
    select type(context=>context_any)
    type is(cdf_context);value=base_h1(context%copula,clamp_probability(x),context%fixed)
    class default;value=0.0_dp
    end select
  end function base_h1_integrand

  real(dp) function pair_h2(cop,u,v) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u,v
    select case(modulo(cop%rotation,360))
    case(0);value=base_h2(cop,u,v)
    case(90);value=1.0_dp-base_h2(cop,1.0_dp-u,v)
    case(180);value=1.0_dp-base_h2(cop,1.0_dp-u,1.0_dp-v)
    case(270);value=base_h2(cop,u,1.0_dp-v)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    value=min(max(value,0.0_dp),1.0_dp)
  end function pair_h2

  real(dp) function pair_h1(cop,u,v) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u,v
    select case(modulo(cop%rotation,360))
    case(0);value=base_h1(cop,u,v)
    case(90);value=base_h1(cop,1.0_dp-u,v)
    case(180);value=1.0_dp-base_h1(cop,1.0_dp-u,1.0_dp-v)
    case(270);value=1.0_dp-base_h1(cop,u,1.0_dp-v)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    value=min(max(value,0.0_dp),1.0_dp)
  end function pair_h1

  real(dp) function base_h2(cop,u_in,v_in) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u_in,v_in
    real(dp)::u,v,rho,nu,zu,zv,s,a,b,t,theta,delta,eu,ev,e1,big
    u=clamp_probability(u_in);v=clamp_probability(v_in)
    select case(cop%family)
    case(cop_indep);value=u
    case(cop_gauss)
      rho=cop%par1;s=sqrt(max(1.0_dp-rho*rho,tiny(1.0_dp)));value=normal_cdf((normal_quantile(u)-rho*normal_quantile(v))/s)
    case(cop_student)
      rho=cop%par1;nu=max(cop%par2,2.01_dp);zu=student_quantile(u,nu);zv=student_quantile(v,nu)
      s=sqrt((nu+1.0_dp)/((nu+zv*zv)*max(1.0_dp-rho*rho,tiny(1.0_dp))))
      value=student_cdf((zu-rho*zv)*s,nu+1.0_dp)
    case(cop_clayton)
      theta=cop%par1;if(theta<=0.0_dp)then;value=u;else;a=u**(-theta)+v**(-theta)-1.0_dp;value=v**(-theta-1.0_dp)*a**(-1.0_dp/theta-1.0_dp);end if
    case(cop_gumbel)
      theta=cop%par1;a=(-log(u))**theta+(-log(v))**theta;t=exp(-a**(1.0_dp/theta));value=t*a**(1.0_dp/theta-1.0_dp)*(-log(v))**(theta-1.0_dp)/v
    case(cop_frank)
      theta=cop%par1;if(abs(theta)<=1.0e-10_dp)then;value=u;else;eu=exp(-theta*u)-1.0_dp;ev=exp(-theta*v)-1.0_dp;e1=exp(-theta)-1.0_dp;value=eu*exp(-theta*v)/(e1+eu*ev);end if
    case(cop_joe)
      theta=cop%par1;a=(1.0_dp-u)**theta;b=(1.0_dp-v)**theta;t=a+b-a*b;value=t**(1.0_dp/theta-1.0_dp)*(1.0_dp-a)*(1.0_dp-v)**(theta-1.0_dp)
    case(cop_bb1)
      theta=cop%par1;delta=cop%par2;a=(u**(-theta)-1.0_dp)**delta+(v**(-theta)-1.0_dp)**delta;big=1.0_dp+a**(1.0_dp/delta)
      value=big**(-1.0_dp/theta-1.0_dp)*a**(1.0_dp/delta-1.0_dp)*(v**(-theta)-1.0_dp)**(delta-1.0_dp)*v**(-theta-1.0_dp)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function base_h2

  real(dp) function base_h1(cop,u_in,v_in) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u_in,v_in
    type(pair_copula)::swapped
    swapped=cop
    value=base_h2(swapped,v_in,u_in)
  end function base_h1

  real(dp) function pair_density(cop,u_in,v_in) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u_in,v_in
    real(dp)::u,v,uu,vv,rho,nu,x,y,q,den,h,lo,hi
    u=clamp_probability(u_in);v=clamp_probability(v_in)
    select case(modulo(cop%rotation,360))
    case(0);uu=u;vv=v
    case(90);uu=1.0_dp-u;vv=v
    case(180);uu=1.0_dp-u;vv=1.0_dp-v
    case(270);uu=u;vv=1.0_dp-v
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan);return
    end select
    select case(cop%family)
    case(cop_indep);value=1.0_dp
    case(cop_gauss)
      rho=cop%par1;x=normal_quantile(uu);y=normal_quantile(vv);den=max(1.0_dp-rho*rho,tiny(1.0_dp))
      value=exp(-(x*x-2.0_dp*rho*x*y+y*y)/(2.0_dp*den)+0.5_dp*(x*x+y*y))/sqrt(den)
    case(cop_student)
      rho=cop%par1;nu=max(cop%par2,2.01_dp);x=student_quantile(uu,nu);y=student_quantile(vv,nu);den=max(1.0_dp-rho*rho,tiny(1.0_dp));q=(x*x-2.0_dp*rho*x*y+y*y)/den
      value=exp(log_gamma(0.5_dp*(nu+2.0_dp))+log_gamma(0.5_dp*nu) &
        -2.0_dp*log_gamma(0.5_dp*(nu+1.0_dp))-0.5_dp*log(den) &
        -0.5_dp*(nu+2.0_dp)*log(1.0_dp+q/nu) &
        +0.5_dp*(nu+1.0_dp)*(log(1.0_dp+x*x/nu)+log(1.0_dp+y*y/nu)))
    case default
      h=1.0e-5_dp*max(0.1_dp,min(uu,1.0_dp-uu));lo=max(uu-h,1.0e-10_dp);hi=min(uu+h,1.0_dp-1.0e-10_dp)
      value=(base_h2(cop,hi,vv)-base_h2(cop,lo,vv))/(hi-lo)
      value=max(value,tiny(1.0_dp))
    end select
  end function pair_density

  real(dp) function pair_log_density(cop,u,v) result(value)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::u,v
    value=log(max(pair_density(cop,u,v),tiny(1.0_dp)))
  end function pair_log_density

  real(dp) function pair_hinv2(cop,p,v) result(u)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::p,v
    real(dp)::lo,hi,mid;integer::iter
    lo=0.0_dp;hi=1.0_dp
    do iter=1,80;mid=0.5_dp*(lo+hi);if(pair_h2(cop,mid,v)<p)then;lo=mid;else;hi=mid;end if;end do
    u=0.5_dp*(lo+hi)
  end function pair_hinv2

  real(dp) function pair_hinv1(cop,p,u) result(v)
    type(pair_copula),intent(in)::cop;real(dp),intent(in)::p,u
    real(dp)::lo,hi,mid;integer::iter
    lo=0.0_dp;hi=1.0_dp
    do iter=1,80;mid=0.5_dp*(lo+hi);if(pair_h1(cop,u,mid)<p)then;lo=mid;else;hi=mid;end if;end do
    v=0.5_dp*(lo+hi)
  end function pair_hinv1

  subroutine pair_simulate(cop,n,values)
    type(pair_copula),intent(in)::cop;integer,intent(in)::n;real(dp),allocatable,intent(out)::values(:,:)
    integer::i;allocate(values(n,2));do i=1,n;values(i,2)=uniform_random();values(i,1)=pair_hinv2(cop,uniform_random(),values(i,2));end do
  end subroutine pair_simulate

  real(dp) function pair_kendall(cop) result(tau)
    type(pair_copula),intent(in)::cop
    type(frank_context)::fc;type(joe_context)::jc;real(dp)::base
    select case(cop%family)
    case(cop_indep);base=0.0_dp
    case(cop_gauss,cop_student);base=2.0_dp*asin(cop%par1)/pi
    case(cop_clayton);base=cop%par1/(cop%par1+2.0_dp)
    case(cop_gumbel);base=1.0_dp-1.0_dp/cop%par1
    case(cop_frank);fc%theta=cop%par1;base=frank_tau(cop%par1)
    case(cop_joe);jc%theta=cop%par1;base=joe_tau(cop%par1)
    case(cop_bb1);base=1.0_dp-2.0_dp/(cop%par2*(cop%par1+2.0_dp))
    case default;base=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    if(modulo(cop%rotation,360)==90.or.modulo(cop%rotation,360)==270)then;tau=-base;else;tau=base;end if
  end function pair_kendall

  real(dp) function kendall_to_parameter(family,tau,aux) result(par)
    integer,intent(in)::family;real(dp),intent(in)::tau;real(dp),intent(in),optional::aux
    real(dp)::lo,hi,mid,t;integer::iter
    select case(family)
    case(cop_indep);par=0.0_dp
    case(cop_gauss,cop_student);par=sin(0.5_dp*pi*tau)
    case(cop_clayton,cop_bb1);if(tau>=1.0_dp)then;par=huge(1.0_dp);else;par=2.0_dp*tau/(1.0_dp-tau);end if
    case(cop_gumbel);par=1.0_dp/(1.0_dp-tau)
    case(cop_frank)
      if(abs(tau)<=1.0e-12_dp)then;par=0.0_dp;return;end if
      if(tau>0.0_dp)then;lo=1.0e-6_dp;hi=60.0_dp;else;lo=-60.0_dp;hi=-1.0e-6_dp;end if
      do iter=1,80;mid=0.5_dp*(lo+hi);t=frank_tau(mid);if(t<tau)then;lo=mid;else;hi=mid;end if;end do;par=0.5_dp*(lo+hi)
    case(cop_joe)
      lo=1.0_dp;hi=80.0_dp;do iter=1,80;mid=0.5_dp*(lo+hi);t=joe_tau(mid);if(t<abs(tau))then;lo=mid;else;hi=mid;end if;end do;par=0.5_dp*(lo+hi)
    case default;par=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    if(present(aux)) par=par+0.0_dp*aux
  end function kendall_to_parameter

  real(dp) function frank_tau(theta) result(tau)
    real(dp),intent(in)::theta;type(frank_context)::context;real(dp)::integral
    if(abs(theta)<=1.0e-8_dp)then;tau=theta/9.0_dp;return;end if
    context%theta=theta
    if(theta>0.0_dp)then;integral=integrate_simpson(frank_integrand,0.0_dp,theta,context,1.0e-9_dp)
    else;integral=-integrate_simpson(frank_integrand,theta,0.0_dp,context,1.0e-9_dp);end if
    tau=1.0_dp-4.0_dp/theta+4.0_dp*integral/(theta*theta)
  end function frank_tau

  real(dp) function frank_integrand(x,context_any) result(value)
    real(dp),intent(in)::x;class(*),intent(inout)::context_any
    if(abs(x)<1.0e-7_dp)then;value=1.0_dp-x/2.0_dp+x*x/12.0_dp;else;value=x/(exp(x)-1.0_dp);end if
    select type(context_any);class default;continue;end select
  end function frank_integrand

  real(dp) function joe_tau(theta) result(tau)
    real(dp),intent(in)::theta;type(joe_context)::context;real(dp)::integral
    if(theta<=1.0_dp)then;tau=0.0_dp;return;end if
    context%theta=theta;integral=integrate_simpson(joe_tau_integrand,1.0e-10_dp,1.0_dp-1.0e-10_dp,context,1.0e-8_dp)
    tau=1.0_dp+4.0_dp*integral
  end function joe_tau

  real(dp) function joe_tau_integrand(u,context_any) result(value)
    real(dp),intent(in)::u;class(*),intent(inout)::context_any
    real(dp)::theta,phi,phip,a
    select type(context=>context_any);type is(joe_context);theta=context%theta;class default;theta=1.0_dp;end select
    a=(1.0_dp-u)**theta;phi=-log(1.0_dp-a);phip=-theta*(1.0_dp-u)**(theta-1.0_dp)/(1.0_dp-a)
    value=phi/phip
  end function joe_tau_integrand
end module tscopula_paircopula
