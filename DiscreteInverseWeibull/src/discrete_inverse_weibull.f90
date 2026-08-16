! SPDX-License-Identifier: GPL-2.0-only
module discrete_inverse_weibull
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, &
      ieee_is_nan, ieee_is_finite
   use rsolnp, only : solnp_problem, solnp_control, solnp_result, solnp, solnp_success
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter :: log2 = log(2.0_dp)

   type, public :: diw_moments
      real(dp) :: ex = 0.0_dp
      real(dp) :: ex2 = 0.0_dp
      integer :: xmax = 0
   end type

   type, public :: diw_control
      real(dp) :: eps = 1.0e-4_dp
      integer :: nmax = 1000
      real(dp) :: beta1 = 1.0_dp
      real(dp) :: z = 0.1_dp
      real(dp) :: r = 0.1_dp
      real(dp) :: leps = 1.0e-4_dp
      integer :: heuristic_max_iter = 10000
      integer :: solnp_max_iter = 400
      real(dp) :: solnp_tol = 1.0e-8_dp
      integer :: trace = 0
   end type

   type, public :: diw_estimate
      real(dp) :: q = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: status = 0
      integer :: iterations = 0
      character(len=160) :: message = 'success'
   end type

   type :: moment_context
      real(dp), allocatable :: x(:)
      real(dp) :: eps = 1.0e-4_dp
      integer :: nmax = 1000
   end type

   public :: ddiweibull, pdiweibull, qdiweibull, rdiweibull
   public :: hrdiweibull, ahrdiweibull, ediweibull
   public :: loglikediw, lossdiw, heuristic, estdiweibull
   public :: set_rng_seed

contains

   pure logical function valid_pars(q,beta)
      real(dp), intent(in) :: q,beta
      valid_pars = q > 0.0_dp .and. q < 1.0_dp .and. beta > 0.0_dp
   end function

   pure logical function valid_x(x)
      real(dp), intent(in) :: x
      valid_x = x >= 1.0_dp .and. abs(x-anint(x)) <= 16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))
   end function

   pure real(dp) function one_minus_exp(z)
      real(dp), intent(in) :: z
      if (abs(z) < 1.0e-5_dp) then
         one_minus_exp = -z*(1.0_dp + z*(0.5_dp + z*(1.0_dp/6.0_dp + z/24.0_dp)))
      else
         one_minus_exp = 1.0_dp-exp(z)
      end if
   end function

   pure real(dp) function log1mexp(z)
      real(dp), intent(in) :: z
      if (z >= 0.0_dp) then
         log1mexp = ieee_value(0.0_dp,ieee_quiet_nan)
      else if (z > -log2) then
         log1mexp = log(-2.0_dp*exp(0.5_dp*z)*sinh(0.5_dp*z))
      else
         log1mexp = log(1.0_dp-exp(z))
      end if
   end function

   pure real(dp) function log_pmf_int(k,q,beta)
      integer, intent(in) :: k
      real(dp), intent(in) :: q,beta
      real(dp) :: lq,a,b
      if (.not. valid_pars(q,beta) .or. k < 1) then
         log_pmf_int = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      lq=log(q)
      if (k == 1) then
         log_pmf_int=lq
      else
         a=lq*real(k,dp)**(-beta)
         b=lq*real(k-1,dp)**(-beta)
         log_pmf_int=a+log1mexp(b-a)
      end if
   end function

   real(dp) function ddiweibull(x,q,beta,log_p)
      real(dp), intent(in) :: x,q,beta
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp=.false.; if(present(log_p)) lp=log_p
      if (.not. valid_x(x) .or. .not. valid_pars(q,beta)) then
         ddiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      lv=log_pmf_int(nint(x),q,beta)
      if(lp) then
         ddiweibull=lv
      else
         ddiweibull=exp(lv)
      end if
   end function

   real(dp) function pdiweibull(x,q,beta,lower_tail,log_p)
      real(dp), intent(in) :: x,q,beta
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lower,lp
      real(dp) :: p
      lower=.true.; lp=.false.
      if(present(lower_tail)) lower=lower_tail
      if(present(log_p)) lp=log_p
      if (.not. valid_pars(q,beta)) then
         pdiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      if(x<1.0_dp) then
         p=0.0_dp
      else
         p=exp(log(q)*real(floor(x),dp)**(-beta))
      end if
      if(.not.lower) p=1.0_dp-p
      if(lp) then
         if(p>0.0_dp) then; pdiweibull=log(p); else; pdiweibull=-ieee_value(0.0_dp,ieee_positive_inf); end if
      else
         pdiweibull=p
      end if
   end function

   real(dp) function qdiweibull(p,q,beta,lower_tail,log_p)
      real(dp), intent(in) :: p,q,beta
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lower,lp
      real(dp) :: pp,v
      lower=.true.; lp=.false.
      if(present(lower_tail)) lower=lower_tail
      if(present(log_p)) lp=log_p
      if(.not.valid_pars(q,beta)) then
         qdiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      pp=p; if(lp) pp=exp(p); if(.not.lower) pp=1.0_dp-pp
      if(pp<=0.0_dp .or. pp>=1.0_dp) then
         qdiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      v=(log(pp)/log(q))**(-1.0_dp/beta)
      qdiweibull=ceiling(v)
   end function

   subroutine rdiweibull(x,q,beta)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: q,beta
      integer :: i
      real(dp) :: u
      if(.not.valid_pars(q,beta)) error stop 'rdiweibull: require 0<q<1 and beta>0'
      do i=1,size(x)
         call random_number(u)
         do while(u<=0.0_dp .or. u>=1.0_dp); call random_number(u); end do
         x(i)=ceiling((log(u)/log(q))**(-1.0_dp/beta))
      end do
   end subroutine

   real(dp) function hrdiweibull(x,q,beta)
      real(dp), intent(in) :: x,q,beta
      integer :: k
      real(dp) :: lp,s
      if(.not.valid_x(x) .or. .not.valid_pars(q,beta)) then
         hrdiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      k=nint(x); lp=log_pmf_int(k,q,beta)
      if(k==1) then
         s=1.0_dp
      else
         s=one_minus_exp(log(q)*real(k-1,dp)**(-beta))
      end if
      hrdiweibull=exp(lp)/s
   end function

   real(dp) function ahrdiweibull(x,q,beta)
      real(dp), intent(in) :: x,q,beta
      integer :: k
      real(dp) :: sprev,snow
      if(.not.valid_x(x) .or. .not.valid_pars(q,beta)) then
         ahrdiweibull=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      k=nint(x)
      if(k==1) then
         sprev=1.0_dp
      else
         sprev=one_minus_exp(log(q)*real(k-1,dp)**(-beta))
      end if
      snow=one_minus_exp(log(q)*real(k,dp)**(-beta))
      ahrdiweibull=log(sprev/snow)
   end function

   function ediweibull(q,beta,eps,nmax) result(ans)
      real(dp), intent(in) :: q,beta
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      type(diw_moments) :: ans
      real(dp) :: ee,qq,surv
      integer :: nm,x,xq
      ee=1.0e-4_dp; if(present(eps)) ee=eps
      nm=1000; if(present(nmax)) nm=nmax
      if(.not.valid_pars(q,beta) .or. ee<=0.0_dp .or. ee>=1.0_dp) then
         ans%ex=ieee_value(0.0_dp,ieee_quiet_nan); ans%ex2=ans%ex; return
      end if
      qq=qdiweibull(1.0_dp-ee,q,beta)
      if(ieee_is_nan(qq) .or. qq>real(huge(1),dp)/4.0_dp) then
         xq=nm
      else
         xq=int(qq)
      end if
      ans%xmax=max(2*xq,nm)
      if(beta>1.0_dp) then
         ans%ex=1.0_dp
         do x=1,ans%xmax
            surv=one_minus_exp(log(q)*real(x,dp)**(-beta))
            ans%ex=ans%ex+surv
         end do
      else
         ans%ex=ieee_value(0.0_dp,ieee_positive_inf)
      end if
      if(beta>2.0_dp) then
         ans%ex2=0.0_dp
         do x=1,ans%xmax
            surv=one_minus_exp(log(q)*real(x,dp)**(-beta))
            ans%ex2=ans%ex2+2.0_dp*real(x,dp)*surv
         end do
         ans%ex2=ans%ex2+ans%ex
      else
         ans%ex2=ieee_value(0.0_dp,ieee_positive_inf)
      end if
   end function

   real(dp) function loglikediw(x,q,beta)
      real(dp), intent(in) :: x(:),q,beta
      integer :: i
      if(.not.valid_pars(q,beta) .or. any(x<1.0_dp) .or. any(abs(x-anint(x))>1.0e-10_dp)) then
         loglikediw=huge(1.0_dp); return
      end if
      loglikediw=0.0_dp
      do i=1,size(x)
         loglikediw=loglikediw-log_pmf_int(nint(x(i)),q,beta)
      end do
   end function

   real(dp) function lossdiw(x,par,eps,nmax)
      real(dp), intent(in) :: x(:),par(2)
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      type(diw_moments) :: e
      real(dp) :: m1,m2
      if(par(1)<=0.0_dp .or. par(1)>=1.0_dp .or. par(2)<=2.0_dp) then
         lossdiw=huge(1.0_dp)/1000.0_dp; return
      end if
      e=ediweibull(par(1),par(2),eps,nmax)
      if(.not.ieee_is_finite(e%ex) .or. .not.ieee_is_finite(e%ex2) .or. &
         e%ex>=huge(1.0_dp)/1000.0_dp .or. e%ex2>=huge(1.0_dp)/1000.0_dp) then
         lossdiw=huge(1.0_dp)/1000.0_dp; return
      end if
      m1=sum(x)/real(size(x),dp); m2=sum(x*x)/real(size(x),dp)
      lossdiw=(m1-e%ex)**2+(m2-e%ex2)**2
   end function

   real(dp) function q_objective(x,beta,q)
      real(dp), intent(in) :: x(:),beta,q
      q_objective=loglikediw(x,q,beta)
   end function

   subroutine optimize_q(x,beta,qbest,fbest)
      real(dp), intent(in) :: x(:),beta
      real(dp), intent(out) :: qbest,fbest
      real(dp) :: a,b,c,d,fc,fd,phi
      integer :: it
      a=1.0e-12_dp; b=1.0_dp-1.0e-12_dp; phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      c=b-phi*(b-a); d=a+phi*(b-a)
      fc=q_objective(x,beta,c); fd=q_objective(x,beta,d)
      do it=1,300
         if(abs(b-a)<1.0e-11_dp) exit
         if(fc<fd) then
            b=d; d=c; fd=fc; c=b-phi*(b-a); fc=q_objective(x,beta,c)
         else
            a=c; c=d; fc=fd; d=a+phi*(b-a); fd=q_objective(x,beta,d)
         end if
      end do
      qbest=0.5_dp*(a+b); fbest=q_objective(x,beta,qbest)
   end subroutine

   function heuristic(x,beta1,z,r,leps,max_iter) result(est)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: beta1,z,r,leps
      integer, intent(in), optional :: max_iter
      type(diw_estimate) :: est
      real(dp) :: b0,zz,rr,tol,beta0,beta_star,lold,lnew
      real(dp) :: vb(3),vq(3),lv(3)
      integer :: i,it,mx,imin
      b0=1.0_dp; zz=0.1_dp; rr=0.1_dp; tol=0.01_dp; mx=10000
      if(present(beta1)) b0=beta1; if(present(z)) zz=z; if(present(r)) rr=r
      if(present(leps)) tol=leps; if(present(max_iter)) mx=max_iter
      if(b0<=0.0_dp .or. zz<=0.0_dp .or. rr<=0.0_dp .or. rr>=1.0_dp) then
         est%status=2; est%message='invalid heuristic controls'; return
      end if
      beta0=b0; beta_star=b0+1.0_dp; lold=huge(1.0_dp); lnew=0.0_dp
      do it=1,mx
         vb=[max(beta0-zz,1.0e-8_dp),beta0,beta0+zz]
         do i=1,3; call optimize_q(x,vb(i),vq(i),lv(i)); end do
         imin=minloc(lv,dim=1); beta_star=vb(imin)
         if(imin==2) then
            zz=zz*rr; lnew=minval(lv)
         end if
         if(imin==2 .and. abs(lnew-lold)<=tol) exit
         lold=lnew; beta0=beta_star
      end do
      call optimize_q(x,beta_star,est%q,est%objective)
      est%beta=beta_star; est%iterations=it
      if(it>mx) then; est%status=1; est%message='heuristic iteration limit reached'; end if
   end function

   subroutine mom_objective(par,value,data)
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if(.not.present(data)) then; value=huge(1.0_dp); return; end if
      select type(ctx=>data)
      type is(moment_context)
         value=lossdiw(ctx%x,par,ctx%eps,ctx%nmax)
      class default
         value=huge(1.0_dp)
      end select
   end subroutine

   function estimate_pp(x) result(est)
      real(dp), intent(in) :: x(:)
      type(diw_estimate) :: est
      real(dp), allocatable :: lx(:),y(:)
      real(dp) :: xb,yb,sxx,sxy,f
      integer :: i,n
      n=size(x); allocate(lx(n),y(n)); lx=log(x)
      do i=1,n
         f=real(count(x<=x(i)),dp)/real(n,dp)-0.5_dp/real(n,dp)
         if(f<=0.0_dp .or. f>=1.0_dp) then; est%status=3; est%message='invalid PP empirical transform'; return; end if
         y(i)=-log(-log(f))
      end do
      xb=sum(lx)/real(n,dp); yb=sum(y)/real(n,dp); sxx=sum((lx-xb)**2)
      if(sxx<=tiny(1.0_dp)) then; est%status=4; est%message='PP regression has zero x variance'; return; end if
      sxy=sum((lx-xb)*(y-yb)); est%beta=sxy/sxx
      est%q=exp(-exp(-(yb-est%beta*xb))); est%objective=loglikediw(x,est%q,est%beta)
   end function

   function estdiweibull(x,method,control) result(est)
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in), optional :: method
      type(diw_control), intent(in), optional :: control
      type(diw_estimate) :: est
      type(diw_control) :: ctl
      character(len=2) :: m
      integer :: n,y,z
      real(dp) :: m1,q0,beta0
      type(solnp_problem) :: prob
      type(solnp_control) :: sc
      type(solnp_result) :: sr
      type(moment_context) :: ctx
      ctl=diw_control(); if(present(control)) ctl=control
      m='P '; if(present(method)) m=adjustl(method)
      n=size(x)
      if(n==0 .or. any(x<1.0_dp) .or. any(abs(x-anint(x))>1.0e-10_dp)) then
         est%status=10; est%message='sample must contain positive integers'; return
      end if
      select case(trim(m))
      case('P','p')
         y=count(nint(x)==1)
         if(y==0) then
            est%q=ieee_value(0.0_dp,ieee_quiet_nan); est%beta=est%q
            est%status=11; est%message='proportion method requires at least one 1'; return
         end if
         est%q=real(y,dp)/real(n,dp); z=count(x<=2.0_dp)
         if(z==y .or. z==n) then
            est%beta=ieee_value(0.0_dp,ieee_quiet_nan)
            est%status=12; est%message='proportion beta estimate unavailable'; return
         end if
         est%beta=-log(log(real(z,dp)/real(n,dp))/log(real(y,dp)/real(n,dp)))/log2
         est%objective=loglikediw(x,est%q,est%beta)
      case('M','m')
         if(count(x<=2.0_dp)==n) then; est%status=13; est%message='moment method requires observations above 2'; return; end if
         m1=sum(x)/real(n,dp); q0=max(1.0e-4_dp,min(0.9999_dp,(m1-1.0_dp)/m1)); beta0=4.0_dp
         prob%name='DiscreteInverseWeibull moments'; prob%n=2; prob%fn=>mom_objective
         prob%start=[q0,beta0]; prob%lower=[0.0_dp,2.0_dp]; prob%upper=[1.0_dp,100.0_dp]
         ctx%x=x; ctx%eps=ctl%eps; ctx%nmax=ctl%nmax; allocate(prob%data,source=ctx)
         sc%max_iter=ctl%solnp_max_iter; sc%tol=ctl%solnp_tol; sc%trace=ctl%trace
         call solnp(prob,sr,sc); est%q=sr%pars(1); est%beta=sr%pars(2); est%objective=sr%objective
         est%iterations=sr%out_iterations; est%status=sr%convergence; est%message=sr%message
      case('H','h')
         est=heuristic(x,ctl%beta1,ctl%z,ctl%r,ctl%leps,ctl%heuristic_max_iter)
      case('PP','pp')
         est=estimate_pp(x)
      case default
         est%status=14; est%message='unknown estimation method'
      end select
   end function

   subroutine set_rng_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: s(:)
      call random_seed(size=n); allocate(s(n))
      do i=1,n; s(i)=mod(abs(seed)+104729*i,2147483646)+1; end do
      call random_seed(put=s)
   end subroutine
end module discrete_inverse_weibull
