! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_arfima
   use rugarch_kinds, only : dp
   use rugarch_rng, only : random_normal
   implicit none
   private

   type, public :: arfima_spec
      real(dp) :: mean = 0.0_dp
      real(dp) :: d = 0.0_dp
      real(dp) :: innovation_sd = 1.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
   end type arfima_spec

   type, public :: arfima_fit_result
      type(arfima_spec) :: spec
      real(dp) :: residual_sd = 0.0_dp
      real(dp) :: css = huge(1.0_dp)
      integer :: status = 1
      real(dp), allocatable :: residuals(:)
   end type arfima_fit_result

   type, public :: arfima_forecast_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: sigma(:)
   end type arfima_forecast_result

   type, public :: arfima_order_result
      type(arfima_fit_result), allocatable :: candidates(:)
      real(dp), allocatable :: criterion(:)
      integer, allocatable :: p(:), q(:)
      integer :: best = 0
   end type arfima_order_result

   public :: make_arfima_spec, fractional_weights
   public :: fractional_difference, fractional_integrate
   public :: arma_residuals, simulate_arfima, fit_arfima
   public :: gph_estimate_d, forecast_arfima, auto_arfima

contains

   function make_arfima_spec(p, q) result(spec)
      integer, intent(in), optional :: p, q
      type(arfima_spec) :: spec
      integer :: pp, qq
      pp=0; qq=0
      if (present(p)) pp=max(0,p)
      if (present(q)) qq=max(0,q)
      allocate(spec%ar(pp),spec%ma(qq))
      if (pp>0) spec%ar=0.0_dp
      if (qq>0) spec%ma=0.0_dp
   end function make_arfima_spec

   pure subroutine fractional_weights(d, weights, inverse)
      real(dp), intent(in) :: d
      real(dp), intent(out) :: weights(:)
      logical, intent(in), optional :: inverse
      logical :: inv
      integer :: k
      real(dp) :: dd
      if (size(weights)==0) return
      inv=.false.
      if (present(inverse)) inv=inverse
      dd=merge(-d,d,inv)
      weights(1)=1.0_dp
      do k=2,size(weights)
         weights(k)=weights(k-1)*(real(k-2,dp)-dd)/real(k-1,dp)
      end do
   end subroutine fractional_weights

   subroutine fractional_difference(x, d, y, truncation)
      real(dp), intent(in) :: x(:), d
      real(dp), intent(out) :: y(size(x))
      integer, intent(in), optional :: truncation
      real(dp), allocatable :: w(:)
      integer :: n, m, i, k
      n=size(x)
      m=min(n,1000)
      if (present(truncation)) m=min(n,max(1,truncation))
      allocate(w(m))
      call fractional_weights(d,w)
      y=0.0_dp
      do i=1,n
         do k=1,min(m,i)
            y(i)=y(i)+w(k)*x(i-k+1)
         end do
      end do
   end subroutine fractional_difference

   subroutine fractional_integrate(x, d, y, truncation)
      real(dp), intent(in) :: x(:), d
      real(dp), intent(out) :: y(size(x))
      integer, intent(in), optional :: truncation
      real(dp), allocatable :: w(:)
      integer :: n, m, i, k
      n=size(x)
      m=min(n,1000)
      if (present(truncation)) m=min(n,max(1,truncation))
      allocate(w(m))
      call fractional_weights(d,w,inverse=.true.)
      y=0.0_dp
      do i=1,n
         do k=1,min(m,i)
            y(i)=y(i)+w(k)*x(i-k+1)
         end do
      end do
   end subroutine fractional_integrate

   subroutine arma_residuals(x, spec, residuals)
      real(dp), intent(in) :: x(:)
      type(arfima_spec), intent(in) :: spec
      real(dp), intent(out) :: residuals(size(x))
      real(dp), allocatable :: xd(:)
      integer :: i, j, p, q
      allocate(xd(size(x)))
      call fractional_difference(x-spec%mean,spec%d,xd)
      p=size_or_zero(spec%ar)
      q=size_or_zero(spec%ma)
      residuals=0.0_dp
      do i=1,size(x)
         residuals(i)=xd(i)
         do j=1,min(p,i-1)
            residuals(i)=residuals(i)-spec%ar(j)*xd(i-j)
         end do
         do j=1,min(q,i-1)
            residuals(i)=residuals(i)-spec%ma(j)*residuals(i-j)
         end do
      end do
   end subroutine arma_residuals

   subroutine simulate_arfima(spec, n, x, burn_in)
      type(arfima_spec), intent(in) :: spec
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      integer, intent(in), optional :: burn_in
      real(dp), allocatable :: innovation(:), arma(:), integrated(:)
      integer :: burn, nt, i, j, p, q
      burn=500
      if (present(burn_in)) burn=max(0,burn_in)
      nt=n+burn
      allocate(innovation(nt),arma(nt),integrated(nt))
      p=size_or_zero(spec%ar)
      q=size_or_zero(spec%ma)
      innovation=0.0_dp
      arma=0.0_dp
      do i=1,nt
         innovation(i)=spec%innovation_sd*random_normal()
         arma(i)=innovation(i)
         do j=1,min(p,i-1)
            arma(i)=arma(i)+spec%ar(j)*arma(i-j)
         end do
         do j=1,min(q,i-1)
            arma(i)=arma(i)+spec%ma(j)*innovation(i-j)
         end do
      end do
      call fractional_integrate(arma,spec%d,integrated)
      x=spec%mean+integrated(burn+1:nt)
   end subroutine simulate_arfima

   function fit_arfima(x, p, q, estimate_d) result(fit)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: p, q
      logical, intent(in), optional :: estimate_d
      type(arfima_fit_result) :: fit
      real(dp), allocatable :: xd(:), r(:)
      real(dp), allocatable :: gamma(:), a(:)
      real(dp) :: meanx
      integer :: pp, qq, j, n
      logical :: use_d

      pp=0; qq=0
      if (present(p)) pp=max(0,p)
      if (present(q)) qq=max(0,q)
      use_d=.true.
      if (present(estimate_d)) use_d=estimate_d
      fit%spec=make_arfima_spec(pp,qq)
      meanx=sum(x)/real(size(x),dp)
      fit%spec%mean=meanx
      if (use_d) fit%spec%d=max(-0.49_dp,min(0.49_dp,gph_estimate_d(x)))
      allocate(xd(size(x)))
      call fractional_difference(x-meanx,fit%spec%d,xd)

      if (pp>0) then
         allocate(gamma(0:pp),a(pp))
         do j=0,pp
            gamma(j)=sum(xd(1:size(xd)-j)*xd(1+j:size(xd)))/real(size(xd),dp)
         end do
         call levinson_yule_walker(gamma,a)
         fit%spec%ar=a
      end if
      if (qq>0) fit%spec%ma=0.0_dp
      allocate(r(size(x)),fit%residuals(size(x)))
      call arma_residuals(x,fit%spec,r)
      fit%residuals=r
      n=max(1,size(x)-max(pp,qq)-1)
      fit%css=sum(r(max(pp,qq)+2:)**2)
      fit%residual_sd=sqrt(fit%css/real(n,dp))
      fit%spec%innovation_sd=fit%residual_sd
      fit%status=0
   end function fit_arfima


   function forecast_arfima(history,spec,horizon,truncation) result(ans)
      real(dp),intent(in)::history(:)
      type(arfima_spec),intent(in)::spec
      integer,intent(in)::horizon
      integer,intent(in),optional::truncation
      type(arfima_forecast_result)::ans
      real(dp),allocatable::xd(:),xdext(:),resid(:),resext(:),integrated(:),weights(:)
      integer::n,h,i,j,p,q,m
      n=size(history);p=size_or_zero(spec%ar);q=size_or_zero(spec%ma)
      m=min(n+horizon,1000);if(present(truncation))m=min(n+horizon,max(1,truncation))
      allocate(ans%mean(horizon),ans%sigma(horizon),xd(n),xdext(n+horizon),resid(n), &
         resext(n+horizon),integrated(n+horizon),weights(max(1,horizon)))
      call fractional_difference(history-spec%mean,spec%d,xd,m)
      call arma_residuals(history,spec,resid)
      xdext(1:n)=xd;resext(1:n)=resid
      do h=1,horizon
         i=n+h;xdext(i)=0.0_dp
         do j=1,min(p,i-1);xdext(i)=xdext(i)+spec%ar(j)*xdext(i-j);end do
         do j=1,min(q,i-1);xdext(i)=xdext(i)+spec%ma(j)*resext(i-j);end do
         resext(i)=0.0_dp
      end do
      call fractional_integrate(xdext,spec%d,integrated,m)
      ans%mean=spec%mean+integrated(n+1:n+horizon)
      call fractional_weights(spec%d,weights,inverse=.true.)
      do h=1,horizon
         ans%sigma(h)=spec%innovation_sd*sqrt(sum(weights(1:h)**2))
      end do
   end function forecast_arfima

   function auto_arfima(x,pmax,qmax,criterion,estimate_d) result(ans)
      real(dp),intent(in)::x(:)
      integer,intent(in)::pmax,qmax
      character(len=*),intent(in),optional::criterion
      logical,intent(in),optional::estimate_d
      type(arfima_order_result)::ans
      character(len=8)::crit
      integer::p,q,idx,ncand,npars
      logical::use_d
      real(dp)::value
      crit='bic';if(present(criterion))crit=adjustl(criterion)
      use_d=.true.;if(present(estimate_d))use_d=estimate_d
      ncand=(max(0,pmax)+1)*(max(0,qmax)+1)
      allocate(ans%candidates(ncand),ans%criterion(ncand),ans%p(ncand),ans%q(ncand))
      idx=0;ans%best=1
      do p=0,max(0,pmax)
         do q=0,max(0,qmax)
            idx=idx+1;ans%p(idx)=p;ans%q(idx)=q
            ans%candidates(idx)=fit_arfima(x,p,q,estimate_d)
            npars=2+p+q+merge(1,0,use_d)
            value=real(size(x),dp)*log(max(ans%candidates(idx)%css/real(size(x),dp),1.0e-30_dp))
            if(crit(1:1)=='a'.or.crit(1:1)=='A')then
               ans%criterion(idx)=value+2.0_dp*real(npars,dp)
            else
               ans%criterion(idx)=value+log(real(size(x),dp))*real(npars,dp)
            end if
            if(ans%criterion(idx)<ans%criterion(ans%best))ans%best=idx
         end do
      end do
   end function auto_arfima

   function gph_estimate_d(x, bandwidth) result(d)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: bandwidth
      real(dp) :: d
      real(dp), allocatable :: xx(:), lx(:), ly(:)
      real(dp) :: pi, freq, re, im, xbar, ybar, denom
      integer :: n, m, j, t
      n=size(x)
      m=max(5,int(sqrt(real(n,dp))))
      if (present(bandwidth)) m=max(3,min(bandwidth,(n-1)/2))
      allocate(xx(n),lx(m),ly(m))
      xx=x-sum(x)/real(n,dp)
      pi=acos(-1.0_dp)
      do j=1,m
         freq=2.0_dp*pi*real(j,dp)/real(n,dp)
         re=0.0_dp; im=0.0_dp
         do t=1,n
            re=re+xx(t)*cos(freq*real(t-1,dp))
            im=im-xx(t)*sin(freq*real(t-1,dp))
         end do
         lx(j)=log(max(4.0_dp*sin(0.5_dp*freq)**2,1.0e-30_dp))
         ly(j)=log(max((re*re+im*im)/real(n,dp),1.0e-30_dp))
      end do
      xbar=sum(lx)/real(m,dp)
      ybar=sum(ly)/real(m,dp)
      denom=sum((lx-xbar)**2)
      if (denom<=0.0_dp) then
         d=0.0_dp
      else
         d=-0.5_dp*sum((lx-xbar)*(ly-ybar))/denom
      end if
   end function gph_estimate_d

   subroutine levinson_yule_walker(gamma,a)
      real(dp), intent(in) :: gamma(0:)
      real(dp), intent(out) :: a(:)
      real(dp), allocatable :: old(:)
      real(dp) :: variance, reflection
      integer :: k, j, p
      p=size(a)
      allocate(old(p))
      a=0.0_dp
      variance=max(gamma(0),1.0e-20_dp)
      do k=1,p
         reflection=gamma(k)
         do j=1,k-1
            reflection=reflection-a(j)*gamma(k-j)
         end do
         reflection=reflection/variance
         old=a
         a(k)=reflection
         do j=1,k-1
            a(j)=old(j)-reflection*old(k-j)
         end do
         variance=max(variance*(1.0_dp-reflection*reflection),1.0e-20_dp)
      end do
   end subroutine levinson_yule_walker

   pure integer function size_or_zero(x) result(n)
      real(dp), allocatable, intent(in) :: x(:)
      if (allocated(x)) then
         n=size(x)
      else
         n=0
      end if
   end function size_or_zero

end module rugarch_arfima
