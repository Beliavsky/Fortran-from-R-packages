! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_curves
  use rq_kinds, only: dp
  use rq_linalg, only: least_squares_solve
  implicit none
  private
  public :: discount_curve_t, curve_fit_result
  public :: make_flat_curve, make_zero_curve, bootstrap_curve
  public :: fit_nelson_siegel, fit_svensson, nelson_siegel_rate, svensson_rate

  type :: discount_curve_t
    real(dp), allocatable :: times(:)
    real(dp), allocatable :: discounts(:)
  contains
    procedure :: discount => curve_discount
    procedure :: zero_rate => curve_zero_rate
    procedure :: forward_rate => curve_forward_rate
  end type discount_curve_t

  type :: curve_fit_result
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: sse=0.0_dp
    integer :: status=1
  end type curve_fit_result
contains
  subroutine make_flat_curve(rate,max_time,dt,curve)
    real(dp),intent(in)::rate,max_time,dt
    type(discount_curve_t),intent(out)::curve
    integer::n,i
    n=max(2,int(max_time/dt)+1); allocate(curve%times(n),curve%discounts(n))
    do i=1,n; curve%times(i)=real(i-1,dp)*dt; curve%discounts(i)=exp(-rate*curve%times(i)); end do
  end subroutine make_flat_curve

  subroutine make_zero_curve(times,zero_rates,curve)
    real(dp),intent(in)::times(:),zero_rates(:)
    type(discount_curve_t),intent(out)::curve
    integer::n
    n=min(size(times),size(zero_rates)); allocate(curve%times(n),curve%discounts(n))
    curve%times=times(1:n); curve%discounts=exp(-zero_rates(1:n)*times(1:n))
    if(curve%times(1)>0.0_dp) then
      curve%discounts(1)=exp(-zero_rates(1)*curve%times(1))
    end if
  end subroutine make_zero_curve

  pure real(dp) function curve_discount(self,t) result(df)
    class(discount_curve_t),intent(in)::self
    real(dp),intent(in)::t
    integer::n,i
    real(dp)::w,ld1,ld2
    n=size(self%times)
    if(t<=0.0_dp) then; df=1.0_dp; return; end if
    if(t<=self%times(1)) then
      if(self%times(1)<=0.0_dp) then; df=self%discounts(1); else; df=exp(log(self%discounts(1))*t/self%times(1)); end if
      return
    end if
    if(t>=self%times(n)) then
      if(n==1) then; df=exp(log(self%discounts(n))*t/self%times(n)); else
        ld1=log(self%discounts(n-1)); ld2=log(self%discounts(n))
        df=exp(ld2+(t-self%times(n))*(ld2-ld1)/(self%times(n)-self%times(n-1)))
      end if
      return
    end if
    do i=1,n-1
      if(t<=self%times(i+1)) then
        w=(t-self%times(i))/(self%times(i+1)-self%times(i))
        df=exp((1.0_dp-w)*log(self%discounts(i))+w*log(self%discounts(i+1))); return
      end if
    end do
    df=self%discounts(n)
  end function curve_discount

  pure real(dp) function curve_zero_rate(self,t) result(r)
    class(discount_curve_t),intent(in)::self
    real(dp),intent(in)::t
    if(t<=1.0e-12_dp) then
      r=-log(self%discount(max(self%times(1),1.0e-6_dp)))/max(self%times(1),1.0e-6_dp)
    else
      r=-log(self%discount(t))/t
    end if
  end function curve_zero_rate

  pure real(dp) function curve_forward_rate(self,t1,t2) result(f)
    class(discount_curve_t),intent(in)::self
    real(dp),intent(in)::t1,t2
    if(t2<=t1) then; f=0.0_dp; else; f=log(self%discount(t1)/self%discount(t2))/(t2-t1); end if
  end function curve_forward_rate

  subroutine bootstrap_curve(deposit_times,deposit_rates,future_starts,future_ends,future_prices, &
                             swap_maturities,swap_rates,swap_frequency,curve,status)
    real(dp),intent(in)::deposit_times(:),deposit_rates(:)
    real(dp),intent(in),optional::future_starts(:),future_ends(:),future_prices(:)
    real(dp),intent(in),optional::swap_maturities(:),swap_rates(:)
    integer,intent(in),optional::swap_frequency
    type(discount_curve_t),intent(out)::curve
    integer,intent(out)::status
    real(dp),allocatable::t(:),d(:)
    integer::maxn,n,i,j,freq,nc
    real(dp)::rate,dt,mat,coupon,sumdf,lastdf
    maxn=size(deposit_times)
    if(present(future_starts)) maxn=maxn+size(future_starts)
    freq=2
    if(present(swap_frequency)) freq=max(1,swap_frequency)
    if(present(swap_maturities)) then
      maxn=maxn+sum(max(1,nint(swap_maturities*real(freq,dp))))
    end if
    allocate(t(maxn+1),d(maxn+1)); n=1; t(1)=0.0_dp; d(1)=1.0_dp
    do i=1,min(size(deposit_times),size(deposit_rates))
      call insert_point(deposit_times(i),1.0_dp/(1.0_dp+deposit_rates(i)*deposit_times(i)))
    end do
    if(present(future_starts).and.present(future_ends).and.present(future_prices)) then
      do i=1,min(size(future_starts),size(future_ends),size(future_prices))
        rate=(100.0_dp-future_prices(i))/100.0_dp; dt=future_ends(i)-future_starts(i)
        call insert_point(future_ends(i),interp_discount(future_starts(i))/(1.0_dp+rate*dt))
      end do
    end if
    if(present(swap_maturities).and.present(swap_rates)) then
      do i=1,min(size(swap_maturities),size(swap_rates))
        mat=swap_maturities(i); nc=max(1,nint(mat*real(freq,dp))); coupon=swap_rates(i)/real(freq,dp); sumdf=0.0_dp
        do j=1,nc-1; sumdf=sumdf+interp_discount(real(j,dp)/real(freq,dp)); end do
        lastdf=(1.0_dp-coupon*sumdf)/(1.0_dp+coupon)
        call insert_point(mat,max(lastdf,1.0e-10_dp))
      end do
    end if
    allocate(curve%times(n),curve%discounts(n)); curve%times=t(1:n); curve%discounts=d(1:n); status=0
  contains
    subroutine insert_point(tt,dd)
      real(dp),intent(in)::tt,dd
      integer::k,pos
      pos=n+1
      do k=1,n
        if(abs(t(k)-tt)<1.0e-12_dp) then; d(k)=dd; return; end if
        if(t(k)>tt) then; pos=k; exit; end if
      end do
      n=n+1
      do k=n,pos+1,-1; t(k)=t(k-1); d(k)=d(k-1); end do
      t(pos)=tt; d(pos)=dd
    end subroutine insert_point
    function interp_discount(tt) result(df)
      real(dp),intent(in)::tt
      real(dp)::df,w
      integer::k
      if(tt<=t(1)) then; df=1.0_dp; return; end if
      if(tt>=t(n)) then; df=exp(log(d(n))*tt/max(t(n),1.0e-12_dp)); return; end if
      do k=1,n-1
        if(tt<=t(k+1)) then; w=(tt-t(k))/(t(k+1)-t(k)); df=exp((1.0_dp-w)*log(d(k))+w*log(d(k+1))); return; end if
      end do
      df=d(n)
    end function interp_discount
  end subroutine bootstrap_curve

  pure elemental real(dp) function nelson_siegel_rate(t,b0,b1,b2,tau) result(r)
    real(dp),intent(in)::t,b0,b1,b2,tau
    real(dp)::x,l1,l2
    if(t<=1.0e-12_dp) then; r=b0+b1; return; end if
    x=t/tau; l1=(1.0_dp-exp(-x))/x; l2=l1-exp(-x); r=b0+b1*l1+b2*l2
  end function nelson_siegel_rate

  pure elemental real(dp) function svensson_rate(t,b0,b1,b2,b3,tau1,tau2) result(r)
    real(dp),intent(in)::t,b0,b1,b2,b3,tau1,tau2
    real(dp)::x1,x2,l1,l2,l3
    if(t<=1.0e-12_dp) then; r=b0+b1; return; end if
    x1=t/tau1; x2=t/tau2; l1=(1.0_dp-exp(-x1))/x1; l2=l1-exp(-x1); l3=(1.0_dp-exp(-x2))/x2-exp(-x2)
    r=b0+b1*l1+b2*l2+b3*l3
  end function svensson_rate

  subroutine fit_nelson_siegel(times,rates,result)
    real(dp),intent(in)::times(:),rates(:)
    type(curve_fit_result),intent(out)::result
    real(dp)::tau,best_tau,beta(3),sse,best_sse,step
    integer::i,k
    best_sse=huge(1.0_dp); best_tau=1.0_dp
    do i=0,200
      tau=exp(log(0.03_dp)+real(i,dp)/200.0_dp*(log(30.0_dp)-log(0.03_dp)))
      call linear_beta_ns(times,rates,tau,beta,sse)
      if(sse<best_sse) then; best_sse=sse; best_tau=tau; end if
    end do
    step=0.2_dp
    do k=1,50
      tau=max(0.005_dp,best_tau*exp(-step)); call linear_beta_ns(times,rates,tau,beta,sse)
      if(sse<best_sse) then; best_sse=sse; best_tau=tau; cycle; end if
      tau=best_tau*exp(step); call linear_beta_ns(times,rates,tau,beta,sse)
      if(sse<best_sse) then; best_sse=sse; best_tau=tau; cycle; end if
      step=step*0.7_dp
    end do
    call linear_beta_ns(times,rates,best_tau,beta,best_sse)
    allocate(result%parameters(4),result%fitted(size(times)),result%residuals(size(times)))
    result%parameters=[beta,best_tau]; result%fitted=nelson_siegel_rate(times,beta(1),beta(2),beta(3),best_tau)
    result%residuals=rates-result%fitted; result%sse=sum(result%residuals**2); result%status=0
  end subroutine fit_nelson_siegel

  subroutine fit_svensson(times,rates,result)
    real(dp),intent(in)::times(:),rates(:)
    type(curve_fit_result),intent(out)::result
    real(dp)::tau1,tau2,best1,best2,beta(4),sse,best_sse,step
    integer::i,j,k
    best_sse=huge(1.0_dp); best1=1.0_dp; best2=3.0_dp
    do i=0,35; tau1=exp(log(0.05_dp)+real(i,dp)/35.0_dp*(log(10.0_dp)-log(0.05_dp)))
      do j=0,35; tau2=exp(log(0.1_dp)+real(j,dp)/35.0_dp*(log(40.0_dp)-log(0.1_dp)))
        if(abs(log(tau1/tau2))<0.05_dp) cycle
        call linear_beta_sv(times,rates,tau1,tau2,beta,sse)
        if(sse<best_sse) then; best_sse=sse; best1=tau1; best2=tau2; end if
      end do
    end do
    step=0.15_dp
    do k=1,60
      call try_pair(best1*exp(-step),best2); call try_pair(best1*exp(step),best2)
      call try_pair(best1,best2*exp(-step)); call try_pair(best1,best2*exp(step)); step=step*0.85_dp
    end do
    call linear_beta_sv(times,rates,best1,best2,beta,best_sse)
    allocate(result%parameters(6),result%fitted(size(times)),result%residuals(size(times)))
    result%parameters=[beta,best1,best2]
    result%fitted=svensson_rate(times,beta(1),beta(2),beta(3),beta(4),best1,best2)
    result%residuals=rates-result%fitted; result%sse=sum(result%residuals**2); result%status=0
  contains
    subroutine try_pair(a,b)
      real(dp),intent(in)::a,b
      real(dp)::bb(4),ss
      if(a<0.005_dp.or.b<0.005_dp.or.abs(log(a/b))<0.02_dp) return
      call linear_beta_sv(times,rates,a,b,bb,ss)
      if(ss<best_sse) then; best_sse=ss; best1=a; best2=b; end if
    end subroutine try_pair
  end subroutine fit_svensson

  subroutine linear_beta_ns(t,y,tau,beta,sse)
    real(dp),intent(in)::t(:),y(:),tau
    real(dp),intent(out)::beta(3),sse
    real(dp),allocatable::x(:,:)
    integer::i
    allocate(x(size(t),3)); x(:,1)=1.0_dp
    do i=1,size(t)
      if(t(i)<=1.0e-12_dp) then; x(i,2)=1.0_dp; x(i,3)=0.0_dp
      else; x(i,2)=(1.0_dp-exp(-t(i)/tau))/(t(i)/tau); x(i,3)=x(i,2)-exp(-t(i)/tau); end if
    end do
    call least_squares(x,y,beta); sse=sum((y-matmul(x,beta))**2)
  end subroutine linear_beta_ns

  subroutine linear_beta_sv(t,y,tau1,tau2,beta,sse)
    real(dp),intent(in)::t(:),y(:),tau1,tau2
    real(dp),intent(out)::beta(4),sse
    real(dp),allocatable::x(:,:)
    integer::i
    allocate(x(size(t),4)); x(:,1)=1.0_dp
    do i=1,size(t)
      if(t(i)<=1.0e-12_dp) then; x(i,2)=1.0_dp; x(i,3)=0.0_dp; x(i,4)=0.0_dp
      else
        x(i,2)=(1.0_dp-exp(-t(i)/tau1))/(t(i)/tau1); x(i,3)=x(i,2)-exp(-t(i)/tau1)
        x(i,4)=(1.0_dp-exp(-t(i)/tau2))/(t(i)/tau2)-exp(-t(i)/tau2)
      end if
    end do
    call least_squares(x,y,beta); sse=sum((y-matmul(x,beta))**2)
  end subroutine linear_beta_sv

  subroutine least_squares(x,y,beta)
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),intent(out)::beta(:)
    call least_squares_solve(x,y,beta)
  end subroutine least_squares
end module rq_curves
