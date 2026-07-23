! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_stats
   use tseries_kinds, only : dp
   implicit none
   private

   public :: mean_value
   public :: variance_value
   public :: standard_deviation
   public :: autocovariance
   public :: long_run_variance
   public :: linear_interpolate
   public :: cumulative_sum
   public :: rank_values
   public :: sample_quantile

contains

   pure real(dp) function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x)/real(size(x),dp)
      end if
   end function mean_value

   pure real(dp) function variance_value(x, sample) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: sample
      logical :: use_sample
      real(dp) :: mu
      integer :: denom
      use_sample=.true.
      if(present(sample)) use_sample=sample
      if(size(x)<2 .and. use_sample) then
         value=0.0_dp
         return
      end if
      if(size(x)==0) then
         value=0.0_dp
         return
      end if
      mu=mean_value(x)
      denom=size(x)
      if(use_sample) denom=denom-1
      value=sum((x-mu)**2)/real(max(1,denom),dp)
   end function variance_value

   pure real(dp) function standard_deviation(x, sample) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: sample
      if(present(sample)) then
         value=sqrt(max(0.0_dp,variance_value(x,sample)))
      else
         value=sqrt(max(0.0_dp,variance_value(x)))
      end if
   end function standard_deviation

   pure real(dp) function autocovariance(x, lag, demean, divisor_n) result(value)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lag
      logical, intent(in), optional :: demean, divisor_n
      logical :: dm, dn
      real(dp) :: mu
      integer :: n
      n=size(x); dm=.true.; dn=.true.
      if(present(demean)) dm=demean
      if(present(divisor_n)) dn=divisor_n
      if(lag<0 .or. lag>=n) then
         value=0.0_dp; return
      end if
      mu=0.0_dp
      if(dm) mu=mean_value(x)
      value=sum((x(1+lag:n)-mu)*(x(1:n-lag)-mu))
      if(dn) then
         value=value/real(n,dp)
      else
         value=value/real(n-lag,dp)
      end if
   end function autocovariance

   pure real(dp) function long_run_variance(u, lag) result(value)
      real(dp), intent(in) :: u(:)
      integer, intent(in) :: lag
      integer :: j,n,l
      n=size(u); l=max(0,min(lag,n-1))
      value=sum(u*u)/real(max(1,n),dp)
      do j=1,l
         value=value+2.0_dp*(1.0_dp-real(j,dp)/real(l+1,dp))* &
            sum(u(1+j:n)*u(1:n-j))/real(n,dp)
      end do
   end function long_run_variance

   pure real(dp) function linear_interpolate(xp, yp, x, clamp) result(y)
      real(dp), intent(in) :: xp(:), yp(:), x
      logical, intent(in), optional :: clamp
      logical :: use_clamp,ascending
      integer :: n,i
      use_clamp=.true.
      if(present(clamp)) use_clamp=clamp
      n=size(xp)
      if(n==0 .or. size(yp)/=n) then
         y=0.0_dp; return
      end if
      if(n==1) then
         y=yp(1); return
      end if
      ascending=xp(n)>=xp(1)
      if((ascending .and. x<=xp(1)) .or. (.not.ascending .and. x>=xp(1))) then
         if(use_clamp) then
            y=yp(1)
         else
            y=yp(1)+(x-xp(1))*(yp(2)-yp(1))/(xp(2)-xp(1))
         end if
         return
      end if
      if((ascending .and. x>=xp(n)) .or. (.not.ascending .and. x<=xp(n))) then
         if(use_clamp) then
            y=yp(n)
         else
            y=yp(n-1)+(x-xp(n-1))*(yp(n)-yp(n-1))/(xp(n)-xp(n-1))
         end if
         return
      end if
      do i=1,n-1
         if((x>=xp(i) .and. x<=xp(i+1)) .or. (x<=xp(i) .and. x>=xp(i+1))) then
            y=yp(i)+(x-xp(i))*(yp(i+1)-yp(i))/(xp(i+1)-xp(i))
            return
         end if
      end do
      y=yp(n)
   end function linear_interpolate

   pure subroutine cumulative_sum(x, s)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: s(:)
      integer :: i
      if(size(s)/=size(x)) return
      if(size(x)==0) return
      s(1)=x(1)
      do i=2,size(x)
         s(i)=s(i-1)+x(i)
      end do
   end subroutine cumulative_sum

   subroutine rank_values(x, ranks)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: ranks(:)
      integer, allocatable :: idx(:)
      integer :: i,j,key,n
      n=size(x)
      if(size(ranks)/=n) return
      allocate(idx(n))
      idx=[(i,i=1,n)]
      do i=2,n
         key=idx(i); j=i-1
         do while(j>=1)
            if(x(idx(j))<=x(key)) exit
            idx(j+1)=idx(j); j=j-1
         end do
         idx(j+1)=key
      end do
      do i=1,n
         ranks(idx(i))=i
      end do
   end subroutine rank_values

   real(dp) function sample_quantile(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp), allocatable :: work(:)
      real(dp) :: h,frac,tmp
      integer :: i,j,n,k
      n=size(x)
      if(n==0) then
         value=0.0_dp; return
      end if
      allocate(work(n)); work=x
      do i=2,n
         tmp=work(i); j=i-1
         do while(j>=1)
            if(work(j)<=tmp) exit
            work(j+1)=work(j); j=j-1
         end do
         work(j+1)=tmp
      end do
      h=1.0_dp+(real(n-1,dp))*max(0.0_dp,min(1.0_dp,probability))
      k=int(floor(h)); frac=h-real(k,dp)
      if(k>=n) then
         value=work(n)
      else
         value=(1.0_dp-frac)*work(k)+frac*work(k+1)
      end if
   end function sample_quantile

end module tseries_stats
