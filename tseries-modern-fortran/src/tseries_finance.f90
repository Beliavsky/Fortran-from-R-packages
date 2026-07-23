! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_finance
   use tseries_kinds, only : dp
   use tseries_types, only : drawdown_result, portfolio_result
   use tseries_stats, only : mean_value, standard_deviation
   use tseries_linalg, only : covariance_matrix, solve_linear
   implicit none
   private

   public :: maximum_drawdown
   public :: sharpe_ratio
   public :: sterling_ratio
   public :: portfolio_optimize

contains

   function maximum_drawdown(x) result(drawdown)
      real(dp), intent(in) :: x(:)
      type(drawdown_result) :: drawdown
      real(dp) :: peak,current
      integer :: peak_index,i
      if(size(x)==0) return
      peak=x(1); peak_index=1
      do i=1,size(x)
         if(x(i)>peak) then
            peak=x(i); peak_index=i
         end if
         current=peak-x(i)
         if(current>drawdown%maximum) then
            drawdown%maximum=current
            drawdown%from_index=peak_index
            drawdown%to_index=i
         end if
      end do
   end function maximum_drawdown

   real(dp) function sharpe_ratio(x,risk_free_change,scale) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: risk_free_change,scale
      real(dp), allocatable :: change(:)
      real(dp) :: rf,scl,sd
      rf=0.0_dp; scl=sqrt(250.0_dp)
      if(present(risk_free_change)) rf=risk_free_change
      if(present(scale)) scl=scale
      if(size(x)<2) then
         value=0.0_dp; return
      end if
      allocate(change(size(x)-1)); change=x(2:)-x(:size(x)-1)
      sd=standard_deviation(change)
      if(sd<=tiny(1.0_dp)) then
         value=0.0_dp
      else
         value=scl*(mean_value(change)-rf)/sd
      end if
   end function sharpe_ratio

   real(dp) function sterling_ratio(x) result(value)
      real(dp), intent(in) :: x(:)
      type(drawdown_result) :: dd
      if(size(x)<2) then
         value=0.0_dp; return
      end if
      dd=maximum_drawdown(x)
      if(dd%maximum<=tiny(1.0_dp)) then
         value=huge(1.0_dp)
      else
         value=(x(size(x))-x(1))/dd%maximum
      end if
   end function sterling_ratio

   function portfolio_optimize(returns,target_return,riskless,shorts,risk_free_rate,covariance) result(result)
      real(dp), intent(in) :: returns(:, :)
      real(dp), intent(in), optional :: target_return,risk_free_rate
      logical, intent(in), optional :: riskless,shorts
      real(dp), intent(in), optional :: covariance(:, :)
      type(portfolio_result) :: result
      real(dp), allocatable :: mu(:),cov(:,:),weights(:)
      logical, allocatable :: active(:)
      logical :: has_riskless,allow_shorts
      real(dp) :: target,rf
      integer :: n,k,status,iteration,negative_index,i

      n=size(returns,1); k=size(returns,2)
      has_riskless=.false.; allow_shorts=.false.; rf=0.0_dp
      if(present(riskless)) has_riskless=riskless
      if(present(shorts)) allow_shorts=shorts
      if(present(risk_free_rate)) rf=risk_free_rate
      allocate(mu(k),cov(k,k),weights(k),active(k))
      do i=1,k
         mu(i)=mean_value(returns(:,i))
      end do
      target=sum(returns)/real(max(1,n*k),dp)
      if(present(target_return)) target=target_return
      if(present(covariance)) then
         if(size(covariance,1)/=k .or. size(covariance,2)/=k) then
            result%status=1; result%message='invalid covariance dimensions'; return
         end if
         cov=covariance
      else
         call covariance_matrix(returns,cov)
      end if
      active=.true.
      do iteration=1,k
         call solve_active(cov,mu,target,rf,has_riskless,active,weights,status)
         if(status/=0) then
            result%status=status; result%message='singular portfolio optimization system'; return
         end if
         if(allow_shorts .or. minval(weights,mask=active)>=-1.0e-10_dp) exit
         negative_index=0
         do i=1,k
            if(active(i) .and. weights(i)<-1.0e-10_dp) then
               if(negative_index==0 .or. weights(i)<weights(negative_index)) negative_index=i
            end if
         end do
         if(negative_index==0) exit
         active(negative_index)=.false.; weights(negative_index)=0.0_dp
         if(count(active)<merge(1,2,has_riskless)) then
            result%status=3; result%message='constraints are infeasible without short sales'; return
         end if
      end do
      where(.not.active) weights=0.0_dp
      allocate(result%weights(k)); result%weights=weights
      result%expected_return=dot_product(mu,weights)+merge(rf*(1.0_dp-sum(weights)),0.0_dp,has_riskless)
      result%standard_deviation=sqrt(max(0.0_dp,dot_product(weights,matmul(cov,weights))))
      result%status=0; result%message='ok'
   end function portfolio_optimize

   subroutine solve_active(cov,mu,target,rf,riskless,active,weights,status)
      real(dp), intent(in) :: cov(:, :),mu(:),target,rf
      logical, intent(in) :: riskless,active(:)
      real(dp), intent(out) :: weights(:)
      integer, intent(out) :: status
      integer, allocatable :: idx(:)
      real(dp), allocatable :: kkt(:,:),rhs(:),solution(:)
      integer :: ka,m,i,j,row
      ka=count(active); m=merge(1,2,riskless)
      allocate(idx(ka)); row=0
      do i=1,size(active)
         if(active(i)) then
            row=row+1; idx(row)=i
         end if
      end do
      allocate(kkt(ka+m,ka+m),rhs(ka+m),solution(ka+m))
      kkt=0.0_dp; rhs=0.0_dp
      do i=1,ka
         do j=1,ka
            kkt(i,j)=cov(idx(i),idx(j))
         end do
      end do
      if(riskless) then
         do i=1,ka
            kkt(i,ka+1)=mu(idx(i))-rf
            kkt(ka+1,i)=kkt(i,ka+1)
         end do
         rhs(ka+1)=target-rf
      else
         do i=1,ka
            kkt(i,ka+1)=1.0_dp; kkt(ka+1,i)=1.0_dp
            kkt(i,ka+2)=mu(idx(i)); kkt(ka+2,i)=mu(idx(i))
         end do
         rhs(ka+1)=1.0_dp; rhs(ka+2)=target
      end if
      call solve_linear(kkt,rhs,solution,status)
      weights=0.0_dp
      if(status==0) then
         do i=1,ka
            weights(idx(i))=solution(i)
         end do
      end if
   end subroutine solve_active

end module tseries_finance
