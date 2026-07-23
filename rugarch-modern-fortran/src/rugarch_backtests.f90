! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_backtests
   use rugarch_kinds, only : dp
   use rugarch_math, only : normal_cdf, regularized_gamma_p
   implicit none
   private

   type, public :: var_test_result
      integer :: observations = 0
      integer :: expected_exceedances = 0
      integer :: actual_exceedances = 0
      real(dp) :: uc_statistic = 0.0_dp
      real(dp) :: uc_p_value = 1.0_dp
      real(dp) :: cc_statistic = 0.0_dp
      real(dp) :: cc_p_value = 1.0_dp
   end type var_test_result

   type, public :: directional_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: directional_accuracy = 0.0_dp
   end type directional_test_result

   type, public :: berkowitz_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: sd = 1.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: lr_statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
   end type berkowitz_result

   public :: var_test, es_test_p_value, directional_accuracy_test
   public :: berkowitz_test, var_loss, quantile_loss

contains

   function var_test(alpha, actual, var_forecast) result(ans)
      real(dp), intent(in) :: alpha, actual(:), var_forecast(:)
      type(var_test_result) :: ans
      integer, allocatable :: hit(:)
      integer :: n, i, n00, n01, n10, n11
      real(dp) :: phat, p01, p11, ll0, ll1, lli, eps

      n=min(size(actual),size(var_forecast))
      allocate(hit(n))
      hit=merge(1,0,actual(1:n)<var_forecast(1:n))
      ans%observations=n
      ans%expected_exceedances=int(alpha*real(n,dp))
      ans%actual_exceedances=sum(hit)
      eps=1.0e-12_dp
      phat=max(eps,min(1.0_dp-eps,real(sum(hit),dp)/real(max(1,n),dp)))
      ll0=real(sum(hit),dp)*log(max(alpha,eps))+real(n-sum(hit),dp)*log(max(1.0_dp-alpha,eps))
      ll1=real(sum(hit),dp)*log(phat)+real(n-sum(hit),dp)*log(1.0_dp-phat)
      ans%uc_statistic=max(0.0_dp,2.0_dp*(ll1-ll0))
      ans%uc_p_value=chi_square_sf(ans%uc_statistic,1.0_dp)

      n00=0; n01=0; n10=0; n11=0
      do i=2,n
         if (hit(i-1)==0 .and. hit(i)==0) n00=n00+1
         if (hit(i-1)==0 .and. hit(i)==1) n01=n01+1
         if (hit(i-1)==1 .and. hit(i)==0) n10=n10+1
         if (hit(i-1)==1 .and. hit(i)==1) n11=n11+1
      end do
      p01=max(eps,min(1.0_dp-eps,real(n01,dp)/real(max(1,n00+n01),dp)))
      p11=max(eps,min(1.0_dp-eps,real(n11,dp)/real(max(1,n10+n11),dp)))
      lli=real(n00,dp)*log(1.0_dp-p01)+real(n01,dp)*log(p01)+ &
          real(n10,dp)*log(1.0_dp-p11)+real(n11,dp)*log(p11)
      ans%cc_statistic=max(0.0_dp,ans%uc_statistic+2.0_dp*(lli-ll1))
      ans%cc_p_value=chi_square_sf(ans%cc_statistic,2.0_dp)
   end function var_test

   function es_test_p_value(actual, es_forecast, var_forecast) result(p_value)
      real(dp), intent(in) :: actual(:), es_forecast(:), var_forecast(:)
      real(dp) :: p_value
      real(dp), allocatable :: z(:)
      real(dp) :: meanz, sdz, statistic
      integer :: i, n, k
      n=min(size(actual),min(size(es_forecast),size(var_forecast)))
      k=count(actual(1:n)<var_forecast(1:n))
      if (k<2) then
         p_value=1.0_dp
         return
      end if
      allocate(z(k))
      k=0
      do i=1,n
         if (actual(i)<var_forecast(i)) then
            k=k+1
            z(k)=es_forecast(i)-actual(i)
         end if
      end do
      meanz=sum(z)/real(k,dp)
      sdz=sqrt(sum((z-meanz)**2)/real(k-1,dp))
      statistic=meanz/max(sdz/sqrt(real(k,dp)),1.0e-20_dp)
      p_value=1.0_dp-normal_cdf(statistic)
   end function es_test_p_value

   function directional_accuracy_test(forecast, actual) result(ans)
      real(dp), intent(in) :: forecast(:), actual(:)
      type(directional_test_result) :: ans
      real(dp) :: px, py, phat, pstar, vhat, vstar
      integer :: n
      n=min(size(forecast),size(actual))
      if (n<2) return
      px=real(count(actual(1:n)>0.0_dp),dp)/real(n,dp)
      py=real(count(forecast(1:n)>0.0_dp),dp)/real(n,dp)
      phat=real(count(forecast(1:n)*actual(1:n)>0.0_dp),dp)/real(n,dp)
      pstar=py*px+(1.0_dp-py)*(1.0_dp-px)
      vhat=pstar*(1.0_dp-pstar)/real(n,dp)
      vstar=((2.0_dp*py-1.0_dp)**2*px*(1.0_dp-px)+(2.0_dp*px-1.0_dp)**2*py*(1.0_dp-py))/real(n,dp) + &
            4.0_dp*px*py*(1.0_dp-px)*(1.0_dp-py)/real(n*n,dp)
      ans%statistic=(phat-pstar)/sqrt(max(vhat-vstar,1.0e-20_dp))
      ans%p_value=1.0_dp-normal_cdf(ans%statistic)
      ans%directional_accuracy=phat
   end function directional_accuracy_test

   function berkowitz_test(z) result(ans)
      real(dp), intent(in) :: z(:)
      type(berkowitz_result) :: ans
      real(dp) :: xbar, denom, numer, variance, ull, rll, pi
      integer :: n
      n=size(z)
      if (n<3) return
      xbar=sum(z)/real(n,dp)
      denom=sum((z(1:n-1)-xbar)**2)
      numer=sum((z(1:n-1)-xbar)*(z(2:n)-xbar))
      ans%rho=numer/max(denom,1.0e-20_dp)
      ans%mean=xbar*(1.0_dp-ans%rho)
      variance=sum((z(2:n)-ans%mean-ans%rho*z(1:n-1))**2)/real(n-1,dp)
      ans%sd=sqrt(max(variance,1.0e-20_dp))
      pi=acos(-1.0_dp)
      ull=-real(n-1,dp)*log(ans%sd)-0.5_dp*real(n-1,dp)*log(2.0_dp*pi) - &
          0.5_dp*sum(((z(2:n)-ans%mean-ans%rho*z(1:n-1))/ans%sd)**2)
      rll=-0.5_dp*real(n-1,dp)*log(2.0_dp*pi)-0.5_dp*sum(z(2:n)**2)
      ans%lr_statistic=max(0.0_dp,2.0_dp*(ull-rll))
      ans%p_value=chi_square_sf(ans%lr_statistic,3.0_dp)
   end function berkowitz_test

   pure elemental function var_loss(actual, var_forecast, alpha) result(value)
      real(dp), intent(in) :: actual, var_forecast, alpha
      real(dp) :: value
      if (actual<var_forecast) then
         value=1.0_dp+(actual-var_forecast)**2
      else
         value=-alpha
      end if
   end function var_loss

   pure elemental function quantile_loss(actual, quantile_forecast, alpha) result(value)
      real(dp), intent(in) :: actual, quantile_forecast, alpha
      real(dp) :: value
      value=(alpha-merge(1.0_dp,0.0_dp,actual<quantile_forecast))*(actual-quantile_forecast)
   end function quantile_loss

   pure elemental function chi_square_sf(x, df) result(value)
      real(dp), intent(in) :: x, df
      real(dp) :: value
      value=1.0_dp-regularized_gamma_p(0.5_dp*df,0.5_dp*max(x,0.0_dp))
      value=max(0.0_dp,min(1.0_dp,value))
   end function chi_square_sf

end module rugarch_backtests
