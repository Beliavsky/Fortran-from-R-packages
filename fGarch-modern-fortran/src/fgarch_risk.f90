! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module fgarch_risk
   use fgarch_kinds, only : dp
   use fgarch_distributions, only : distribution_quantile
   implicit none
   private

   public :: value_at_risk, expected_shortfall
   public :: jarque_bera_statistic, ljung_box_statistic

contains

   pure elemental function value_at_risk(probability, location, scale, kind, shape, skew) result(value)
      real(dp), intent(in) :: probability, location, scale, shape, skew
      integer, intent(in) :: kind
      real(dp) :: value

      value = location+scale*distribution_quantile(probability,kind,shape,skew)
   end function value_at_risk

   pure function expected_shortfall(probability, location, scale, kind, shape, skew) result(value)
      real(dp), intent(in) :: probability, location, scale, shape, skew
      integer, intent(in) :: kind
      real(dp) :: value, p, dp_grid, sumq
      integer, parameter :: ngrid = 2000
      integer :: i, weight

      if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
         value = 0.0_dp
         return
      end if
      dp_grid = probability/real(ngrid,dp)
      sumq = 0.0_dp
      do i = 0, ngrid
         p = max(1.0e-12_dp,real(i,dp)*dp_grid)
         if (i == 0 .or. i == ngrid) then
            weight = 1
         else if (mod(i,2) == 0) then
            weight = 2
         else
            weight = 4
         end if
         sumq = sumq+real(weight,dp)*distribution_quantile(p,kind,shape,skew)
      end do
      value = location+scale*(sumq*dp_grid/3.0_dp)/probability
   end function expected_shortfall

   pure function jarque_bera_statistic(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, meanx, m2, m3, m4, skewness, kurtosis
      integer :: n

      n = size(x)
      meanx = sum(x)/real(n,dp)
      m2 = sum((x-meanx)**2)/real(n,dp)
      if (m2 <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      m3 = sum((x-meanx)**3)/real(n,dp)
      m4 = sum((x-meanx)**4)/real(n,dp)
      skewness = m3/m2**1.5_dp
      kurtosis = m4/(m2*m2)
      value = real(n,dp)/6.0_dp*(skewness**2+0.25_dp*(kurtosis-3.0_dp)**2)
   end function jarque_bera_statistic

   pure function ljung_box_statistic(x, lags) result(value)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lags
      real(dp) :: value, meanx, denom, rho
      integer :: n, k

      n = size(x)
      meanx = sum(x)/real(n,dp)
      denom = sum((x-meanx)**2)
      value = 0.0_dp
      if (denom <= 0.0_dp) return
      do k = 1, min(lags,n-1)
         rho = sum((x(1:n-k)-meanx)*(x(1+k:n)-meanx))/denom
         value = value+rho*rho/real(n-k,dp)
      end do
      value = real(n*(n+2),dp)*value
   end function ljung_box_statistic

end module fgarch_risk
