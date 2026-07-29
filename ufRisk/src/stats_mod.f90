! SPDX-License-Identifier: GPL-3.0-only
module stats_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: variance, normal_quantile
contains
   pure real(dp) function variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: center
      integer :: n
      n = size(x)
      if (n < 2) then
         value = 0.0_dp
      else
         center = sum(x)/real(n, dp)
         value = sum((x-center)**2)/real(n-1, dp)
      end if
   end function variance

   pure real(dp) function normal_quantile(p) result(x)
      !! Acklam inverse-normal approximation with one Newton correction.
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, &
          2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, &
          2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r, err, density
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp); return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp); return
      end if
      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp; r=q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      err = 0.5_dp*erfc(-x/sqrt(2.0_dp))-p
      density = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
      if (density > tiny(1.0_dp)) x = x-err/density
   end function normal_quantile
end module stats_mod
