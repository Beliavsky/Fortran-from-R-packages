! Computational surface corresponding to DPQ's bundled TOMS Algorithm 1006 code.
! SPDX-License-Identifier: GPL-3.0-or-later
!
! The Pugh 11-term coefficients and generalized incomplete-gamma definition
! originate in Abergel & Moisan, Algorithm 1006, ACM TOMS 46(1), 2020.
module dpq_toms1006
   use r_compat, only: dp, pgamma, r_lgamma
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
      ieee_positive_inf, ieee_is_finite
   implicit none
   private

   type, public :: dltgamma_result
      real(dp) :: rho = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      integer :: method = 0
   contains
      procedure :: value => dltgamma_value
   end type dltgamma_result

   public :: lgamma_p11, dltgamma_inc

contains

   pure elemental real(dp) function lgamma_p11(p) result(v)
      real(dp), intent(in) :: p
      real(dp), parameter :: d(0:10) = [ &
         2.48574089138753565546e-5_dp, 1.05142378581721974210e0_dp, &
        -3.45687097222016235469e0_dp, 4.51227709466894823700e0_dp, &
        -2.98285225323576655721e0_dp, 1.05639711577126713077e0_dp, &
        -1.95428773191645869583e-1_dp, 1.70970543404441224307e-2_dp, &
        -5.71926117404305781283e-4_dp, 4.63399473359905636708e-6_dp, &
        -2.71994908488607703910e-9_dp ]
      real(dp) :: sumv, z
      integer :: k
      if (p <= 0.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      z = p - 1.0_dp
      sumv = d(0)
      do k = 1, 10
         sumv = sumv + d(k)/(z+real(k,dp))
      end do
      v = log(1.860382734205265717_dp*sumv) - (z+0.5_dp) &
         + (z+0.5_dp)*log(z+11.400511_dp)
   end function lgamma_p11

   pure elemental real(dp) function dltgamma_value(self) result(v)
      class(dltgamma_result), intent(in) :: self
      if (self%rho == 0.0_dp) then
         v = 0.0_dp
      else if (self%sigma > log(huge(1.0_dp))) then
         v = sign(ieee_value(0.0_dp,ieee_positive_inf),self%rho)
      else
         v = self%rho*exp(self%sigma)
      end if
   end function dltgamma_value

   pure function dltgamma_inc(x,y,mu,p) result(res)
      real(dp), intent(in) :: x,y,mu,p
      type(dltgamma_result) :: res
      real(dp) :: xx, yy, a, plo, phi, val, scale, fx, fy
      integer :: n
      if (p <= 0.0_dp .or. mu == 0.0_dp .or. .not.ieee_is_finite(mu) &
          .or. x < 0.0_dp .or. y < x) then
         res%rho = ieee_value(0.0_dp,ieee_quiet_nan)
         res%sigma = 0.0_dp
         res%method = -1
         return
      end if
      if (x == y) then
         res%rho=0.0_dp
         res%sigma=0.0_dp
         res%method=0
         return
      end if

      if (mu > 0.0_dp) then
         xx=mu*x
         if (ieee_is_finite(y)) then
            yy=mu*y
            ! Select lower- or upper-tail subtraction according to location.
            if (yy <= p+1.0_dp) then
               plo=pgamma(xx,p,1.0_dp)
               phi=pgamma(yy,p,1.0_dp)
               val=phi-plo
               res%method=1
            else
               plo=1.0_dp-pgamma(xx,p,1.0_dp)
               phi=1.0_dp-pgamma(yy,p,1.0_dp)
               val=plo-phi
               res%method=2
            end if
         else
            val=1.0_dp-pgamma(xx,p,1.0_dp)
            res%method=2
         end if
         if (val <= 0.0_dp) then
            res%rho=max(0.0_dp,val)
            res%sigma=0.0_dp
         else
            res%rho=1.0_dp
            res%sigma=r_lgamma(p)-p*log(mu)+log(val)
         end if
      else
         ! For mu < 0 the Algorithm 1006 domain requires integer p and finite y.
         if (.not.ieee_is_finite(y) .or. abs(p-anint(p))>16.0_dp*epsilon(1.0_dp)) then
            res%rho=ieee_value(0.0_dp,ieee_quiet_nan)
            res%sigma=0.0_dp
            res%method=-2
            return
         end if
         n=nint(p)
         a=-mu
         xx=a*x
         yy=a*y
         fx=expint_pos_integer(xx,n)
         fy=expint_pos_integer(yy,n)
         val=(fy-fx)/(a**n)
         res%method=3
         if (val == 0.0_dp) then
            res%rho=0.0_dp
            res%sigma=0.0_dp
         else if (ieee_is_finite(val)) then
            res%rho=val
            res%sigma=0.0_dp
         else
            ! Recompute in log scale when the direct recurrence overflows.
            scale=max(xx,yy)
            fx=expint_pos_integer_scaled(xx,n,scale)
            fy=expint_pos_integer_scaled(yy,n,scale)
            res%rho=(fy-fx)/(a**n)
            res%sigma=scale
         end if
      end if
   end function dltgamma_inc

   pure real(dp) function expint_pos_integer(t,n) result(v)
      real(dp), intent(in) :: t
      integer, intent(in) :: n
      real(dp) :: poly
      integer :: j
      if (n <= 0) then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      poly=1.0_dp
      do j=2,n
         poly=t**(j-1)-real(j-1,dp)*poly
      end do
      v=exp(t)*poly
   end function expint_pos_integer

   pure real(dp) function expint_pos_integer_scaled(t,n,scale) result(v)
      real(dp), intent(in) :: t,scale
      integer, intent(in) :: n
      real(dp) :: poly
      integer :: j
      poly=1.0_dp
      do j=2,n
         poly=t**(j-1)-real(j-1,dp)*poly
      end do
      v=exp(t-scale)*poly
   end function expint_pos_integer_scaled

end module dpq_toms1006
