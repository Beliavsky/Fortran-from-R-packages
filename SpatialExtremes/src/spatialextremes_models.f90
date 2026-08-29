module spatialextremes_models
   use spatialextremes_base, only: dp
   use r_compat, only: normal_cdf,pt
   implicit none
   private
   public :: extremal_coefficient_smith,extremal_coefficient_schlather
   public :: extremal_coefficient_schlather_ind,extremal_coefficient_extremalt
   public :: concurrence_from_extcoeff
contains
   pure elemental real(dp) function extremal_coefficient_smith(a) result(theta)
      real(dp),intent(in)::a
      theta=2.0_dp*normal_cdf(0.5_dp*a)
   end function
   pure elemental real(dp) function extremal_coefficient_schlather(rho) result(theta)
      real(dp),intent(in)::rho
      theta=1.0_dp+sqrt(max(0.0_dp,(1.0_dp-rho)/2.0_dp))
   end function
   pure elemental real(dp) function extremal_coefficient_schlather_ind(rho,alpha) result(theta)
      real(dp),intent(in)::rho,alpha
      theta=1.0_dp+alpha+(1.0_dp-alpha)*sqrt(max(0.0_dp,(1.0_dp-rho)/2.0_dp))
   end function
   pure elemental real(dp) function extremal_coefficient_extremalt(rho,nu) result(theta)
      real(dp),intent(in)::rho,nu
      theta=2.0_dp*pt(sqrt((nu+1.0_dp)*(1.0_dp-rho)/(1.0_dp+rho)),nu+1.0_dp)
   end function
   pure elemental real(dp) function concurrence_from_extcoeff(theta) result(p)
      real(dp),intent(in)::theta
      ! For simple max-stable pairs, extremal concurrence p = 2 - theta.
      p=max(0.0_dp,min(1.0_dp,2.0_dp-theta))
   end function
end module spatialextremes_models
