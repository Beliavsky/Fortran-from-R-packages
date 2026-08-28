module lavaan_robust_difference
   use lavaan_kinds, only : dp
   implicit none
   private

   type, public :: robust_difference_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: statistic_unscaled = 0.0_dp
      real(dp) :: scaling_factor = 1.0_dp
      real(dp) :: shift = 0.0_dp
      real(dp) :: df = 0.0_dp
      integer :: status = 0
   end type robust_difference_result

   public :: satorra_bentler_difference_2001
   public :: satorra_bentler_difference_2010
   public :: scaled_shifted_difference

contains

   subroutine satorra_bentler_difference_2001(t_restricted, df_restricted, c_restricted, &
                                               t_unrestricted, df_unrestricted, c_unrestricted, result)
      real(dp), intent(in) :: t_restricted, df_restricted, c_restricted
      real(dp), intent(in) :: t_unrestricted, df_unrestricted, c_unrestricted
      type(robust_difference_result), intent(out) :: result
      real(dp) :: ddf, cd

      ddf = df_restricted - df_unrestricted
      result%df = ddf
      result%statistic_unscaled = t_restricted - t_unrestricted
      if (ddf < 0.0_dp) then
         result%status = -1
         return
      else if (abs(ddf) <= epsilon(1.0_dp)) then
         result%scaling_factor = huge(1.0_dp)
         result%statistic = result%statistic_unscaled
         result%status = 1
         return
      end if
      cd = (df_restricted * c_restricted - df_unrestricted * c_unrestricted) / ddf
      if (cd <= tiny(1.0_dp)) then
         result%status = 2
         result%scaling_factor = cd
         result%statistic = huge(1.0_dp)
         return
      end if
      result%scaling_factor = cd
      result%statistic = result%statistic_unscaled / cd
      result%status = 0
   end subroutine satorra_bentler_difference_2001

   subroutine satorra_bentler_difference_2010(t_restricted, df_restricted, c_restricted, &
                                               t_unrestricted, df_unrestricted, c_cross, result)
      real(dp), intent(in) :: t_restricted, df_restricted, c_restricted
      real(dp), intent(in) :: t_unrestricted, df_unrestricted, c_cross
      type(robust_difference_result), intent(out) :: result
      real(dp) :: ddf, cd

      ddf = df_restricted - df_unrestricted
      result%df = ddf
      result%statistic_unscaled = t_restricted - t_unrestricted
      if (ddf < 0.0_dp) then
         result%status = -1
         return
      else if (abs(ddf) <= epsilon(1.0_dp)) then
         result%scaling_factor = huge(1.0_dp)
         result%statistic = result%statistic_unscaled
         result%status = 1
         return
      end if
      ! Satorra & Bentler (2010): c_cross is the scale factor of M10,
      ! i.e. the less restricted model evaluated at the restricted estimates.
      cd = (df_restricted * c_restricted - df_unrestricted * c_cross) / ddf
      if (cd <= tiny(1.0_dp)) then
         result%status = 2
         result%scaling_factor = cd
         result%statistic = huge(1.0_dp)
         return
      end if
      result%scaling_factor = cd
      result%statistic = result%statistic_unscaled / cd
      result%status = 0
   end subroutine satorra_bentler_difference_2010

   subroutine scaled_shifted_difference(t_restricted, df_restricted, a_restricted, b_restricted, &
                                        t_unrestricted, df_unrestricted, a_unrestricted, b_unrestricted, result)
      real(dp), intent(in) :: t_restricted, df_restricted, a_restricted, b_restricted
      real(dp), intent(in) :: t_unrestricted, df_unrestricted, a_unrestricted, b_unrestricted
      type(robust_difference_result), intent(out) :: result
      real(dp) :: ddf, ad, bd

      ddf = df_restricted - df_unrestricted
      result%df = ddf
      result%statistic_unscaled = t_restricted - t_unrestricted
      if (ddf <= 0.0_dp) then
         result%status = -1
         return
      end if
      ! If T* = a*T + b has E(T*)=df for each model, the difference
      ! transformation is obtained by differencing the corresponding
      ! first two moment corrections. This is useful for lavaan's T3/MV
      ! bookkeeping even when an estimator-specific wrapper supplies a,b.
      ad = (df_restricted * a_restricted - df_unrestricted * a_unrestricted) / ddf
      bd = b_restricted - b_unrestricted
      if (abs(ad) <= tiny(1.0_dp)) then
         result%status = 2
         return
      end if
      result%scaling_factor = 1.0_dp / ad
      result%shift = bd
      result%statistic = ad * result%statistic_unscaled + bd
      result%status = 0
   end subroutine scaled_shifted_difference

end module lavaan_robust_difference
