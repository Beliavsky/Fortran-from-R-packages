program demo_imputets
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use imputets_api, only : dp, na_interpolation
   implicit none

   real(dp) :: x(8)
   real(dp), allocatable :: y(:)

   x = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
        ieee_value(0.0_dp, ieee_quiet_nan), 7.0_dp, 8.0_dp]
   y = na_interpolation(x)
   write (*, '(*(f8.3,1x))') y
end program demo_imputets
