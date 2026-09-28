program car_basic_example
   use car, only: bc_power, dp, durbin_watson
   implicit none

   real(dp) :: y(5), transformed(5)
   real(dp), allocatable :: dw(:)
   integer :: status

   y = [1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, 5.0_dp]
   call bc_power(y, 0.0_dp, transformed, status)
   if (status /= 0) error stop 'bc_power failed'
   write(*, '(a,5(1x,f9.5))') 'log-transformed:', transformed

   call durbin_watson([1.0_dp, -0.4_dp, 0.2_dp, -0.1_dp, 0.05_dp], 1, dw, status=status)
   if (status /= 0) error stop 'durbin_watson failed'
   write(*, '(a,f9.5)') 'Durbin-Watson: ', dw(1)
end program car_basic_example
