program test_rollcast
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use quarks
   implicit none
   real(dp), parameter :: tol = 2.0e-13_dp
   real(dp) :: x(30)
   type(rollcast_result) :: rolled, smoothed
   type(risk_result) :: direct
   integer :: i, first, last

   do i = 1, size(x)
      x(i) = 0.01_dp * sin(0.41_dp * real(i, dp)) + &
         0.003_dp * cos(0.13_dp * real(i, dp))
   end do
   rolled = rollcast(x, p=0.80_dp, method=method_plain, nout=6, nwin=12)
   if (size(rolled%var) /= 6 .or. size(rolled%xout) /= 6) error stop 'wrong sizes'
   do i = 1, 6
      last = 24 + i - 1
      first = last - 11
      direct = hs(x(first:last), 0.80_dp, method_plain)
      call assert_close(rolled%var(i), direct%var, tol, 'rolling VaR')
      call assert_close(rolled%es(i), direct%es, tol, 'rolling ES')
   end do

   smoothed = rollcast(x, p=0.80_dp, method=method_age, lambda=0.97_dp, &
      nout=4, nwin=15, smoothing=smooth_lpr, smoothing_bandwidth=0.20_dp)
   if (any(.not. ieee_is_finite(smoothed%var)) .or. &
       any(.not. ieee_is_finite(smoothed%es))) then
      error stop 'smoothed rolling forecasts are nonfinite'
   end if
   print *, 'test_rollcast: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_rollcast
