program test_volatility_methods
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use iso_fortran_env, only : int64
   use quarks
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: x(8), y(160)
   type(risk_result) :: vw, boot1, boot2, garch
   type(rng_state) :: rng1, rng2
   integer :: i

   x = [0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.005_dp, -0.01_dp, &
      0.02_dp, -0.025_dp]
   vw = vwhs(x, 0.75_dp, volatility_ewma, 0.94_dp)
   call assert_close(vw%var, 0.0211009520802234_dp, tol, 'vwhs VaR')
   call assert_close(vw%es, 0.0274688952172985_dp, tol, 'vwhs ES')

   call seed_rng(rng1, 1234567_int64)
   call seed_rng(rng2, 1234567_int64)
   boot1 = fhs(x, 0.75_dp, volatility_ewma, 0.94_dp, 5000, rng1)
   boot2 = fhs(x, 0.75_dp, volatility_ewma, 0.94_dp, 5000, rng2)
   call assert_close(boot1%var, boot2%var, 0.0_dp, 'seeded FHS VaR')
   call assert_close(boot1%es, boot2%es, 0.0_dp, 'seeded FHS ES')
   if (boot1%es < boot1%var) error stop 'FHS ES below VaR'

   do i = 1, size(y)
      y(i) = 0.008_dp * sin(0.17_dp * real(i, dp)) + &
         0.004_dp * cos(0.071_dp * real(i, dp))
   end do
   garch = vwhs(y, 0.95_dp, volatility_garch, 0.94_dp, 300)
   if (.not. ieee_is_finite(garch%var) .or. .not. ieee_is_finite(garch%es)) then
      error stop 'GARCH volatility path returned nonfinite values'
   end if
   if (garch%volatility_model /= volatility_garch) error stop 'wrong model tag'
   print *, 'test_volatility_methods: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_volatility_methods
