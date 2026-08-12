program test_continuous
   use ceoptim
   implicit none
   type(ce_control) :: ctl
   type(ce_continuous_control) :: cc
   type(ce_result) :: res

   allocate(cc%mean(2), cc%sd(2))
   cc%mean = [-3.0_dp, -3.0_dp]
   cc%sd = [10.0_dp, 10.0_dp]
   cc%smooth_mean = 1.0_dp
   cc%smooth_sd = 1.0_dp
   cc%sd_thr = 1.0e-4_dp

   ctl%n = 500
   ctl%rho = 0.1_dp
   ctl%iter_thr = 150
   ctl%no_improve_thr = 12
   ctl%maximize = .true.
   ctl%seed = 1234_i64

   call ce_optimize(peaks, res, ctl, continuous=cc)
   if (res%status /= 0) error stop res%message
   if (res%optimum < 7.9_dp) error stop 'continuous CE optimum too low'
   if (abs(res%continuous(1) + 0.0093_dp) > 0.20_dp) error stop 'continuous x1 inaccurate'
   if (abs(res%continuous(2) - 1.5814_dp) > 0.20_dp) error stop 'continuous x2 inaccurate'
   if (res%actual_nfe /= (res%niter + 1) * ctl%n) error stop 'actual_nfe mismatch'
   if (res%nfe /= res%niter * ctl%n) error stop 'package-compatible nfe mismatch'
   print *, 'test_continuous: PASS', res%optimum, res%continuous

contains

   function peaks(xc, xd) result(v)
      real(dp), intent(in) :: xc(:)
      integer, intent(in) :: xd(:)
      real(dp) :: v
      real(dp) :: x, y
      if (size(xd) /= 0) error stop 'unexpected discrete argument'
      x = xc(1)
      y = xc(2)
      v = 3.0_dp*(1.0_dp-x)**2*exp(-x*x-(y+1.0_dp)**2) &
          - 10.0_dp*(x/5.0_dp-x**3-y**5)*exp(-x*x-y*y) &
          - (1.0_dp/3.0_dp)*exp(-(x+1.0_dp)**2-y*y)
   end function peaks
end program test_continuous
