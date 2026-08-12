program peaks_example
   use ceoptim
   implicit none
   type(ce_control) :: ctl
   type(ce_continuous_control) :: continuous
   type(ce_result) :: res

   continuous%mean = [-3.0_dp, -3.0_dp]
   continuous%sd = [10.0_dp, 10.0_dp]
   ctl%n = 300
   ctl%rho = 0.1_dp
   ctl%maximize = .true.
   ctl%seed = 1234_i64

   call ce_optimize(peaks, res, ctl, continuous=continuous)
   if (res%status /= 0) error stop res%message
   print '(a,f12.6)', 'maximum = ', res%optimum
   print '(a,2f12.6)', 'at x = ', res%continuous

contains
   function peaks(xc, xd) result(v)
      real(dp), intent(in) :: xc(:)
      integer, intent(in) :: xd(:)
      real(dp) :: v, x, y
      x = xc(1)
      y = xc(2)
      v = 3.0_dp*(1.0_dp-x)**2*exp(-x*x-(y+1.0_dp)**2) &
          - 10.0_dp*(x/5.0_dp-x**3-y**5)*exp(-x*x-y*y) &
          - (1.0_dp/3.0_dp)*exp(-(x+1.0_dp)**2-y*y)
      if (size(xd) /= 0) error stop 'unexpected discrete variables'
   end function peaks
end program peaks_example
