program example_gamm4
   use gamm4
   implicit none
   integer, parameter :: n = 10
   type(gamm4_smooth_t) :: smooths(1)
   type(gamm4_result_t) :: fit
   real(dp) :: x(n), y(n), fixed(n, 1)
   integer :: i, status

   do i = 1, n
      x(i) = real(i - 1, dp) / real(n - 1, dp)
      y(i) = 1.0_dp + 1.5_dp * x(i) + 0.4_dp * x(i) * x(i) + 0.02_dp * sin(real(i, dp))
   end do
   fixed(:, 1) = 1.0_dp
   allocate(smooths(1)%basis(n, 2), smooths(1)%spec%penalties(2, 2, 1))
   smooths(1)%basis(:, 1) = x
   smooths(1)%basis(:, 2) = x * x
   smooths(1)%spec%penalties = 0.0_dp
   smooths(1)%spec%penalties(2, 2, 1) = 1.0_dp
   smooths(1)%label = 's(x)'
   call gamm4_fit(y, fixed, smooths, fit, status, reml=.false.)
   if (status /= 0) error stop 'gamm4 example fit failed'
   write(*, '(a,*(1x,f10.5))') 'coefficients:', fit%coefficients
   write(*, '(a,es13.5)') 'smoothing parameter:', fit%sp(1)
   write(*, '(a,f10.5)') 'residual scale:', fit%scale
end program example_gamm4
