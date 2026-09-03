program compar_lynch_example
   use ape
   implicit none
   type(compar_lynch_result) :: fit
   real(dp) :: g(4, 4)
   real(dp) :: x(4, 2)
   integer :: info

   x(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp]
   x(:, 2) = [2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp]
   g = reshape([ &
      2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [4, 4])

   call compar_lynch_fit(x, g, fit, info)
   if (info /= 0) error stop 'compar.lynch example failed'
   print '(a,f12.6)', 'log likelihood: ', fit%log_likelihood
   print '(a,*(f10.5,1x))', 'trait means: ', fit%mean
end program compar_lynch_example
