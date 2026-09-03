program simple
   use tmb, only: dp, dnorm, mvnorm_nll
   implicit none
   real(dp) :: sigma(2, 2)
   sigma = reshape([1.0_dp, 0.25_dp, 0.25_dp, 1.0_dp], [2, 2])
   print '(a,f12.8)', 'standard normal density at zero = ', dnorm(0.0_dp, 0.0_dp, 1.0_dp, .false.)
   print '(a,f12.8)', 'bivariate Gaussian NLL at [0,0] = ', mvnorm_nll([0.0_dp, 0.0_dp], sigma)
end program simple
