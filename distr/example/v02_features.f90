program v02_features
   use distr
   implicit none
   type(distribution_t) :: x, y, z
   real(dp) :: q

   x = gamma_dist(2.0_dp, 1.0_dp)
   y = lognormal_dist(0.0_dp, 0.5_dp)
   z = convolve_fft(x, y, grid_points=4096)
   print '(a,f12.8)', 'P(X+Y <= 3) ~= ', z%cdf(3.0_dp)

   x = normal_dist()
   q = x%quantile(log(1.0e-30_dp), lower_tail=.false., log_p=.true.)
   print '(a,es16.8)', 'upper 1e-30 normal quantile = ', q
   print '(a,es16.8)', 'log survival at quantile = ', x%logsf(q)

   y = power_transform(x, 2.0_dp)
   print '(a,f12.8)', 'P(N(0,1)^2 <= 2) = ', y%cdf(2.0_dp)
end program v02_features
