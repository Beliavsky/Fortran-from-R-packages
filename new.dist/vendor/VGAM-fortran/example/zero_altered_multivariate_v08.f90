program zero_altered_multivariate_v08
   use vgam
   implicit none
   integer, parameter :: n = 500
   integer :: y(n), i, nseed, stat
   integer, allocatable :: seed(:)
   real(dp) :: xmat(n, 1), tri(n, 3), draw(3)
   real(dp), parameter :: lambda0 = 2.4_dp, pzero0 = 0.22_dp
   real(dp), parameter :: means(3) = [0.2_dp, -0.5_dp, 0.7_dp]
   real(dp), parameter :: sds(3) = [1.0_dp, 0.8_dp, 1.2_dp]
   type(zero_altered_count_result_t) :: zfit
   type(trivariate_normal_result_t) :: tfit

   call random_seed(size=nseed)
   allocate(seed(nseed)); seed = 808008; call random_seed(put=seed)
   xmat = 1.0_dp
   do i = 1, n
      y(i) = rzapois_v(lambda0, pzero0)
      call random_trivariate_normal(means, sds, 0.40_dp, -0.20_dp, 0.25_dp, draw, stat)
      if (stat /= 0) error stop 'trivariate RNG failed'
      tri(i, :) = draw
   end do

   call fit_zero_altered_poisson(y, xmat, xmat, zfit)
   call fit_trivariate_normal(tri, tfit)

   print '(a,f10.4)', 'Zero-altered Poisson lambda: ', zfit%fitted_parent_parameter(1)
   print '(a,f10.4)', 'Observed zero probability:    ', zfit%fitted_zero_probability(1)
   print '(a,3f10.4)', 'Trivariate normal means:     ', tfit%mean
   print '(a,3f10.4)', 'Trivariate normal correlations: ', tfit%rho12, tfit%rho13, tfit%rho23
end program zero_altered_multivariate_v08
