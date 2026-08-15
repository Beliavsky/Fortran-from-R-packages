program censored_bivariate_v07
   use vgam
   implicit none
   integer, parameter :: n = 300
   real(dp) :: low(n), up(n), x(n, 1), xs(n, 1), z
   real(dp) :: y1(n), y2(n), xb(n, 1)
   integer :: ct(n), i, nseed
   integer, allocatable :: seed(:)
   type(censored_regression_result_t) :: cfit
   type(bivariate_normal_result_t) :: bfit

   call random_seed(size=nseed)
   allocate(seed(nseed)); seed = 7070707
   call random_seed(put=seed)

   x = 1.0_dp; xs = 1.0_dp
   do i = 1, n
      z = rnorm_v(0.5_dp, 1.1_dp)
      up(i) = 0.0_dp
      if (z <= -0.75_dp) then
         low(i) = -0.75_dp; ct(i) = censor_left
      else if (z >= 1.35_dp) then
         low(i) = 1.35_dp; ct(i) = censor_right
      else
         low(i) = z; ct(i) = censor_exact
      end if
   end do
   call fit_censored_normal(low, up, ct, x, cfit, x_sd=xs)
   print '(a,f9.4)', 'Censored-normal fitted mean: ', cfit%coefficients(1)
   print '(a,f9.4)', 'Censored-normal fitted sd:   ', exp(cfit%scale_coefficients(1))

   xb = 1.0_dp
   do i = 1, n
      call random_bivariate_normal(0.3_dp, -0.4_dp, 1.2_dp, 0.8_dp, 0.5_dp, y1(i), y2(i))
   end do
   call fit_bivariate_normal(y1, y2, xb, xb, xb, xb, xb, bfit)
   print '(a,2f9.4)', 'Bivariate-normal means:      ', bfit%fitted_mean1(1), bfit%fitted_mean2(1)
   print '(a,2f9.4)', 'Bivariate-normal sds:        ', bfit%fitted_sd1(1), bfit%fitted_sd2(1)
   print '(a,f9.4)', 'Bivariate-normal rho:        ', bfit%fitted_rho(1)
end program censored_bivariate_v07
