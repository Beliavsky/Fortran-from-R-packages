program student_mix_v06
   use vgam
   implicit none
   integer, parameter :: n = 400
   integer :: i, nseed
   integer, allocatable :: seed(:)
   real(dp) :: y1(n), y2(n), xd(n, 1), xr(n, 1)
   type(gaitd_distribution_t) :: gd
   type(bivariate_student_t_result_t) :: fit

   call gaitd_mix_poisson(2.0_dp, 40, gd, a_mix=[0, 2], i_mix=[4, 5], &
      pobs_mix=0.10_dp, pstr_mix=0.04_dp, lambda_a=0.7_dp, lambda_i=4.5_dp, truncate=[6])
   print '(a,f10.6)', 'GAITD mix P(Y=0): ', gd%probability(0)
   print '(a,f10.6)', 'GAITD mix P(Y=4): ', gd%probability(4)
   print '(a,f10.6)', 'GAITD mix mean:   ', gd%mean

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = 6062026
   call random_seed(put=seed)
   xd = 1.0_dp
   xr = 1.0_dp
   do i = 1, n
      call random_bivariate_student_t(7.0_dp, 0.50_dp, y1(i), y2(i))
   end do
   call fit_bivariate_student_t(y1, y2, xd, xr, fit, max_iter=250, tol=2.0e-6_dp)
   print '(a,f10.5)', 'Bivariate-t fitted df:  ', fit%fitted_df(1)
   print '(a,f10.5)', 'Bivariate-t fitted rho: ', fit%fitted_rho(1)
end program student_mix_v06
