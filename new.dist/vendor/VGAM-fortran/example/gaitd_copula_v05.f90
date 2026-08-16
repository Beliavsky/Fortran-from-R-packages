program gaitd_copula_v05
   use vgam
   implicit none
   integer, parameter :: n = 250
   integer :: i, nseed
   integer, allocatable :: seed(:)
   real(dp) :: u(n), v(n), x(n, 1)
   type(gaitd_distribution_t) :: gd
   type(copula_regression_result_t) :: cfit

   call gaitd_mlm_poisson(2.0_dp, 30, gd, altered_points=[0], &
      altered_probabilities=[0.12_dp], inflated_points=[3], &
      inflation_probabilities=[0.04_dp], deflated_points=[1], &
      deflation_probabilities=[0.015_dp])
   print '(a,f9.5)', 'GAITD P(Y=0): ', gd%probability(0)
   print '(a,f9.5)', 'GAITD P(Y=1): ', gd%probability(1)
   print '(a,f9.5)', 'GAITD mean:   ', gd%mean

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = 314159
   call random_seed(put=seed)
   x = 1.0_dp
   do i = 1, n
      call random_fgm_copula(0.55_dp, u(i), v(i))
   end do
   call fit_copula_regression(u, v, x, copula_fgm, cfit, max_iter=250, tol=2.0e-6_dp)
   print '(a,f9.5)', 'FGM true parameter:   ', 0.55_dp
   print '(a,f9.5)', 'FGM fitted parameter: ', cfit%fitted_parameter(1)
end program gaitd_copula_v05
