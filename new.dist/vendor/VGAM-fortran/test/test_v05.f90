program test_v05
   use vgam
   implicit none
   call test_gaitd_mlm_distribution()
   call test_gaitd_mlm_regression()
   call test_copula_identities()
   call test_copula_regression()
   print '(a)', 'test_v05: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_gaitd_mlm_distribution()
      type(gaitd_distribution_t) :: dist
      real(dp) :: lambda, p0, p1, p3, delta, denom, tmp6
      lambda = 2.0_dp
      call gaitd_mlm_poisson(lambda, 30, dist, truncate=[5], &
         altered_points=[0], altered_probabilities=[0.12_dp], &
         inflated_points=[3], inflation_probabilities=[0.05_dp], &
         deflated_points=[1], deflation_probabilities=[0.02_dp])
      call assert_true(dist%status == 0, 'GAITD MLM distribution status')
      call assert_true(abs(sum(dist%pmf) - 1.0_dp) < 1.0e-12_dp, 'GAITD MLM PMF normalizes')
      call assert_true(abs(dist%probability(0) - 0.12_dp) < 1.0e-12_dp, 'GAITD altered mass')
      denom = 1.0_dp - dpois_v(5, lambda) - dpois_v(0, lambda)
      tmp6 = 1.0_dp - 0.12_dp - 0.05_dp + 0.02_dp
      delta = tmp6/denom
      p0 = 0.12_dp
      p1 = delta*dpois_v(1, lambda) - 0.02_dp
      p3 = delta*dpois_v(3, lambda) + 0.05_dp
      call assert_true(abs(dist%probability(0) - p0) < 1.0e-12_dp, 'GAITD altered oracle')
      call assert_true(abs(dist%probability(1) - p1) < 1.0e-12_dp, 'GAITD deflated oracle')
      call assert_true(abs(dist%probability(3) - p3) < 1.0e-12_dp, 'GAITD inflated oracle')
      call assert_true(dist%probability(5) == 0.0_dp, 'GAITD truncation')
      call gaitd_mlm_negative_binomial(2.2_dp, 1.7_dp, 35, dist, &
         inflated_points=[0], inflation_probabilities=[0.06_dp], &
         deflated_points=[2], deflation_probabilities=[0.015_dp])
      call assert_true(dist%status == 0, 'GAITD MLM NB distribution status')
      call assert_true(abs(sum(dist%pmf) - 1.0_dp) < 1.0e-12_dp, 'GAITD MLM NB normalizes')
   end subroutine test_gaitd_mlm_distribution

   subroutine test_gaitd_mlm_regression()
      integer, parameter :: n = 600
      integer :: y(n), i
      real(dp) :: xm(n, 1), xz(n, 1), target(3)
      type(gaitd_distribution_t) :: dist
      type(gaitd_mlm_regression_result_t) :: fit
      xm = 1.0_dp
      xz = 1.0_dp
      target = [0.10_dp, 0.035_dp, 0.015_dp]
      call gaitd_mlm_poisson(2.3_dp, 35, dist, altered_points=[0], &
         altered_probabilities=[target(1)], inflated_points=[4], &
         inflation_probabilities=[target(2)], deflated_points=[1], &
         deflation_probabilities=[target(3)])
      call assert_true(dist%status == 0, 'GAITD MLM synthetic distribution')
      do i = 1, n
         y(i) = dist%quantile((real(i, dp) - 0.5_dp)/real(n, dp))
      end do
      call fit_gaitd_mlm_poisson_regression(y, xm, xz, [0, 4, 1], &
         [gaitd_altered, gaitd_inflated, gaitd_deflated], fit, max_iter=350, tol=2.0e-6_dp)
      call assert_true(allocated(fit%fitted_mean), 'GAITD MLM regression allocated')
      call assert_true(fit%loglik > -huge(1.0_dp)/10.0_dp, 'GAITD MLM regression finite likelihood')
      call assert_true(abs(exp(fit%mean_coefficients(1)) - 2.3_dp) < 0.35_dp, &
         'GAITD MLM parent mean recovery')
      call assert_true(abs(fit%special_probabilities(1, 1) - target(1)) < 0.035_dp, &
         'GAITD altered regression recovery')
      call assert_true(abs(fit%special_probabilities(1, 2) - target(2)) < 0.025_dp, &
         'GAITD inflated regression recovery')
      call assert_true(abs(fit%special_probabilities(1, 3) - target(3)) < 0.012_dp, &
         'GAITD deflated regression recovery')
   end subroutine test_gaitd_mlm_regression

   subroutine test_copula_identities()
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp) :: exact
      exact = 0.25_dp + asin(0.5_dp)/(2.0_dp*pi)
      call assert_true(abs(clayton_copula_cdf(0.3_dp, 0.7_dp, 0.0_dp) - 0.21_dp) < 1.0e-14_dp, &
         'Clayton independence CDF')
      call assert_true(abs(frank_copula_cdf(0.3_dp, 0.7_dp, 1.0_dp) - 0.21_dp) < 1.0e-14_dp, &
         'Frank independence CDF')
      call assert_true(abs(fgm_copula_cdf(0.3_dp, 0.7_dp, 0.0_dp) - 0.21_dp) < 1.0e-14_dp, &
         'FGM independence CDF')
      call assert_true(abs(plackett_copula_cdf(0.3_dp, 0.7_dp, 1.0_dp) - 0.21_dp) < 1.0e-14_dp, &
         'Plackett independence CDF')
      call assert_true(abs(gaussian_copula_cdf(0.5_dp, 0.5_dp, 0.5_dp) - exact) < 2.0e-6_dp, &
         'Gaussian copula quadrant probability')
      call assert_true(abs(clayton_copula_pdf(0.4_dp, 0.6_dp, 0.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'Clayton independence density')
      call assert_true(abs(frank_copula_pdf(0.4_dp, 0.6_dp, 1.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'Frank independence density')
      call assert_true(abs(fgm_copula_pdf(0.4_dp, 0.6_dp, 0.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'FGM independence density')
      call assert_true(abs(plackett_copula_pdf(0.4_dp, 0.6_dp, 1.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'Plackett independence density')
      call assert_true(abs(gaussian_copula_pdf(0.4_dp, 0.6_dp, 0.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'Gaussian independence density')
   end subroutine test_copula_identities

   subroutine test_copula_regression()
      integer, parameter :: n = 500
      real(dp) :: u(n), v(n), x(n, 1), apar, mean_ap
      integer :: i, nseed
      integer, allocatable :: seed(:)
      type(copula_regression_result_t) :: fit
      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed = 271828
      call random_seed(put=seed)
      apar = 0.62_dp
      x = 1.0_dp
      do i = 1, n
         call random_fgm_copula(apar, u(i), v(i))
      end do
      call fit_copula_regression(u, v, x, copula_fgm, fit, max_iter=250, tol=2.0e-6_dp)
      call assert_true(allocated(fit%coefficients), 'copula regression allocated')
      call assert_true(fit%loglik > -huge(1.0_dp)/10.0_dp, 'copula regression finite likelihood')
      mean_ap = sum(fit%fitted_parameter)/real(n, dp)
      call assert_true(abs(mean_ap - apar) < 0.18_dp, 'FGM parameter recovery')
   end subroutine test_copula_regression

end program test_v05
