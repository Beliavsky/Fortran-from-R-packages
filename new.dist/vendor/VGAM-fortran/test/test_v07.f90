program test_v07
   use vgam
   implicit none
   call test_information_tools()
   call test_censored_models()
   call test_mix_regression()
   call test_bivariate_normal_fit()
   call test_bivariate_logistic()
   call test_freund_fit()
   print '(a)', 'test_v07: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_information_tools()
      real(dp) :: scores(3, 2), full(3, 3), cmat(3, 2), freecov(2, 2)
      real(dp), allocatable :: info(:, :), cinfo(:, :), lifted(:, :)
      scores = reshape([1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 3.0_dp], [3, 2])
      call score_outer_information(scores, info)
      call assert_true(abs(info(1, 1) - sum(scores(:, 1)**2)) < 1.0e-14_dp, 'score outer information')
      full = reshape([4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 0.5_dp, &
         0.0_dp, 0.5_dp, 2.0_dp], [3, 3])
      cmat = reshape([1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, -1.0_dp], [3, 2])
      call constrained_information(full, cmat, cinfo)
      call assert_true(maxval(abs(cinfo - matmul(transpose(cmat), matmul(full, cmat)))) < 1.0e-14_dp, &
         'constrained information projection')
      freecov = reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2])
      call lift_constrained_covariance(freecov, cmat, lifted)
      call assert_true(maxval(abs(lifted - matmul(cmat, matmul(freecov, transpose(cmat))))) < 1.0e-14_dp, &
         'constraint covariance lift')
   end subroutine test_information_tools

   subroutine test_censored_models()
      integer, parameter :: n = 420
      real(dp) :: low(n), up(n), x(n, 1), xs(n, 1), z, mu0, sd0
      integer :: ct(n), i, nseed
      integer, allocatable :: seed(:)
      type(censored_regression_result_t) :: fit
      mu0 = 0.45_dp; sd0 = 1.15_dp; x = 1.0_dp; xs = 1.0_dp
      call assert_true(abs(exp(censored_poisson_logprob(3.0_dp, 0.0_dp, censor_right, 2.0_dp)) - &
         (1.0_dp - ppois_v(2, 2.0_dp))) < 1.0e-14_dp, 'censored Poisson right tail')
      call assert_true(abs(exp(censored_exponential_logprob(1.2_dp, 0.0_dp, censor_right, 0.8_dp)) - &
         exp(-0.8_dp*1.2_dp)) < 1.0e-14_dp, 'censored exponential right tail')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 170717; call random_seed(put=seed)
      do i = 1, n
         z = rnorm_v(mu0, sd0)
         up(i) = 0.0_dp
         if (z <= -0.75_dp) then
            low(i) = -0.75_dp; ct(i) = censor_left
         else if (z >= 1.35_dp) then
            low(i) = 1.35_dp; ct(i) = censor_right
         else
            low(i) = z; ct(i) = censor_exact
         end if
      end do
      call fit_censored_normal(low, up, ct, x, fit, x_sd=xs, max_iter=250, tol=2.0e-6_dp)
      call assert_true(allocated(fit%coefficients), 'censored normal fit allocated')
      call assert_true(abs(fit%coefficients(1) - mu0) < 0.18_dp, 'censored normal mean recovery')
      call assert_true(abs(exp(fit%scale_coefficients(1)) - sd0) < 0.18_dp, 'censored normal sd recovery')
   end subroutine test_censored_models

   subroutine test_mix_regression()
      integer, parameter :: n = 650
      integer :: y(n), i, k, nseed
      integer, allocatable :: seed(:)
      real(dp) :: xp(n, 1), xm(n, 1), xo(n, 1), u, cum, true_parent, true_mass, true_outer
      type(gaitd_distribution_t) :: dist
      type(gaitd_mix_regression_result_t) :: fit
      true_parent = 2.4_dp; true_mass = 0.17_dp; true_outer = 0.75_dp
      xp = 1.0_dp; xm = 1.0_dp; xo = 1.0_dp
      call gaitd_mix_poisson(true_parent, 45, dist, a_mix=[0, 3], pobs_mix=true_mass, lambda_a=true_outer)
      call assert_true(dist%status == 0, 'GAITD mix regression source distribution')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 707071; call random_seed(put=seed)
      do i = 1, n
         call random_number(u); cum = 0.0_dp; y(i) = 45
         do k = 0, 45
            cum = cum + dist%probability(k)
            if (u <= cum) then
               y(i) = k; exit
            end if
         end do
      end do
      call fit_gaitd_mix_poisson_regression(y, xp, xm, xo, fit, a_mix=[0, 3], &
         max_iter=320, tol=2.0e-6_dp)
      call assert_true(allocated(fit%parent_coefficients), 'GAITD outer-mix regression allocated')
      call assert_true(abs(exp(fit%parent_coefficients(1)) - true_parent) < 0.35_dp, &
         'GAITD parent mean recovery')
      call assert_true(abs(fit%fitted_mass(1, 1) - true_mass) < 0.08_dp, 'GAITD mix mass recovery')
      call assert_true(abs(fit%fitted_outer_mean(1, 1) - true_outer) < 0.35_dp, &
         'GAITD outer mean recovery')
   end subroutine test_mix_regression

   subroutine test_bivariate_normal_fit()
      integer, parameter :: n = 480
      real(dp) :: y1(n), y2(n), x(n, 1), mu1, mu2, sd1, sd2, rho
      type(bivariate_normal_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      mu1 = 0.4_dp; mu2 = -0.7_dp; sd1 = 1.2_dp; sd2 = 0.8_dp; rho = 0.52_dp; x = 1.0_dp
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 707072; call random_seed(put=seed)
      do i = 1, n
         call random_bivariate_normal(mu1, mu2, sd1, sd2, rho, y1(i), y2(i))
      end do
      call fit_bivariate_normal(y1, y2, x, x, x, x, x, fit, max_iter=250, tol=2.0e-6_dp)
      call assert_true(abs(fit%fitted_mean1(1) - mu1) < 0.16_dp, 'bivariate normal mean1 recovery')
      call assert_true(abs(fit%fitted_mean2(1) - mu2) < 0.14_dp, 'bivariate normal mean2 recovery')
      call assert_true(abs(fit%fitted_sd1(1) - sd1) < 0.14_dp, 'bivariate normal sd1 recovery')
      call assert_true(abs(fit%fitted_sd2(1) - sd2) < 0.12_dp, 'bivariate normal sd2 recovery')
      call assert_true(abs(fit%fitted_rho(1) - rho) < 0.10_dp, 'bivariate normal rho recovery')
   end subroutine test_bivariate_normal_fit

   subroutine test_bivariate_logistic()
      integer, parameter :: n = 420
      real(dp) :: y1(n), y2(n), x(n, 1), l1, l2, s1, s2
      type(bivariate_logistic_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      l1 = 0.3_dp; l2 = -0.4_dp; s1 = 0.9_dp; s2 = 1.2_dp; x = 1.0_dp
      call assert_true(abs(bivariate_logistic_cdf(l1, l2, l1, s1, l2, s2) - 1.0_dp/3.0_dp) < 1.0e-14_dp, &
         'bivariate logistic CDF at locations')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 707073; call random_seed(put=seed)
      do i = 1, n
         call random_bivariate_logistic(l1, s1, l2, s2, y1(i), y2(i))
      end do
      call fit_bivariate_logistic(y1, y2, x, fit, max_iter=260, tol=3.0e-6_dp)
      call assert_true(abs(fit%fitted_location1(1) - l1) < 0.22_dp, 'bivariate logistic location1 recovery')
      call assert_true(abs(fit%fitted_location2(1) - l2) < 0.25_dp, 'bivariate logistic location2 recovery')
      call assert_true(abs(fit%fitted_scale1(1) - s1) < 0.20_dp, 'bivariate logistic scale1 recovery')
      call assert_true(abs(fit%fitted_scale2(1) - s2) < 0.25_dp, 'bivariate logistic scale2 recovery')
   end subroutine test_bivariate_logistic

   subroutine test_freund_fit()
      integer, parameter :: n = 600
      real(dp) :: y1(n), y2(n), x(n, 1), a, ap, b, bp
      type(freund61_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      a = 1.1_dp; ap = 1.7_dp; b = 0.8_dp; bp = 1.35_dp; x = 1.0_dp
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 707074; call random_seed(put=seed)
      do i = 1, n
         call random_freund61(a, ap, b, bp, y1(i), y2(i))
      end do
      call fit_freund61(y1, y2, x, fit, max_iter=300, tol=2.0e-6_dp)
      call assert_true(allocated(fit%fitted_parameters), 'Freund fit allocated')
      call assert_true(abs(fit%fitted_parameters(1, 1) - a) < 0.20_dp, 'Freund a recovery')
      call assert_true(abs(fit%fitted_parameters(1, 2) - ap) < 0.35_dp, 'Freund ap recovery')
      call assert_true(abs(fit%fitted_parameters(1, 3) - b) < 0.17_dp, 'Freund b recovery')
      call assert_true(abs(fit%fitted_parameters(1, 4) - bp) < 0.30_dp, 'Freund bp recovery')
   end subroutine test_freund_fit

end program test_v07
