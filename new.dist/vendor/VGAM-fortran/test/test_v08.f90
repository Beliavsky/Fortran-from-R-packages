program test_v08
   use vgam
   implicit none
   call test_positive_and_zero_altered()
   call test_zero_altered_regression()
   call test_trivariate_normal()
   call test_fgm_exponential()
   call test_gaitd_nb_dispersion()
   print '(a)', 'test_v08: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_positive_and_zero_altered()
      real(dp) :: p0, s
      integer :: k
      call assert_true(abs(dposnorm_v(0.5_dp, 0.0_dp, 1.0_dp) - &
         2.0_dp*dnorm_v(0.5_dp, 0.0_dp, 1.0_dp)) < 1.0e-14_dp, 'positive normal density')
      call assert_true(abs(pposnorm_v(qposnorm_v(0.37_dp, -0.4_dp, 1.2_dp), -0.4_dp, 1.2_dp) - &
         0.37_dp) < 2.0e-12_dp, 'positive normal quantile inversion')
      s = 0.0_dp
      do k = 0, 80
         s = s + dzapois_v(k, 2.3_dp, 0.18_dp)
      end do
      call assert_true(abs(s - 1.0_dp) < 1.0e-12_dp, 'zero-altered Poisson normalization')
      call assert_true(abs(dzapois_v(0, 2.3_dp, 0.18_dp) - 0.18_dp) < 1.0e-14_dp, &
         'zero-altered Poisson observed zero')
      p0 = dnbinom_v(0, 1.7_dp, 1.7_dp/(1.7_dp + 2.8_dp))
      call assert_true(abs(dzanbinom_v(0, 2.8_dp, 1.7_dp, 0.12_dp) - 0.12_dp) < 1.0e-14_dp, &
         'zero-altered negative-binomial observed zero')
      call assert_true(abs(pzanbinom_v(1, 2.8_dp, 1.7_dp, p0) - &
         pnbinom_v(1, 1.7_dp, 1.7_dp/(1.7_dp + 2.8_dp))) < 2.0e-14_dp, &
         'zero-altered NB reduces to parent at parent zero mass')
      call assert_true(abs(sum([(dzabinom_v(k, 7, 0.35_dp, 0.21_dp), k=0,7)]) - 1.0_dp) < 2.0e-14_dp, &
         'zero-altered binomial normalization')
      call assert_true(abs(sum([(dzigeom_v(k, 0.42_dp, -0.20_dp), k=0,160)]) - 1.0_dp) < 1.0e-12_dp, &
         'zero-deflated geometric normalization')
   end subroutine test_positive_and_zero_altered


   subroutine test_zero_altered_regression()
      integer, parameter :: n = 900
      integer :: y(n), i, nseed
      integer, allocatable :: seed(:)
      real(dp) :: x(n, 1), lambda0, pzero0
      type(zero_altered_count_result_t) :: fit
      lambda0 = 2.35_dp; pzero0 = 0.24_dp; x = 1.0_dp
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 808080; call random_seed(put=seed)
      do i = 1, n
         y(i) = rzapois_v(lambda0, pzero0)
      end do
      call fit_zero_altered_poisson(y, x, x, fit, max_iter=250, tol=2.0e-6_dp)
      call assert_true(abs(fit%fitted_parent_parameter(1) - lambda0) < 0.20_dp, &
         'zero-altered Poisson mean recovery')
      call assert_true(abs(fit%fitted_zero_probability(1) - pzero0) < 0.06_dp, &
         'zero-altered Poisson zero-probability recovery')
   end subroutine test_zero_altered_regression

   subroutine test_trivariate_normal()
      integer, parameter :: n = 650
      real(dp) :: y(n, 3), mu(3), sd(3), x(3), r12, r13, r23, ld0
      type(trivariate_normal_result_t) :: fit
      integer :: i, nseed, stat
      integer, allocatable :: seed(:)
      mu = [0.4_dp, -0.6_dp, 0.8_dp]; sd = [1.1_dp, 0.75_dp, 1.35_dp]
      r12 = 0.45_dp; r13 = -0.25_dp; r23 = 0.30_dp
      ld0 = trivariate_normal_logpdf(mu(1), mu(2), mu(3), mu(1), mu(2), mu(3), &
         sd(1), sd(2), sd(3), r12, r13, r23)
      call assert_true(abs(exp(ld0) - trivariate_normal_pdf(mu(1), mu(2), mu(3), &
         mu(1), mu(2), mu(3), sd(1), sd(2), sd(3), r12, r13, r23)) < 1.0e-15_dp, &
         'trivariate normal pdf/logpdf identity')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 808081; call random_seed(put=seed)
      do i = 1, n
         call random_trivariate_normal(mu, sd, r12, r13, r23, x, stat)
         call assert_true(stat == 0, 'trivariate normal RNG status')
         y(i, :) = x
      end do
      call fit_trivariate_normal(y, fit, max_iter=300, tol=2.0e-6_dp)
      call assert_true(maxval(abs(fit%mean - mu)) < 0.14_dp, 'trivariate normal means')
      call assert_true(maxval(abs(fit%sd - sd)) < 0.15_dp, 'trivariate normal sds')
      call assert_true(abs(fit%rho12 - r12) < 0.10_dp, 'trivariate normal rho12')
      call assert_true(abs(fit%rho13 - r13) < 0.10_dp, 'trivariate normal rho13')
      call assert_true(abs(fit%rho23 - r23) < 0.10_dp, 'trivariate normal rho23')
   end subroutine test_trivariate_normal

   subroutine test_fgm_exponential()
      integer, parameter :: n = 700
      real(dp) :: y1(n), y2(n), x(n, 1), alpha
      type(bifgm_exponential_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      alpha = 0.55_dp; x = 1.0_dp
      call assert_true(abs(bifgm_exponential_cdf(1.0_dp, 1.0_dp, 0.0_dp) - &
         (1.0_dp - exp(-1.0_dp))**2) < 1.0e-14_dp, 'FGM exponential independence CDF')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 808082; call random_seed(put=seed)
      do i = 1, n
         call random_bifgm_exponential(alpha, y1(i), y2(i))
      end do
      call fit_bifgm_exponential(y1, y2, x, fit, max_iter=250, tol=2.0e-6_dp)
      call assert_true(abs(fit%fitted_alpha(1) - alpha) < 0.16_dp, 'FGM exponential alpha recovery')
   end subroutine test_fgm_exponential

   subroutine test_gaitd_nb_dispersion()
      integer, parameter :: n = 1300, kmax = 60
      integer :: y(n), i, k, nseed
      integer, allocatable :: seed(:)
      real(dp) :: xp(n, 1), xm(n, 1), xo(n, 1), xs(n, 1), u, cum
      type(gaitd_distribution_t) :: dist
      type(gaitd_mix_nb_dispersion_result_t) :: fit
      xp = 1.0_dp; xm = 1.0_dp; xo = 1.0_dp; xs = 1.0_dp
      call gaitd_mix_negative_binomial(2.6_dp, 1.4_dp, kmax, dist, a_mix=[0, 3, 6], i_mix=[1, 4, 7], &
         d_mix=[2, 5, 8], pobs_mix=0.10_dp, pstr_mix=0.07_dp, pdip_mix=0.012_dp, &
         mu_a=0.85_dp, size_a=4.5_dp, mu_i=3.4_dp, size_i=0.9_dp, mu_d=5.2_dp, size_d=6.0_dp)
      call assert_true(dist%status == 0, 'GAITD separate-dispersion source distribution')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 808083; call random_seed(put=seed)
      do i = 1, n
         call random_number(u); cum = 0.0_dp; y(i) = kmax
         do k = 0, kmax
            cum = cum + dist%probability(k)
            if (u <= cum) then; y(i) = k; exit; end if
         end do
      end do
      call fit_gaitd_mix_nb_dispersion_regression(y, xp, xm, xo, xs, fit, a_mix=[0, 3, 6], i_mix=[1, 4, 7], &
         d_mix=[2, 5, 8], max_iter=450, tol=4.0e-6_dp)
      call assert_true(allocated(fit%fitted_size), 'GAITD separate dispersion fit allocated')
      call assert_true(all(fit%fitted_size(1, :) > 0.0_dp), 'GAITD separate dispersions positive')
      call assert_true(abs(exp(fit%parent_coefficients(1)) - 2.6_dp) < 0.65_dp, &
         'GAITD separate-dispersion parent mean recovery')
      call assert_true(abs(fit%fitted_size(1, 1) - 1.4_dp) < 1.2_dp, &
         'GAITD parent dispersion recovery')
      call assert_true(maxval(fit%fitted_size(1, :)) - minval(fit%fitted_size(1, :)) > 0.3_dp, &
         'GAITD outer dispersions are independently estimated')
   end subroutine test_gaitd_nb_dispersion
end program test_v08
