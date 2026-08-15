program test_v06
   use vgam
   implicit none
   call test_student_t_functions()
   call test_amh_copula()
   call test_gaitd_mix()
   call test_bivariate_student_t_fit()
   call test_student_t_copula_fit()
   call test_inflated_families()
   print '(a)', 'test_v06: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_student_t_functions()
      real(dp) :: ld, oracle
      call assert_true(abs(student_t_pdf(0.0_dp, 1.0_dp) - 1.0_dp/pi) < 2.0e-14_dp, &
         'Student-t Cauchy density')
      call assert_true(abs(student_t_cdf(1.0_dp, 1.0_dp) - 0.75_dp) < 2.0e-13_dp, &
         'Student-t Cauchy CDF')
      call assert_true(abs(student_t_quantile(0.75_dp, 1.0_dp) - 1.0_dp) < 2.0e-11_dp, &
         'Student-t Cauchy quantile')
      oracle = -(7.0_dp/2.0_dp + 1.0_dp)*log(1.0_dp + &
         (0.4_dp**2 + 0.8_dp**2 - 2.0_dp*0.35_dp*0.4_dp*0.8_dp)/(7.0_dp*(1.0_dp - 0.35_dp**2))) &
         - log(2.0_dp*pi) - 0.5_dp*log(1.0_dp - 0.35_dp**2)
      ld = bivariate_student_t_logpdf(0.4_dp, 0.8_dp, 7.0_dp, 0.35_dp)
      call assert_true(abs(ld - oracle) < 2.0e-14_dp, 'bivariate Student-t upstream density formula')
   end subroutine test_student_t_functions

   subroutine test_amh_copula()
      real(dp) :: u, v
      call assert_true(abs(amh_copula_cdf(0.3_dp, 0.7_dp, 0.0_dp) - 0.21_dp) < 1.0e-14_dp, &
         'AMH independence CDF')
      call assert_true(abs(amh_copula_pdf(0.3_dp, 0.7_dp, 0.0_dp) - 1.0_dp) < 1.0e-14_dp, &
         'AMH independence PDF')
      call random_amh_copula(0.6_dp, u, v)
      call assert_true(u > 0.0_dp .and. u < 1.0_dp .and. v > 0.0_dp .and. v < 1.0_dp, &
         'AMH RNG support')
      call fit_amh_recovery()
   end subroutine test_amh_copula

   subroutine fit_amh_recovery()
      integer, parameter :: n = 420
      real(dp) :: u(n), v(n), x(n, 1), apar
      type(copula_regression_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 818181; call random_seed(put=seed)
      apar = 0.55_dp; x = 1.0_dp
      do i = 1, n
         call random_amh_copula(apar, u(i), v(i))
      end do
      call fit_copula_regression(u, v, x, copula_amh, fit, max_iter=220, tol=2.0e-6_dp)
      call assert_true(allocated(fit%fitted_parameter), 'AMH copula regression allocated')
      call assert_true(abs(fit%fitted_parameter(1) - apar) < 0.20_dp, 'AMH copula parameter recovery')
   end subroutine fit_amh_recovery

   subroutine test_gaitd_mix()
      type(gaitd_distribution_t) :: dist
      real(dp) :: lp, la, li, ld, pobs, pstr, pdip, cdfmax, suma, sumt, delta
      real(dp) :: wa0, wa2, wi4, wi5, wd1, wd3, oracle
      integer :: k
      lp = 2.0_dp; la = 0.7_dp; li = 4.5_dp; ld = 1.2_dp
      pobs = 0.10_dp; pstr = 0.04_dp; pdip = 0.02_dp
      call gaitd_mix_poisson(lp, 40, dist, a_mix=[0, 2], i_mix=[4, 5], d_mix=[1, 3], &
         pobs_mix=pobs, pstr_mix=pstr, pdip_mix=pdip, lambda_a=la, lambda_i=li, lambda_d=ld, &
         truncate=[6])
      call assert_true(dist%status == 0, 'GAITD mix Poisson status')
      call assert_true(abs(sum(dist%pmf) - 1.0_dp) < 2.0e-12_dp, 'GAITD mix normalizes')
      cdfmax = 0.0_dp
      do k = 0, 40
         cdfmax = cdfmax + dpois_v(k, lp)
      end do
      sumt = dpois_v(6, lp)
      suma = dpois_v(0, lp) + dpois_v(2, lp)
      delta = (1.0_dp - pobs - pstr + pdip)/(cdfmax - sumt - suma)
      wa0 = dpois_v(0, la)/(dpois_v(0, la) + dpois_v(2, la))
      wa2 = 1.0_dp - wa0
      wi4 = dpois_v(4, li)/(dpois_v(4, li) + dpois_v(5, li))
      wi5 = 1.0_dp - wi4
      wd1 = dpois_v(1, ld)/(dpois_v(1, ld) + dpois_v(3, ld))
      wd3 = 1.0_dp - wd1
      call assert_true(abs(dist%probability(0) - pobs*wa0) < 2.0e-13_dp, 'GAITD a.mix weight 0')
      call assert_true(abs(dist%probability(2) - pobs*wa2) < 2.0e-13_dp, 'GAITD a.mix weight 2')
      oracle = delta*dpois_v(4, lp) + pstr*wi4
      call assert_true(abs(dist%probability(4) - oracle) < 2.0e-13_dp, 'GAITD i.mix weight 4')
      oracle = delta*dpois_v(5, lp) + pstr*wi5
      call assert_true(abs(dist%probability(5) - oracle) < 2.0e-13_dp, 'GAITD i.mix weight 5')
      oracle = delta*dpois_v(1, lp) - pdip*wd1
      call assert_true(abs(dist%probability(1) - oracle) < 2.0e-13_dp, 'GAITD d.mix weight 1')
      oracle = delta*dpois_v(3, lp) - pdip*wd3
      call assert_true(abs(dist%probability(3) - oracle) < 2.0e-13_dp, 'GAITD d.mix weight 3')
      call assert_true(dist%probability(6) == 0.0_dp, 'GAITD mix truncation')
      call gaitd_mix_negative_binomial(2.5_dp, 1.8_dp, 50, dist, a_mix=[0, 2], &
         pobs_mix=0.08_dp, mu_a=0.9_dp, size_a=1.4_dp)
      call assert_true(dist%status == 0, 'GAITD mix NB status')
      call assert_true(abs(sum(dist%pmf) - 1.0_dp) < 2.0e-12_dp, 'GAITD mix NB normalizes')
   end subroutine test_gaitd_mix

   subroutine test_bivariate_student_t_fit()
      integer, parameter :: n = 500
      real(dp) :: y1(n), y2(n), xd(n, 1), xr(n, 1), df0, rho0
      type(bivariate_student_t_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 606060; call random_seed(put=seed)
      df0 = 7.0_dp; rho0 = 0.55_dp; xd = 1.0_dp; xr = 1.0_dp
      do i = 1, n
         call random_bivariate_student_t(df0, rho0, y1(i), y2(i))
      end do
      call fit_bivariate_student_t(y1, y2, xd, xr, fit, max_iter=250, tol=2.0e-6_dp)
      call assert_true(allocated(fit%fitted_df), 'bivariate Student-t fit allocated')
      call assert_true(fit%loglik > -huge(1.0_dp)/10.0_dp, 'bivariate Student-t finite likelihood')
      call assert_true(abs(fit%fitted_rho(1) - rho0) < 0.13_dp, 'bivariate Student-t rho recovery')
      call assert_true(fit%fitted_df(1) > 2.0_dp .and. fit%fitted_df(1) < 20.0_dp, &
         'bivariate Student-t df recovery range')
   end subroutine test_bivariate_student_t_fit

   subroutine test_student_t_copula_fit()
      integer, parameter :: n = 260
      real(dp) :: u(n), v(n), x(n, 1), df0, rho0
      type(student_t_copula_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 707070; call random_seed(put=seed)
      df0 = 6.0_dp; rho0 = 0.45_dp; x = 1.0_dp
      do i = 1, n
         call random_student_t_copula(df0, rho0, u(i), v(i))
      end do
      call fit_student_t_copula(u, v, x, fit, estimate_df=.false., initial_df=df0, &
         max_iter=180, tol=2.0e-6_dp)
      call assert_true(allocated(fit%fitted_rho), 'Student-t copula fit allocated')
      call assert_true(fit%loglik > -huge(1.0_dp)/10.0_dp, 'Student-t copula finite likelihood')
      call assert_true(abs(fit%fitted_rho(1) - rho0) < 0.15_dp, 'Student-t copula rho recovery')
      call assert_true(abs(fit%df - df0) < 1.0e-14_dp, 'Student-t copula fixed df')
   end subroutine test_student_t_copula_fit

   subroutine test_inflated_families()
      real(dp) :: p0, p1, a, b, x, oracle, psum
      integer :: k, nsize
      p0 = 0.12_dp; p1 = 0.08_dp; a = 2.3_dp; b = 4.1_dp
      call assert_true(abs(dzoabeta(0.0_dp, a, b, p0, p1) - p0) < 1.0e-14_dp, &
         'zero-one altered beta mass zero')
      call assert_true(abs(dzoabeta(1.0_dp, a, b, p0, p1) - p1) < 1.0e-14_dp, &
         'zero-one altered beta mass one')
      x = 0.37_dp
      oracle = (1.0_dp - p0 - p1)*dbeta_v(x, a, b)
      call assert_true(abs(dzoabeta(x, a, b, p0, p1) - oracle) < 2.0e-14_dp, &
         'zero-one altered beta interior density')
      call assert_true(abs(pzoabeta(qzoabeta(0.63_dp, a, b, p0, p1), a, b, p0, p1) - 0.63_dp) < 2.0e-10_dp, &
         'zero-one altered beta quantile inversion')
      nsize = 8; psum = 0.0_dp
      do k = 0, nsize
         psum = psum + dzoibetabinom_ab(k, nsize, a, b, 0.10_dp, 0.06_dp)
      end do
      call assert_true(abs(psum - 1.0_dp) < 2.0e-13_dp, 'zero-one inflated beta-binomial normalizes')
      oracle = (1.0_dp - 0.10_dp - 0.06_dp)*dbetabinom_ab(0, nsize, a, b) + 0.10_dp
      call assert_true(abs(dzoibetabinom_ab(0, nsize, a, b, 0.10_dp, 0.06_dp) - oracle) < 2.0e-14_dp, &
         'zero-one inflated beta-binomial zero mass')
   end subroutine test_inflated_families

end program test_v06
