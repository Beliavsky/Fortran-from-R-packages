program test_tmb
   use tmb, only: dp, dnorm, pnorm, qnorm, dexp, pexp, qexp, pweibull, qweibull, dbinom_robust, &
                  kronecker_product, find_interval, sort_real, interpolate2d, mvnorm_nll, ar1_nll, matrix_exponential, &
        romberg_integrate, gradient_fd, hessian_fd, df_density, dt_density, dsn, dmultinom, dshasho, pshasho, qshasho, &
                  unstructured_corr, unstructured_corr_nll, ar1_mvn_nll, n01_nll
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(dp) :: a(2, 2), b(2, 1), k(4, 2), data(3, 3), sigma(2, 2), x2(2), e(2, 2)
   real(dp) :: g(2), h(2, 2), x(2), val, corr(3, 3), xm(2, 3)
   integer :: info
   integer :: failures
   failures = 0
   call check_close(dnorm(0.0_dp, 0.0_dp, 1.0_dp, .false.), 0.3989422804014327_dp, 1.0e-14_dp, 'dnorm', failures)
   call check_close(pnorm(0.0_dp, 0.0_dp, 1.0_dp), 0.5_dp, 1.0e-15_dp, 'pnorm', failures)
   call check_close(qnorm(0.975_dp, 0.0_dp, 1.0_dp), 1.959963986120195_dp, 3.0e-8_dp, 'qnorm', failures)
   call check_close(pexp(2.0_dp, 0.5_dp), 1.0_dp - exp(-1.0_dp), 1.0e-15_dp, 'pexp', failures)
   call check_close(qexp(1.0_dp - exp(-1.0_dp), 0.5_dp), 2.0_dp, 1.0e-14_dp, 'qexp', failures)
   call check_close(dexp(0.0_dp, 2.0_dp, .false.), 2.0_dp, 1.0e-15_dp, 'dexp', failures)
   call check_close(qweibull(pweibull(3.0_dp, 2.0_dp, 4.0_dp), 2.0_dp, 4.0_dp), 3.0_dp, 1.0e-14_dp, 'weibull', failures)
   call check_close(dbinom_robust(2.0_dp, 4.0_dp, 0.0_dp, .false.), 0.375_dp, 1.0e-14_dp, 'dbinom robust', failures)
call check_close(df_density(1.0_dp, 5.0_dp, 10.0_dp, .false.), 0.4954797834866388_dp, 2.0e-15_dp, 'F density', failures)
   call check_close(dt_density(0.0_dp, 5.0_dp, .false.), 0.3796066898224944_dp, 2.0e-15_dp, 't density', failures)
   call check_close(dsn(0.0_dp, 3.0_dp, .false.), dnorm(0.0_dp, 0.0_dp, 1.0_dp, .false.), 2.0e-15_dp, 'skew normal', failures)
   call check_close(dmultinom([1.0_dp, 1.0_dp], [0.25_dp, 0.75_dp], .false.), 0.375_dp, 2.0e-15_dp, 'multinomial', failures)
   val = pshasho(0.3_dp, 0.2_dp, 1.4_dp, -0.1_dp, 1.2_dp, .false.)
  call check_close(qshasho(val, 0.2_dp, 1.4_dp, -0.1_dp, 1.2_dp, .false.), 0.3_dp, 5.0e-8_dp, 'SHASH inverse', failures)
   if (dshasho(0.3_dp, 0.2_dp, 1.4_dp, -0.1_dp, 1.2_dp, .false.) <= 0.0_dp) failures = failures + 1
   a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2, 2])
   b(:, 1) = [5.0_dp, 6.0_dp]
   k = kronecker_product(a, b)
   call check_close(k(4, 2), 24.0_dp, 1.0e-15_dp, 'kronecker', failures)
   if (find_interval(1.5_dp, [0.0_dp, 1.0_dp, 2.0_dp]) /= 2) failures = failures + 1
   if (any(abs(sort_real([3.0_dp, 1.0_dp, 2.0_dp]) - [1.0_dp, 2.0_dp, 3.0_dp]) > 0.0_dp)) failures = failures + 1
   data = reshape([0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [3, 3])
   val = interpolate2d(data, [0.0_dp, 2.0_dp], [0.0_dp, 2.0_dp], 1.0_dp, 1.0_dp, 1.0_dp)
   call check_close(val, 2.0_dp, 1.0e-12_dp, 'interpolate2d', failures)
   sigma = reshape([1.0_dp, 0.3_dp, 0.3_dp, 2.0_dp], [2, 2])
   x2 = [0.0_dp, 0.0_dp]
   call check_close(mvnorm_nll(x2, sigma), log(2.0_dp * acos(-1.0_dp)) + 0.5_dp * log(1.91_dp), 1.0e-13_dp, 'mvnorm', failures)
   call check_close(ar1_nll([0.0_dp, 0.0_dp], 0.5_dp), &
                    log(2.0_dp * acos(-1.0_dp)) + 0.5_dp * log(0.75_dp), 1.0e-13_dp, 'ar1', failures)
   call check_close(n01_nll(0.0_dp), 0.5_dp * log(2.0_dp * acos(-1.0_dp)), 1.0e-15_dp, 'N01', failures)
   call unstructured_corr([0.2_dp, -0.3_dp, 0.4_dp], corr, info)
   if (info /= 0 .or. maxval(abs([(corr(info, info), info=1, 3)] - 1.0_dp)) > 2.0e-15_dp) failures = failures + 1
   if (unstructured_corr_nll([0.0_dp, 0.0_dp, 0.0_dp], [0.2_dp, -0.3_dp, 0.4_dp]) <= 0.0_dp) failures = failures + 1
   xm = 0.0_dp
   call check_close(ar1_mvn_nll(xm, 0.5_dp, sigma), &
                    3.0_dp * mvnorm_nll([0.0_dp, 0.0_dp], sigma) + 4.0_dp * log(sqrt(0.75_dp)), &
                    1.0e-13_dp, 'multivariate ar1', failures)
   a = reshape([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [2, 2])
   e = matrix_exponential(a * (0.5_dp * acos(-1.0_dp)))
   call check_close(e(2, 1), 1.0_dp, 1.0e-13_dp, 'matrix exponential', failures)
   call check_close(romberg_integrate(square_fun, 0.0_dp, 1.0_dp), 1.0_dp / 3.0_dp, 1.0e-13_dp, 'romberg', failures)
   x = [1.0_dp, -2.0_dp]
   call gradient_fd(quadratic, x, g, 1.0e-5_dp)
   call hessian_fd(quadratic, x, h, 1.0e-3_dp)
   if (maxval(abs(g - [2.0_dp, -8.0_dp])) > 1.0e-8_dp) failures = failures + 1
   if (maxval(abs(h - reshape([2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2]))) > 1.0e-7_dp) failures = failures + 1
   if (ieee_is_nan(val)) failures = failures + 1
   if (failures /= 0) error stop 'TMB tests failed'
   print '(a)', 'All TMB Fortran tests passed.'
contains
   pure real(dp) function square_fun(z) result(y)
      real(dp), intent(in) :: z !! Scalar at which z**2 is evaluated.
      y = z * z
   end function square_fun
   pure real(dp) function quadratic(z) result(y)
      real(dp), intent(in) :: z(:) !! Two-component point for a deterministic quadratic test objective.
      y = z(1)**2 + 2.0_dp * z(2)**2
   end function quadratic
   subroutine check_close(got, expected, tol, label, failures)
      real(dp), intent(in) :: got !! Computed scalar result.
      real(dp), intent(in) :: expected !! Reference scalar result.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error.
      character(*), intent(in) :: label !! Human-readable test label.
      integer, intent(inout) :: failures !! Running number of failed checks.
      if (abs(got - expected) > tol .or. ieee_is_nan(got)) then
         print '(a,2(1x,es24.16))', trim(label), got, expected
         failures = failures + 1
      end if
   end subroutine check_close
end program test_tmb
