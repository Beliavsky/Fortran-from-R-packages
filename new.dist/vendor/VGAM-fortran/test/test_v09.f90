program test_v09
   use vgam
   implicit none
   call test_dirichlet()
   call test_tobit_and_folded()
   call test_positive_nb()
   call test_zoa_beta()
   print '(a)', 'test_v09: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_dirichlet()
      integer, parameter :: n = 650
      real(dp) :: y(n, 3), x(n, 1), alpha(3), one(3), v(3)
      real(dp), allocatable :: info(:, :)
      type(dirichlet_regression_result_t) :: fit
      integer :: i, stat, nseed
      integer, allocatable :: seed(:)
      alpha = [2.0_dp, 3.0_dp, 5.0_dp]; one = 1.0_dp; x = 1.0_dp
      call assert_true(abs(dirichlet_pdf([0.2_dp, 0.3_dp, 0.5_dp], one) - 2.0_dp) < 1.0e-14_dp, &
         'Dirichlet uniform density')
      call dirichlet_eim_logshape(alpha, info)
      call assert_true(size(info, 1) == 3 .and. maxval(abs(info - transpose(info))) < 1.0e-14_dp, &
         'Dirichlet EIM symmetry')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 909090; call random_seed(put=seed)
      do i = 1, n
         call random_dirichlet(alpha, v, stat)
         call assert_true(stat == 0, 'Dirichlet RNG status')
         y(i, :) = v
      end do
      call fit_dirichlet_regression(y, x, fit, max_iter=300, tol=2.0e-6_dp)
      call assert_true(fit%converged, 'Dirichlet regression convergence')
      call assert_true(maxval(abs(fit%fitted_shape(1, :) - alpha)) < 0.55_dp, &
         'Dirichlet shape recovery')
      call assert_true(abs(sum(fit%fitted_mean(1, :)) - 1.0_dp) < 1.0e-14_dp, &
         'Dirichlet fitted mean normalization')
   end subroutine test_dirichlet

   subroutine test_tobit_and_folded()
      integer, parameter :: n = 750
      real(dp) :: y(n), x(n, 2), z, mu, lower, upper, pleft
      type(tobit_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      lower = 0.0_dp; upper = 1.5_dp
      pleft = pnorm_v(lower, 0.4_dp, 0.9_dp)
      call assert_true(abs(dtobit_v(lower, 0.4_dp, 0.9_dp, lower, upper) - pleft) < 1.0e-14_dp, &
         'Tobit lower point mass')
      call assert_true(qtobit_v(0.5_dp*pleft, 0.4_dp, 0.9_dp, lower, upper) == lower, &
         'Tobit lower quantile mass')
      call assert_true(abs(dfoldnorm_v(0.6_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp) - &
         2.0_dp*dnorm_v(0.6_dp, 0.0_dp, 1.0_dp)) < 1.0e-14_dp, 'folded-normal half-normal density')
      call assert_true(abs(pfoldnorm_v(qfoldnorm_v(0.37_dp, 0.3_dp, 1.1_dp, 1.3_dp, 0.8_dp), &
         0.3_dp, 1.1_dp, 1.3_dp, 0.8_dp) - 0.37_dp) < 2.0e-12_dp, 'folded-normal quantile inversion')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 909091; call random_seed(put=seed)
      do i = 1, n
         call random_number(z); z = 2.0_dp*z - 1.0_dp
         x(i, :) = [1.0_dp, z]
         mu = 0.35_dp + 0.75_dp*z
         y(i) = rtobit_v(mu, 0.85_dp, lower, upper)
      end do
      call fit_tobit(y, x, fit, lower, upper, max_iter=350, tol=2.0e-6_dp)
      call assert_true(fit%converged, 'Tobit regression convergence')
      call assert_true(maxval(abs(fit%mean_coefficients - [0.35_dp, 0.75_dp])) < 0.13_dp, &
         'Tobit mean coefficient recovery')
      call assert_true(abs(exp(fit%scale_coefficients(1)) - 0.85_dp) < 0.10_dp, &
         'Tobit scale recovery')
      call test_folded_fit()
   end subroutine test_tobit_and_folded

   subroutine test_folded_fit()
      integer, parameter :: n = 900
      real(dp) :: y(n), x(n, 1)
      type(folded_normal_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      x = 1.0_dp
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 909094; call random_seed(put=seed)
      do i = 1, n
         y(i) = rfoldnorm_v(0.6_dp, 0.8_dp, 1.4_dp, 0.7_dp)
      end do
      call fit_folded_normal(y, x, fit, 1.4_dp, 0.7_dp, max_iter=350, tol=2.0e-6_dp)
      call assert_true(fit%converged, 'folded-normal regression convergence')
      call assert_true(abs(fit%mean_coefficients(1) - 0.6_dp) < 0.14_dp, 'folded-normal mean recovery')
      call assert_true(abs(exp(fit%scale_coefficients(1)) - 0.8_dp) < 0.14_dp, 'folded-normal sd recovery')
   end subroutine test_folded_fit

   subroutine test_positive_nb()
      integer, parameter :: n = 800
      integer :: y(n), i, k, nseed
      integer, allocatable :: seed(:)
      real(dp) :: x(n, 1), s, mu0, size0
      type(positive_nb_result_t) :: fit
      mu0 = 2.5_dp; size0 = 1.8_dp; x = 1.0_dp
      s = 0.0_dp
      do k = 1, 120
         s = s + dposnbinom_v(k, mu0, size0)
      end do
      call assert_true(abs(s - 1.0_dp) < 1.0e-11_dp, 'positive NB normalization')
      k = qposnbinom_v(0.61_dp, mu0, size0)
      call assert_true(pposnbinom_v(k, mu0, size0) >= 0.61_dp, 'positive NB quantile inversion')
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 909092; call random_seed(put=seed)
      do i = 1, n
         y(i) = rposnbinom_v(mu0, size0)
      end do
      call fit_positive_negative_binomial(y, x, fit, max_iter=300, tol=2.0e-6_dp)
      call assert_true(fit%converged, 'positive NB regression convergence')
      call assert_true(abs(fit%fitted_parent_mean(1) - mu0) < 0.35_dp, 'positive NB parent mean recovery')
      call assert_true(abs(fit%size - size0) < 0.65_dp, 'positive NB size recovery')
   end subroutine test_positive_nb

   subroutine test_zoa_beta()
      integer, parameter :: n = 850
      real(dp) :: y(n), x(n, 1), mu0, phi0, p0, p1
      type(zoa_beta_result_t) :: fit
      integer :: i, nseed
      integer, allocatable :: seed(:)
      mu0 = 0.40_dp; phi0 = 7.0_dp; p0 = 0.11_dp; p1 = 0.14_dp; x = 1.0_dp
      call random_seed(size=nseed); allocate(seed(nseed)); seed = 909093; call random_seed(put=seed)
      do i = 1, n
         y(i) = rzoabeta(mu0*phi0, (1.0_dp - mu0)*phi0, p0, p1)
      end do
      call fit_zoa_beta_regression(y, x, x, fit, max_iter=350, tol=2.0e-6_dp)
      call assert_true(fit%converged, 'zero/one-altered beta convergence')
      call assert_true(abs(fit%fitted_beta_mean(1) - mu0) < 0.05_dp, 'zero/one-altered beta mean recovery')
      call assert_true(abs(fit%fitted_precision(1) - phi0) < 1.2_dp, 'zero/one-altered beta precision recovery')
      call assert_true(abs(fit%fitted_pzero(1) - p0) < 0.04_dp, 'zero/one-altered beta p0 recovery')
      call assert_true(abs(fit%fitted_pone(1) - p1) < 0.04_dp, 'zero/one-altered beta p1 recovery')
   end subroutine test_zoa_beta
end program test_v09
