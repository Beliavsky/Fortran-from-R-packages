program test_kde1d
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use kde1d_api, only : dkde1d, dp, equi_jitter, kde1d_fit, kde1d_loglik
   use kde1d_api, only : kde1d_model, kde1d_summary, pkde1d, qkde1d, rkde1d
   implicit none

   call test_jitter()
   call test_continuous_fixed()
   call test_continuous_automatic()
   call test_automatic_scenarios()
   call test_bounded()
   call test_boundary_repair_effect()
   call test_scale_equivariance()
   call test_discrete()
   call test_zero_inflated()
   call test_zero_inflated_bounded_component()
   call test_zero_inflated_all_zero()
   call test_weighted()
   print '(a)', 'All kde1d tests passed.'

contains

   pure real(dp) function d1(x, model) result(value)
      real(dp), intent(in) :: x !! Single point at which the fitted density or mass is evaluated.
      type(kde1d_model), intent(in) :: model !! Fitted model passed to dkde1d.

      real(dp) :: work(1)

      work = dkde1d([x], model)
      value = work(1)
   end function d1

   pure real(dp) function p1(x, model) result(value)
      real(dp), intent(in) :: x !! Single point at which the fitted CDF is evaluated.
      type(kde1d_model), intent(in) :: model !! Fitted model passed to pkde1d.

      real(dp) :: work(1)

      work = pkde1d([x], model)
      value = work(1)
   end function p1

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must hold for the test to continue.
      character(len=*), intent(in) :: message !! Diagnostic printed if the assertion fails.

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Computed scalar to compare with the expected value.
      real(dp), intent(in) :: expected !! Reference scalar used by the assertion.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference.
      character(len=*), intent(in) :: message !! Diagnostic printed if the assertion fails.

      call assert_true(abs(actual - expected) <= tolerance, message)
   end subroutine assert_close

   subroutine test_jitter()
      real(dp) :: y(5)

      y = equi_jitter([1, 1, 2, 2, 2])
      call assert_close(y(1), 5.0_dp / 6.0_dp, 1.0e-12_dp, 'first factor jitter')
      call assert_close(y(2), 7.0_dp / 6.0_dp, 1.0e-12_dp, 'second factor jitter')
      call assert_close(y(3), 1.75_dp, 1.0e-12_dp, 'third factor jitter')
      call assert_close(y(4), 2.0_dp, 1.0e-12_dp, 'fourth factor jitter')
      call assert_close(y(5), 2.25_dp, 1.0e-12_dp, 'fifth factor jitter')
      call assert_true(all(equi_jitter([1.0_dp, 2.0_dp]) == [1.0_dp, 2.0_dp]), &
         'numeric equi_jitter must be unchanged')
   end subroutine test_jitter

   subroutine test_continuous_fixed()
      type(kde1d_model) :: model
      real(dp) :: cdf(5)
      real(dp) :: dens(5)
      real(dp) :: draws(16)
      real(dp) :: q(5)
      real(dp) :: s(5)
      real(dp) :: x(81)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         x(i) = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(size(x) - 1, dp) &
            + 0.08_dp * sin(real(i, dp))
      end do
      call kde1d_fit(x, model, ierr, bw=0.35_dp, deg=2, boundary_repair=.false.)
      call assert_true(ierr == 0, 'fixed-bandwidth continuous fit')
      dens = dkde1d([-1.5_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.5_dp], model)
      call assert_true(all(dens >= 0.0_dp), 'continuous density is nonnegative')
      call assert_true(all(ieee_is_finite(dens)), 'continuous density is finite')
      cdf = pkde1d([-10.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 10.0_dp], model)
      call assert_true(all(cdf(2:) >= cdf(:4)), 'continuous CDF is nondecreasing')
      call assert_close(cdf(1), 0.0_dp, 1.0e-10_dp, 'continuous left CDF tail')
      call assert_close(cdf(5), 1.0_dp, 1.0e-10_dp, 'continuous right CDF tail')
      q = qkde1d([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], model)
      call assert_true(all(q(2:) >= q(:4)), 'continuous quantiles are nondecreasing')
      call assert_close(p1(q(3), model), 0.5_dp, 3.0e-3_dp, 'continuous median inversion')
      call rkde1d(size(draws), model, draws, seed=17)
      call assert_true(all(ieee_is_finite(draws)), 'pseudorandom draws are finite')
      call rkde1d(size(draws), model, draws, quasi=.true., seed=17)
      call assert_true(all(ieee_is_finite(draws)), 'quasi-random draws are finite')
      s = kde1d_summary(model)
      call assert_close(s(1), real(size(x), dp), 0.0_dp, 'summary nobs')
      call assert_close(kde1d_loglik(model), s(4), 0.0_dp, 'summary loglik')
   end subroutine test_continuous_fixed

   subroutine test_continuous_automatic()
      type(kde1d_model) :: model
      real(dp) :: x(61)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         x(i) = cos(0.17_dp * real(i, dp)) + 0.35_dp * sin(0.63_dp * real(i, dp))
      end do
      call kde1d_fit(x, model, ierr, deg=2, boundary_repair=.false.)
      call assert_true(ierr == 0, 'automatic-bandwidth continuous fit')
      call assert_true(model%bw > 0.0_dp .and. ieee_is_finite(model%bw), 'automatic bandwidth is positive')
      call assert_true(ieee_is_finite(model%loglik), 'automatic fit loglik is finite')
   end subroutine test_continuous_automatic

   subroutine test_automatic_scenarios()
      type(kde1d_model) :: model
      real(dp) :: x(48)
      integer :: degree
      integer :: i
      integer :: ierr

      do degree = 0, 2
         do i = 1, size(x)
            x(i) = sin(0.19_dp * real(i, dp)) + 0.02_dp * real(i, dp)
         end do
         call kde1d_fit(x, model, ierr, deg=degree, boundary_repair=.false.)
         call assert_true(ierr == 0, 'automatic unbounded fit for degrees 0 through 2')

         do i = 1, size(x)
            x(i) = 0.01_dp + 2.0_dp * (real(i, dp) / real(size(x), dp))**2
         end do
         call kde1d_fit(x, model, ierr, xmin=0.0_dp, deg=degree, boundary_repair=.false.)
         call assert_true(ierr == 0, 'automatic left-bounded fit for degrees 0 through 2')

         x = -x
         call kde1d_fit(x, model, ierr, xmax=0.0_dp, deg=degree, boundary_repair=.false.)
         call assert_true(ierr == 0, 'automatic right-bounded fit for degrees 0 through 2')

         do i = 1, size(x)
            x(i) = 0.01_dp + 0.98_dp * real(i - 1, dp) / real(size(x) - 1, dp)
         end do
         call kde1d_fit(x, model, ierr, xmin=0.0_dp, xmax=1.0_dp, deg=degree, boundary_repair=.false.)
         call assert_true(ierr == 0, 'automatic two-sided fit for degrees 0 through 2')
      end do
   end subroutine test_automatic_scenarios

   subroutine test_bounded()
      type(kde1d_model) :: left_model
      type(kde1d_model) :: both_model
      real(dp) :: xleft(40)
      real(dp) :: xunit(40)
      real(dp) :: p(4)
      integer :: i
      integer :: ierr

      do i = 1, 40
         xleft(i) = 0.015_dp + 2.5_dp * (real(i, dp) / 40.0_dp)**2
         xunit(i) = 0.005_dp + 0.99_dp * real(i - 1, dp) / 39.0_dp
      end do
      call kde1d_fit(xleft, left_model, ierr, xmin=0.0_dp, bw=0.32_dp, deg=1, boundary_repair=.true.)
      call assert_true(ierr == 0, 'one-sided bounded fit with boundary repair')
      call assert_close(d1(-0.1_dp, left_model), 0.0_dp, 0.0_dp, 'left support density')
      p = pkde1d([-0.1_dp, 0.0_dp, 1.0_dp, 10.0_dp], left_model)
      call assert_close(p(1), 0.0_dp, 1.0e-12_dp, 'left support CDF')
      call assert_close(p(4), 1.0_dp, 1.0e-10_dp, 'left bounded CDF tail')

      call kde1d_fit(xunit, both_model, ierr, xmin=0.0_dp, xmax=1.0_dp, bw=0.30_dp, deg=2, &
         boundary_repair=.true.)
      call assert_true(ierr == 0, 'two-sided bounded fit with boundary repair')
      call assert_close(d1(-0.1_dp, both_model), 0.0_dp, 0.0_dp, 'two-sided lower support density')
      call assert_close(d1(1.1_dp, both_model), 0.0_dp, 0.0_dp, 'two-sided upper support density')
      p = pkde1d([-0.1_dp, 0.0_dp, 1.0_dp, 1.1_dp], both_model)
      call assert_close(p(1), 0.0_dp, 1.0e-12_dp, 'two-sided CDF below support')
      call assert_close(p(4), 1.0_dp, 1.0e-10_dp, 'two-sided CDF above support')
   end subroutine test_bounded

   subroutine test_boundary_repair_effect()
      type(kde1d_model) :: bulk
      type(kde1d_model) :: repaired
      real(dp) :: p
      real(dp) :: x(80)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         p = (real(i, dp) - 0.5_dp) / real(size(x), dp)
         x(i) = -log(1.0_dp - p)
      end do
      call kde1d_fit(x, repaired, ierr, xmin=0.0_dp, bw=0.35_dp, deg=2, boundary_repair=.true.)
      call assert_true(ierr == 0, 'repaired endpoint fit')
      call kde1d_fit(x, bulk, ierr, xmin=0.0_dp, bw=0.35_dp, deg=2, boundary_repair=.false.)
      call assert_true(ierr == 0, 'bulk endpoint fit')
      call assert_true(maxval(abs(repaired%grid%values - bulk%grid%values)) > 1.0e-8_dp, &
         'boundary repair changes eligible endpoint estimate')
   end subroutine test_boundary_repair_effect

   subroutine test_scale_equivariance()
      type(kde1d_model) :: base
      type(kde1d_model) :: scaled
      real(dp), parameter :: factor = 1.0e4_dp
      real(dp) :: d1v(7)
      real(dp) :: d2v(7)
      real(dp) :: eval(7)
      real(dp) :: p(7)
      real(dp) :: q1(7)
      real(dp) :: q2(7)
      real(dp) :: x(64)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         x(i) = -log(1.0_dp - (real(i, dp) - 0.5_dp) / real(size(x), dp))
      end do
      p = [0.08_dp, 0.20_dp, 0.35_dp, 0.50_dp, 0.65_dp, 0.80_dp, 0.92_dp]
      eval = -log(1.0_dp - p)
      call kde1d_fit(x, base, ierr, xmin=0.0_dp, deg=1, boundary_repair=.false.)
      call assert_true(ierr == 0, 'base scale-equivariance fit')
      call kde1d_fit(factor * x, scaled, ierr, xmin=0.0_dp, deg=1, boundary_repair=.false.)
      call assert_true(ierr == 0, 'scaled scale-equivariance fit')
      d1v = dkde1d(eval, base)
      d2v = factor * dkde1d(factor * eval, scaled)
      call assert_true(maxval(abs(d1v - d2v)) < 2.0e-8_dp, 'left-bound density scale equivariance')
      call assert_true(maxval(abs(pkde1d(eval, base) - pkde1d(factor * eval, scaled))) < 2.0e-8_dp, &
         'left-bound CDF scale equivariance')
      q1 = qkde1d(p, base)
      q2 = qkde1d(p, scaled) / factor
      call assert_true(maxval(abs(q1 - q2)) < 2.0e-8_dp, 'left-bound quantile scale equivariance')

      call kde1d_fit(-x, base, ierr, xmax=0.0_dp, deg=1, boundary_repair=.false.)
      call assert_true(ierr == 0, 'right-bound scale-equivariance fit')
      call kde1d_fit(-factor * x, scaled, ierr, xmax=0.0_dp, deg=1, boundary_repair=.false.)
      call assert_true(ierr == 0, 'scaled right-bound scale-equivariance fit')
      d1v = dkde1d(-eval, base)
      d2v = factor * dkde1d(-factor * eval, scaled)
      call assert_true(maxval(abs(d1v - d2v)) < 2.0e-8_dp, 'right-bound density scale equivariance')
      call assert_true(maxval(abs(pkde1d(-eval, base) - pkde1d(-factor * eval, scaled))) < 2.0e-8_dp, &
         'right-bound CDF scale equivariance')
   end subroutine test_scale_equivariance

   subroutine test_discrete()
      type(kde1d_model) :: model
      real(dp) :: mass(5)
      real(dp) :: q(5)
      real(dp) :: x(30)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         x(i) = real(mod(i - 1, 5), dp)
      end do
      call kde1d_fit(x, model, ierr, type_name='discrete', xmin=0.0_dp, xmax=4.0_dp, &
         bw=0.45_dp, deg=1, boundary_repair=.false.)
      call assert_true(ierr == 0, 'discrete fit')
      mass = dkde1d([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], model)
      call assert_close(sum(mass), 1.0_dp, 2.0e-12_dp, 'discrete masses normalize')
      call assert_close(d1(0.5_dp, model), 0.0_dp, 0.0_dp, 'noninteger discrete mass')
      call assert_close(p1(-1.0_dp, model), 0.0_dp, 0.0_dp, 'discrete CDF below support')
      call assert_close(p1(4.0_dp, model), 1.0_dp, 1.0e-12_dp, 'discrete CDF at upper support')
      q = qkde1d([0.0_dp, 0.2_dp, 0.5_dp, 0.8_dp, 1.0_dp], model)
      call assert_true(all(q == anint(q)), 'discrete quantiles are integer levels')
      call assert_true(all(q >= 0.0_dp .and. q <= 4.0_dp), 'discrete quantiles respect support')
   end subroutine test_discrete

   subroutine test_zero_inflated()
      type(kde1d_model) :: model
      real(dp) :: x(50)
      real(dp) :: p(4)
      real(dp) :: q(5)
      integer :: i
      integer :: ierr

      x(:10) = 0.0_dp
      do i = 11, size(x)
         x(i) = 0.05_dp + real(i - 10, dp) / 20.0_dp
      end do
      call kde1d_fit(x, model, ierr, type_name='zero-inflated', xmin=0.0_dp, bw=0.30_dp, deg=1, &
         boundary_repair=.false.)
      call assert_true(ierr == 0, 'zero-inflated fit')
      call assert_close(model%prob0, 0.2_dp, 1.0e-12_dp, 'zero mass estimate')
      call assert_close(d1(0.0_dp, model), 0.2_dp, 1.0e-12_dp, 'zero-inflated density at atom')
      p = pkde1d([-0.1_dp, 0.0_dp, 0.5_dp, 10.0_dp], model)
      call assert_close(p(1), 0.0_dp, 1.0e-10_dp, 'zero-inflated CDF below zero')
      call assert_true(p(2) >= 0.2_dp - 1.0e-12_dp, 'zero-inflated CDF includes atom')
      call assert_close(p(4), 1.0_dp, 1.0e-10_dp, 'zero-inflated CDF right tail')
      q = qkde1d([0.05_dp, 0.1_dp, 0.2_dp, 0.5_dp, 0.9_dp], model)
      call assert_true(all(q(:2) == 0.0_dp), 'probabilities strictly inside the hurdle mass quantile to zero')
      call assert_true(q(3) >= 0.0_dp, 'hurdle-boundary quantile is nonnegative')
      call assert_true(q(4) > 0.0_dp .and. q(5) > q(4), 'continuous hurdle quantiles are positive')
   end subroutine test_zero_inflated

   subroutine test_zero_inflated_bounded_component()
      type(kde1d_model) :: model
      real(dp) :: x(50)
      integer :: i
      integer :: ierr

      x(:10) = 0.0_dp
      do i = 11, size(x)
         x(i) = 1.01_dp + 0.98_dp * real(i - 11, dp) / real(size(x) - 11, dp)
      end do
      call kde1d_fit(x, model, ierr, type_name='zero_inflated', xmin=1.0_dp, xmax=2.0_dp, &
         bw=0.25_dp, boundary_repair=.false.)
      call assert_true(ierr == 0, 'zero-inflated fit with bounded continuous component')
      call assert_close(d1(0.0_dp, model), 0.2_dp, 1.0e-12_dp, 'bounded zero-inflated atom')
      call assert_close(d1(0.5_dp, model), 0.0_dp, 0.0_dp, 'continuous density outside bounded component')
      call assert_close(p1(0.5_dp, model), 0.2_dp, 1.0e-12_dp, 'CDF between atom and continuous support')
   end subroutine test_zero_inflated_bounded_component

   subroutine test_zero_inflated_all_zero()
      type(kde1d_model) :: model
      real(dp) :: q(3)
      integer :: ierr

      call kde1d_fit([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], model, ierr, type_name='zero-inflated')
      call assert_true(ierr == 0, 'all-zero zero-inflated fit')
      call assert_close(model%prob0, 1.0_dp, 0.0_dp, 'all-zero atom probability')
      call assert_true(ieee_is_nan(model%bw), 'all-zero bandwidth is NaN')
      call assert_close(d1(0.0_dp, model), 1.0_dp, 0.0_dp, 'all-zero atom density')
      call assert_close(p1(-1.0_dp, model), 0.0_dp, 0.0_dp, 'all-zero left CDF')
      call assert_close(p1(0.0_dp, model), 1.0_dp, 0.0_dp, 'all-zero CDF at atom')
      q = qkde1d([0.0_dp, 0.5_dp, 1.0_dp], model)
      call assert_true(all(q == 0.0_dp), 'all-zero quantiles')
   end subroutine test_zero_inflated_all_zero

   subroutine test_weighted()
      type(kde1d_model) :: model
      real(dp) :: w(24)
      real(dp) :: x(24)
      integer :: i
      integer :: ierr

      do i = 1, size(x)
         x(i) = 0.1_dp * real(i, dp)
      end do
      w = 1.0_dp
      w(1:4) = 0.0_dp
      w(21:24) = 2.0_dp
      call kde1d_fit(x, model, ierr, bw=0.25_dp, deg=0, weights=w, boundary_repair=.false.)
      call assert_true(ierr == 0, 'weighted fit')
      call assert_true(ieee_is_finite(model%loglik), 'weighted loglik is finite')
      call assert_true(model%edf >= 0.0_dp, 'weighted EDF is nonnegative')
   end subroutine test_weighted

end program test_kde1d
