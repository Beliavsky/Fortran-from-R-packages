program test_gamm4
   use gamm4
   implicit none
   integer :: failures

   failures = 0
   call test_get_vb(failures)
   call test_gaussian_fit(failures)
   call test_poisson_fit(failures)
   call test_multi_penalty_rejection(failures)
   if (failures /= 0) then
      write(*, '(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*, '(a)') 'All gamm4 tests passed.'

contains

   subroutine test_get_vb(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions.
      type(gamm4_vb_result_t) :: ans
      real(dp) :: v(4), z(4, 1), phi(1, 1), xf(4, 2), xfp(4, 2), sp(2), b(2, 2)

      v = [1.0_dp, 1.2_dp, 0.8_dp, 1.1_dp]
      z(:, 1) = [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
      phi(1, 1) = 0.25_dp
      xf(:, 1) = 1.0_dp
      xf(:, 2) = [-1.0_dp, -0.2_dp, 0.4_dp, 1.0_dp]
      xfp = xf
      sp = [0.0_dp, 0.7_dp]
      b = 0.0_dp
      b(1, 1) = 1.0_dp
      b(2, 2) = 1.0_dp
      call gamm4_get_vb(v, z, phi, 1.5_dp, xf, xfp, sp, b, ans)
      call check(ans%status == 0, 'getVb status', failures)
      call check(maxval(abs(ans%vb - transpose(ans%vb))) < 1.0e-11_dp, 'getVb symmetry', failures)
      call check(all([(ans%vb(1, 1) > 0.0_dp), (ans%vb(2, 2) > 0.0_dp)]), 'getVb positive diagonal', failures)
   end subroutine test_get_vb

   subroutine test_gaussian_fit(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions.
      integer, parameter :: n = 16
      type(gamm4_smooth_t) :: smooths(1)
      type(random_term_t) :: terms(1)
      type(gamm4_result_t) :: fit
      type(gamm4_control_t) :: control
      real(dp) :: y(n), x(n), xf(n, 1), group_effect, err
      integer :: i, status

      do i = 1, n
         x(i) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
         group_effect = merge(-0.35_dp, 0.35_dp, i <= n / 2)
         err = 0.03_dp * sin(1.7_dp * real(i, dp))
         y(i) = 1.2_dp + 1.4_dp * x(i) + 0.65_dp * x(i) * x(i) + group_effect + err
      end do
      xf(:, 1) = 1.0_dp
      allocate(smooths(1)%basis(n, 2), smooths(1)%spec%penalties(2, 2, 1))
      smooths(1)%basis(:, 1) = x
      smooths(1)%basis(:, 2) = x * x
      smooths(1)%spec%penalties = 0.0_dp
      smooths(1)%spec%penalties(2, 2, 1) = 1.0_dp
      smooths(1)%label = 's(x)'
      allocate(terms(1)%z(n, 1), terms(1)%group(n))
      terms(1)%z = 1.0_dp
      terms(1)%group(1:n / 2) = 1
      terms(1)%group(n / 2 + 1:n) = 2
      terms(1)%n_levels = 2
      terms(1)%covariance_structure = covariance_diagonal
      terms(1)%name = 'group'
      control%max_outer = 6
      control%tolerance = 2.0e-4_dp
      call gamm4_fit(y, xf, smooths, fit, status, random_terms=terms, reml=.false., control=control)
      call check(status == 0 .and. fit%converged, 'Gaussian GAMM convergence', failures)
      if (status == 0) then
         call check(size(fit%coefficients) == 3, 'Gaussian coefficient dimension', failures)
         call check(fit%sp(1) > 0.0_dp, 'Gaussian smoothing parameter', failures)
         call check(fit%scale > 0.0_dp, 'Gaussian scale', failures)
         call check(sum(fit%residuals * fit%residuals) < 0.20_dp, 'Gaussian conditional residuals', failures)
         call check(maxval(abs(fit%conditional_fitted - fit%fitted)) > 0.05_dp, 'ordinary random effects retained', failures)
         call check(all(fit%edf >= -1.0e-8_dp), 'Gaussian nonnegative EDF contributions', failures)
      end if
   end subroutine test_gaussian_fit

   subroutine test_poisson_fit(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions.
      integer, parameter :: n = 12
      type(gamm4_smooth_t) :: smooths(1)
      type(gamm4_result_t) :: fit
      type(gamm4_control_t) :: control
      real(dp) :: y(n), x(n), xf(n, 1), lambda
      integer :: i, status

      do i = 1, n
         x(i) = -0.9_dp + 1.8_dp * real(i - 1, dp) / real(n - 1, dp)
         lambda = exp(0.7_dp + 0.4_dp * x(i) + 0.25_dp * x(i) * x(i))
         y(i) = real(max(0, nint(lambda + 0.35_dp * sin(real(i, dp)))), dp)
      end do
      xf(:, 1) = 1.0_dp
      allocate(smooths(1)%basis(n, 2), smooths(1)%spec%penalties(2, 2, 1))
      smooths(1)%basis(:, 1) = x
      smooths(1)%basis(:, 2) = x * x
      smooths(1)%spec%penalties = 0.0_dp
      smooths(1)%spec%penalties(2, 2, 1) = 1.0_dp
      smooths(1)%label = 's(x)'
      control%max_outer = 5
      control%max_pirls = 80
      control%tolerance = 5.0e-4_dp
      call gamm4_fit(y, xf, smooths, fit, status, family=family_poisson, control=control)
      call check(status == 0 .and. fit%converged, 'Poisson GAMM convergence', failures)
      if (status == 0) then
         call check(all(fit%fitted > 0.0_dp), 'Poisson positive fitted means', failures)
         call check(fit%sp(1) > 0.0_dp, 'Poisson smoothing parameter', failures)
         call check(fit%method == 'glmer.Laplace', 'Poisson method label', failures)
      end if
   end subroutine test_poisson_fit

   subroutine test_multi_penalty_rejection(failures)
      integer, intent(inout) :: failures !! Running count of failed assertions.
      type(gamm4_smooth_t) :: smooths(1)
      type(gamm4_result_t) :: fit
      real(dp) :: y(4), xf(4, 1)
      integer :: status

      y = [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp]
      xf = 1.0_dp
      allocate(smooths(1)%basis(4, 2), smooths(1)%spec%penalties(2, 2, 2))
      smooths(1)%basis(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
      smooths(1)%basis(:, 2) = smooths(1)%basis(:, 1) ** 2
      smooths(1)%spec%penalties = 0.0_dp
      smooths(1)%spec%penalties(1, 1, 1) = 1.0_dp
      smooths(1)%spec%penalties(2, 2, 2) = 1.0_dp
      call gamm4_fit(y, xf, smooths, fit, status)
      call check(status /= 0, 'multi-penalty smooth rejection', failures)
   end subroutine test_multi_penalty_rejection

   subroutine check(condition, label, failures)
      logical, intent(in) :: condition !! Assertion condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed on failure.
      integer, intent(inout) :: failures !! Running number of failed assertions.

      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // trim(label)
      end if
   end subroutine check

end program test_gamm4
