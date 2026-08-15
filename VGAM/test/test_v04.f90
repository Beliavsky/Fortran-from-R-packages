program test_v04
   use vgam
   implicit none
   integer, parameter :: n = 50
   call test_drr_constraints()
   call test_rrar_nested()
   call test_gaitd_regression()
   call test_cqo_calibration()
   call test_cao_rank1()
   print '(a)', 'test_v04: PASS'
contains

   subroutine assert_true(ok, msg)
      logical, intent(in) :: ok
      character(*), intent(in) :: msg
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_drr_constraints()
      real(dp) :: x(n, 3), y(n, 3), h_a(3, 1, 2), h_c(2, 1, 2)
      real(dp) :: t, z1, z2
      integer :: i
      integer :: fam(3)
      logical :: nr(3)
      type(drrvglm_result_t) :: fit
      fam = family_gaussian
      nr = [.true., .false., .false.]
      h_a = 0.0_dp
      h_a(:, 1, 1) = [1.0_dp, 2.0_dp, 0.0_dp]
      h_a(:, 1, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
      h_c = 0.0_dp
      h_c(:, 1, 1) = [1.0_dp, 0.0_dp]
      h_c(:, 1, 2) = [0.0_dp, 1.0_dp]
      do i = 1, n
         t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
         x(i, :) = [1.0_dp, t, t*t - 0.3_dp]
         z1 = 0.7_dp*x(i, 2)
         z2 = -0.4_dp*x(i, 3)
         y(i, :) = [0.2_dp + z1, -0.3_dp + 2.0_dp*z1, 0.5_dp + z2]
      end do
      call fit_drrvglm(y, x, 2, fam, h_a, h_c, fit, no_rrr=nr, max_iter=300, tol=1.0e-9_dp)
      call assert_true(allocated(fit%loadings), 'DRR fit allocated')
      call assert_true(fit%deviance < 1.0e-12_dp, 'DRR exact constrained fit')
      call assert_true(abs(fit%loadings(2, 1) - 2.0_dp*fit%loadings(1, 1)) < 1.0e-10_dp, &
                       'DRR H.A loading constraint')
      call assert_true(maxval(abs(fit%latent_coefficients(1, :) - [1.0_dp, 0.0_dp])) < 1.0e-10_dp, &
                       'DRR first H.C constraint')
      call assert_true(maxval(abs(fit%latent_coefficients(2, :) - [0.0_dp, 1.0_dp])) < 1.0e-10_dp, &
                       'DRR second H.C constraint')
   end subroutine test_drr_constraints

   subroutine test_rrar_nested()
      integer, parameter :: nt = 100
      real(dp) :: y(nt, 2), phi0(2, 2), eps1, eps2
      real(dp), allocatable :: fc(:, :)
      integer :: i
      type(rrar_result_t) :: fit
      phi0 = reshape([0.50_dp, 0.25_dp, 0.10_dp, 0.05_dp], [2, 2])
      y = 0.0_dp
      y(1, :) = [0.6_dp, -0.2_dp]
      do i = 2, nt
         eps1 = 0.025_dp*sin(0.71_dp*real(i, dp))
         eps2 = 0.020_dp*cos(0.53_dp*real(i, dp))
         y(i, :) = matmul(phi0, y(i - 1, :)) + [eps1, eps2]
      end do
      call fit_rrar(y, [1], fit, max_iter=250, tol=1.0e-8_dp, compute_covariance=.true.)
      call assert_true(allocated(fit%phi), 'RRAR phi allocated')
      call assert_true(matrix_rank(fit%phi(:, :, 1)) == 1, 'RRAR lag rank')
      call assert_true(maxval(abs(fit%phi(:, :, 1) - phi0)) < 0.12_dp, 'RRAR coefficient recovery')
      call fit%forecast(y, 4, fc)
      call assert_true(all(abs(fc) < huge(1.0_dp)), 'RRAR forecast finite')
   end subroutine test_rrar_nested

   subroutine test_gaitd_regression()
      integer, parameter :: ng = 60
      integer :: y(ng), i
      real(dp) :: xm(ng, 2), xz(ng, 1), t
      type(gaitd_regression_result_t) :: fit
      do i = 1, ng
         t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(ng - 1, dp)
         xm(i, :) = [1.0_dp, t]
         xz(i, 1) = 1.0_dp
         if (mod(i, 4) == 0 .or. mod(i, 9) == 0) then
            y(i) = 0
         else
            y(i) = max(0, nint(exp(0.35_dp + 0.45_dp*t) + 0.35_dp*sin(real(i, dp))))
         end if
      end do
      call fit_gaitd_poisson_regression(y, xm, xz, [0], [gaitd_inflated], fit, &
                                         max_iter=250, tol=1.0e-7_dp)
      call assert_true(allocated(fit%fitted_mean), 'GAITD regression allocated')
      call assert_true(all(fit%fitted_mean >= 0.0_dp), 'GAITD fitted mean nonnegative')
      call assert_true(all(abs(fit%baseline_probability + fit%special_probabilities(:, 1) &
                           - 1.0_dp) < 1.0e-10_dp), 'GAITD mixture probabilities normalize')
      call assert_true(fit%loglik > -huge(1.0_dp)/10.0_dp, 'GAITD finite likelihood')
   end subroutine test_gaitd_regression

   subroutine test_cqo_calibration()
      integer, parameter :: nc = 45
      real(dp) :: x(nc, 2), y(nc, 2), z, t
      real(dp), allocatable :: scores(:, :), ycal(:, :), x1(:, :), dev(:)
      integer :: i
      integer :: fam(2)
      logical :: nr(2)
      type(qrrvglm_result_t) :: fit
      fam = family_gaussian
      nr = [.true., .false.]
      do i = 1, nc
         t = -1.2_dp + 2.4_dp*real(i - 1, dp)/real(nc - 1, dp)
         x(i, :) = [1.0_dp, t]
         z = 0.8_dp*t
         y(i, 1) = 0.3_dp + 1.1_dp*z - 0.7_dp*z*z
         y(i, 2) = -0.2_dp - 0.5_dp*z - 0.4_dp*z*z
      end do
      call fit_cqo(y, x, 1, fam, fit, no_rrr=nr, max_iter=35, tol=1.0e-8_dp)
      call assert_true(allocated(fit%quadratic), 'CQO fit allocated')
      allocate(x1(nc, 1)); x1 = 1.0_dp
      call cqo_calibrate(fit, y, x1, scores, max_iter=100, tol=1.0e-9_dp, deviance=dev)
      call assert_true(size(scores, 2) == 1, 'CQO calibrated scores')
      call cqo_response_surface(fit, scores, x1, ycal)
      call assert_true(maxval(abs(ycal - y)) < 5.0e-3_dp, 'CQO calibration reproduces responses')
   end subroutine test_cqo_calibration

   subroutine test_cao_rank1()
      integer, parameter :: na = 42
      real(dp) :: x(na, 2), y(na, 2), z, t, u, baseline, dev0
      real(dp), allocatable :: pred(:, :)
      integer :: i, j
      integer :: fam(2)
      type(cao_result_t) :: fit
      fam = family_gaussian
      do i = 1, na
         t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(na - 1, dp)
         u = sin(0.6_dp*real(i, dp))
         x(i, :) = [t, u]
         z = t + 0.55_dp*u
         y(i, 1) = sin(1.7_dp*z) + 0.15_dp*z
         y(i, 2) = cos(1.2_dp*z) - 0.10_dp*z*z
      end do
      dev0 = 0.0_dp
      do j = 1, 2
         baseline = sum(y(:, j))/real(na, dp)
         dev0 = dev0 + sum((y(:, j) - baseline)**2)
      end do
      call fit_cao_rank1(y, x, fam, fit, df=6, lambda=0.05_dp, max_iter=8, tol=2.0e-5_dp)
      call assert_true(allocated(fit%canonical_coefficients), 'CAO coefficients allocated')
      call assert_true(fit%deviance < 0.35_dp*dev0, 'CAO improves strongly over intercept-only')
      call fit%predict(x, pred)
      call assert_true(size(pred, 1) == na .and. size(pred, 2) == 2, 'CAO prediction shape')
   end subroutine test_cao_rank1

end program test_v04
