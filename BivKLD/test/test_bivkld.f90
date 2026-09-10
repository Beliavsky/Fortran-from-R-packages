program test_bivkld
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use bivkld, only : BIVKLD_INVALID_INPUT, BIVKLD_SUCCESS, biv_kld, biv_kld_discrete, &
                      biv_kld_independent_weibull, biv_kld_matrix, biv_kld_normal, biv_kld_pareto2, &
                      biv_sample, bivkld_estimate, dp, hscv_bandwidth
   implicit none

   call test_exact()
   call test_kernel()
   call test_pairwise()
   print '(a)', 'BivKLD tests passed.'

contains

   subroutine test_exact()
      real(dp) :: p(2, 2), q(2, 2), s1(2, 2), s2(2, 2), mean1(2), mean2(2)
      real(dp) :: alpha(2), beta(2), value

      p = reshape([0.2_dp, 0.3_dp, 0.1_dp, 0.4_dp], [2, 2])
      q = reshape([0.1_dp, 0.4_dp, 0.2_dp, 0.3_dp], [2, 2])
      call assert_close(biv_kld_discrete(p, p), 0.0_dp, 1.0e-14_dp, 'discrete identity')
      call assert_true(biv_kld_discrete(p, q) > 0.0_dp, 'discrete direction is positive')
      call assert_close(biv_kld_discrete([2.0_dp, 2.0_dp], [1.0_dp, 3.0_dp], normalize=.true.), &
                        0.5_dp*log(4.0_dp/3.0_dp), 1.0e-14_dp, 'discrete normalization')
      call assert_true(.not. ieee_is_finite(biv_kld_discrete([1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp])), &
                       'discrete zero support gives infinity')
      call assert_true(ieee_is_nan(biv_kld_discrete([0.7_dp, 0.2_dp], [0.5_dp, 0.5_dp])), &
                       'invalid probability sum gives NaN')

      mean1 = [-2.0_dp, 2.0_dp]
      mean2 = mean1
      s1 = reshape([1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
      s2 = reshape([3.0_dp, 3.0_dp, 3.0_dp, 6.0_dp], [2, 2])
      call assert_close(biv_kld_normal(mean1, s1, mean2, s2), 0.431945622001443_dp, 1.0e-12_dp, &
                        'normal reference value')
      call assert_close(biv_kld_normal([0.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
                                       [0.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])), &
                        0.0_dp, 1.0e-14_dp, 'normal identity')
      call assert_close(biv_kld_pareto2(1.0_dp, 2.0_dp), 0.401387711331890_dp, 1.0e-12_dp, 'Pareto reference')
      call assert_close(biv_kld_pareto2(2.0_dp, 2.0_dp), 0.0_dp, 1.0e-14_dp, 'Pareto identity')
      alpha = [1.0_dp, 2.0_dp]
      beta = alpha
      value = biv_kld_independent_weibull(alpha, beta)
      call assert_close(value, 0.0_dp, 1.0e-14_dp, 'Weibull identity')
   end subroutine test_exact

   subroutine test_kernel()
      type(bivkld_estimate) :: details
      real(dp) :: x(8, 2), y(8, 2), H(2, 2), Hscv(2, 2), estimate
      integer :: info

      x(:, 1) = [-1.4_dp, -1.0_dp, -0.5_dp, -0.1_dp, 0.3_dp, 0.7_dp, 1.1_dp, 1.6_dp]
      x(:, 2) = [-0.8_dp, -0.2_dp, 0.4_dp, 0.9_dp, 0.2_dp, -0.5_dp, 0.7_dp, 1.2_dp]
      y(:, 1) = x(:, 1) + 0.35_dp
      y(:, 2) = x(:, 2) - 0.20_dp
      H = reshape([0.4_dp, 0.0_dp, 0.0_dp, 0.5_dp], [2, 2])

      call biv_kld(x, x, estimate, Hx=H, Hy=H, grid_size=[30, 30], details=details, info=info)
      call assert_true(info == BIVKLD_SUCCESS, 'fixed-bandwidth identical call succeeds')
      call assert_close(estimate, 0.0_dp, 1.0e-13_dp, 'fixed-bandwidth identical KDE divergence')
      call assert_true(size(details%axis1) == 30 .and. size(details%axis2) == 30, 'details grid sizes')
      call assert_true(details%range(1, 1) < details%range(1, 2), 'details range is increasing')

      call biv_kld(x, x, estimate, Hx=H, Hy=H, grid_size=[20, 20], standardize='pooled', info=info)
      call assert_true(info == BIVKLD_SUCCESS, 'pooled-standardization call succeeds')
      call assert_close(estimate, 0.0_dp, 1.0e-13_dp, 'pooled-standardization identity')

      call biv_kld(x, y, estimate, Hx=H, Hy=H, grid_size=[9, 20], info=info)
      call assert_true(info == BIVKLD_INVALID_INPUT, 'invalid grid size is reported')
      call assert_true(ieee_is_nan(estimate), 'invalid kernel input leaves a NaN estimate')

      call biv_kld(x, y, estimate, bandwidth='normal', grid_size=[30, 30], info=info)
      call assert_true(info == BIVKLD_SUCCESS, 'normal bandwidth call succeeds')
      call assert_true(ieee_is_finite(estimate) .and. estimate >= 0.0_dp, 'normal bandwidth result is finite')

      call hscv_bandwidth(x, Hscv, info)
      call assert_true(info == BIVKLD_SUCCESS, 'SCV bandwidth selector succeeds')
      call assert_true(Hscv(1, 1) > 0.0_dp .and. Hscv(2, 2) > 0.0_dp, 'SCV diagonal is positive')
      call assert_true(Hscv(1, 1)*Hscv(2, 2) - Hscv(1, 2)*Hscv(2, 1) > 0.0_dp, 'SCV bandwidth is SPD')
      call biv_kld(x, y, estimate, bandwidth='scv', grid_size=[20, 20], info=info)
      call assert_true(info == BIVKLD_SUCCESS, 'SCV kernel divergence call succeeds')
      call assert_true(ieee_is_finite(estimate), 'SCV kernel divergence is finite')
   end subroutine test_kernel

   subroutine test_pairwise()
      type(biv_sample) :: samples(3)
      real(dp), allocatable :: divergence(:, :)
      real(dp) :: H(2, 2, 3), base(6, 2)
      integer :: info, i

      base(:, 1) = [-1.2_dp, -0.7_dp, -0.1_dp, 0.4_dp, 0.9_dp, 1.3_dp]
      base(:, 2) = [-0.4_dp, 0.5_dp, -0.8_dp, 0.9_dp, 0.1_dp, 1.1_dp]
      do i = 1, 3
         allocate(samples(i)%values(6, 2))
      end do
      samples(1)%values = base
      samples(2)%values = base
      samples(2)%values(:, 1) = samples(2)%values(:, 1) + 0.3_dp
      samples(3)%values = base
      samples(3)%values(:, 2) = samples(3)%values(:, 2) - 0.4_dp
      do i = 1, 3
         H(:, :, i) = reshape([0.35_dp, 0.0_dp, 0.0_dp, 0.45_dp], [2, 2])
      end do

      call biv_kld_matrix(samples, divergence, H=H, grid_size=[25, 25], info=info)
      call assert_true(info == BIVKLD_SUCCESS, 'pairwise call succeeds')
      call assert_true(all(shape(divergence) == [3, 3]), 'pairwise result shape')
      do i = 1, 3
         call assert_close(divergence(i, i), 0.0_dp, 1.0e-14_dp, 'pairwise diagonal')
      end do
      call assert_true(all(ieee_is_finite(divergence)), 'pairwise entries are finite')
   end subroutine test_pairwise

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed value being checked.
      real(dp), intent(in) :: expected !! Reference value expected by the test.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Short test description printed on assertion failure.

      if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', 'FAIL '//trim(label)//':', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must hold for the test to pass.
      character(len=*), intent(in) :: label !! Short test description printed on assertion failure.

      if (.not. condition) then
         print '(a)', 'FAIL '//trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_bivkld
