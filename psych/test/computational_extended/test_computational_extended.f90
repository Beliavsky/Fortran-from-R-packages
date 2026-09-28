program test_computational_extended
   use r_kinds, only : dp
   use psych_computational_extended
   use psych_omega_sem, only : omega_sem_bifactor
   implicit none
   integer :: failures

   failures = 0
   call test_multilevel(failures)
   call test_factor_poly(failures)
   call test_factor_multi(failures)
   call test_irt_factor(failures)
   call test_omega_cfa(failures)
   call test_omega_sem(failures)
   call test_circadian_reliability(failures)
   call test_cross_validation_bootstrap(failures)
   call test_factor_solution_correlation(failures)
   call test_organize_factor_loadings(failures)
   call test_parallel_select(failures)
   call test_vss_select(failures)
   if (failures /= 0) error stop 1
   print *, 'computational extended tests passed'

contains

   subroutine check(condition, message, failures)
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message identifying the assertion.
      integer, intent(inout) :: failures !! Running assertion-failure count.
      if (.not. condition) then
         failures = failures + 1
         print *, 'FAIL: ', trim(message)
      end if
   end subroutine check

   subroutine test_multilevel(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: values(3, 3, 3)
      real(dp), allocatable :: components(:), reliabilities(:), mean_squares(:)
      integer :: i, j, k, info

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               values(i, j, k) = 10.0_dp + 2.0_dp * real(i - 2, dp) + 1.5_dp * real(j - 2, dp) + &
                  0.8_dp * real(k - 2, dp) + 0.4_dp * real((i - 2) * (j - 2), dp) + &
                  0.3_dp * real((i - 2) * (k - 2), dp) + 0.2_dp * real((j - 2) * (k - 2), dp) + &
                  0.1_dp * real((i - 2) * (j - 2) * (k - 2), dp)
            end do
         end do
      end do
      call multilevel_reliability_balanced(values, components, reliabilities, mean_squares, info)
      call check(info == 0, 'multilevel info', failures)
      call check(maxval(abs(mean_squares - [36.0_dp, 20.25_dp, 5.76_dp, 0.48_dp, 0.27_dp, 0.12_dp, &
         0.01_dp])) < 1.0e-10_dp, 'multilevel mean squares', failures)
      call check(abs(reliabilities(1) - 0.99971855_dp) < 1.0e-7_dp, 'multilevel RkF', failures)
      call check(abs(reliabilities(4) - 0.97916667_dp) < 1.0e-7_dp, 'multilevel Rc', failures)
   end subroutine test_multilevel

   subroutine test_factor_poly(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: data(40, 5), x
      real(dp), allocatable :: rho(:, :), loadings(:, :), phi(:, :), uniqueness(:)
      real(dp) :: objective
      integer :: i, info

      do i = 1, 40
         x = real(i, dp)
         data(i, 1) = sin(0.17_dp * x) + 0.03_dp * x
         data(i, 2) = 0.8_dp * data(i, 1) + 0.2_dp * cos(0.31_dp * x)
         data(i, 3) = cos(0.11_dp * x) - 0.02_dp * x
         data(i, 4) = 0.75_dp * data(i, 3) + 0.25_dp * sin(0.23_dp * x)
         data(i, 5) = 0.4_dp * data(i, 1) + 0.4_dp * data(i, 3) + 0.2_dp * sin(0.41_dp * x)
      end do
      call factor_analysis_poly_data(data, 2, rho, loadings, phi, uniqueness, objective, info)
      call check(info == 0, 'fa.poly info', failures)
      call check(size(loadings, 1) == 5 .and. size(loadings, 2) == 2, 'fa.poly loading shape', failures)
      call check(maxval(abs(rho - transpose(rho))) < 1.0e-12_dp, 'fa.poly rho symmetry', failures)
      call check(maxval(abs([(rho(i, i) - 1.0_dp, i = 1, 5)])) < 1.0e-12_dp, 'fa.poly rho diagonal', failures)
   end subroutine test_factor_poly

   subroutine test_factor_multi(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: l(6, 2), phi0(2, 2), r(6, 6)
      real(dp), allocatable :: first(:, :), phi(:, :), second(:, :)
      integer :: i, info

      l = 0.0_dp
      l(1:3, 1) = [0.80_dp, 0.75_dp, 0.70_dp]
      l(4:6, 2) = [0.78_dp, 0.73_dp, 0.68_dp]
      phi0 = reshape([1.0_dp, 0.35_dp, 0.35_dp, 1.0_dp], [2, 2])
      r = matmul(l, matmul(phi0, transpose(l)))
      do i = 1, 6
         r(i, i) = 1.0_dp
      end do
      call factor_analysis_multi(r, 2, 1, first, phi, second, info)
      call check(info == 0, 'fa.multi info', failures)
      call check(size(first, 1) == 6 .and. size(first, 2) == 2, 'fa.multi first shape', failures)
      call check(size(second, 1) == 2 .and. size(second, 2) == 1, 'fa.multi second shape', failures)
   end subroutine test_factor_multi

   subroutine test_irt_factor(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      integer :: items(120, 4)
      real(dp), allocatable :: rho(:, :), thresholds(:, :), loadings(:, :), discrimination(:, :)
      real(dp), allocatable :: difficulty(:, :, :), uniqueness(:)
      real(dp) :: z, objective
      integer :: i, j, info

      do i = 1, 120
         z = (real(i, dp) - 60.5_dp) / 25.0_dp
         do j = 1, 4
            items(i, j) = merge(1, 0, z + 0.25_dp * sin(real(i * j, dp)) > -0.35_dp + 0.22_dp * real(j - 1, dp))
         end do
      end do
      call irt_factor_analysis(items, 0, 1, 1, rho, thresholds, loadings, discrimination, difficulty, uniqueness, &
         objective, info)
      call check(info == 0, 'irt.fa info', failures)
      call check(size(thresholds, 2) == 1, 'irt.fa threshold count', failures)
      call check(size(discrimination, 1) == 4 .and. size(discrimination, 2) == 1, 'irt.fa discrimination shape', failures)
      call check(all(abs(discrimination) < huge(1.0_dp)), 'irt.fa finite discrimination', failures)
   end subroutine test_irt_factor

   subroutine test_omega_cfa(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: l(6, 3), r(6, 6)
      real(dp), allocatable :: r2(:), groups(:, :)
      real(dp) :: omega_h, omega_total
      integer :: i, info

      l = 0.0_dp
      l(:, 1) = 0.50_dp
      l(1:3, 2) = [0.45_dp, 0.40_dp, 0.35_dp]
      l(4:6, 3) = [0.42_dp, 0.38_dp, 0.34_dp]
      r = matmul(l, transpose(l))
      do i = 1, 6
         r(i, i) = 1.0_dp
      end do
      call omega_from_cfa_loadings(r, l, omega_h, omega_total, r2, groups, info)
      call check(info == 0, 'omegaFromSem info', failures)
      call check(abs(omega_h - sum(l(:, 1))**2 / sum(r)) < 1.0e-10_dp, 'omegaFromSem omega_h', failures)
      call check(omega_total >= omega_h .and. omega_total <= 1.0_dp + 1.0e-12_dp, 'omegaFromSem omega_total', failures)
      call check(all(r2 > 0.0_dp), 'omegaFromSem score R2', failures)
   end subroutine test_omega_cfa

   subroutine test_omega_sem(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: r1(6, 6), r2model(6, 6), l1(6, 1), l2(6, 3), psi(6)
      real(dp), allocatable :: fitted(:, :), uniqueness(:), score_r2(:), omega_group(:, :)
      real(dp) :: omega_h, omega_total, fit, chi_square
      integer :: i, info, iterations

      l1(:, 1) = [0.62_dp, 0.57_dp, 0.53_dp, 0.66_dp, 0.59_dp, 0.55_dp]
      r1 = matmul(l1, transpose(l1))
      do i = 1, 6
         r1(i, i) = 1.0_dp
      end do
      call omega_sem_bifactor(r1, 1, fitted, uniqueness, omega_h, omega_total, score_r2, omega_group, fit, &
         iterations, info, n_obs=400, chi_square=chi_square, max_iter=700, tolerance=1.0e-8_dp)
      call check(info == 0, 'omegaSem one-factor convergence', failures)
      call check(fit < 1.0e-7_dp, 'omegaSem one-factor exact ML fit', failures)
      call check(size(fitted, 2) == 1, 'omegaSem one-factor loading shape', failures)
      call check(all(uniqueness > 0.0_dp), 'omegaSem positive uniqueness', failures)
      call check(abs(chi_square - 399.0_dp * fit) < 1.0e-10_dp, 'omegaSem chi-square scaling', failures)

      l2 = 0.0_dp
      l2(:, 1) = [0.55_dp, 0.50_dp, 0.45_dp, 0.60_dp, 0.52_dp, 0.48_dp]
      l2(1:3, 2) = [0.45_dp, 0.50_dp, 0.55_dp]
      l2(4:6, 3) = [0.40_dp, 0.47_dp, 0.53_dp]
      do i = 1, 6
         psi(i) = 1.0_dp - sum(l2(i, :) ** 2)
      end do
      r2model = matmul(l2, transpose(l2))
      do i = 1, 6
         r2model(i, i) = r2model(i, i) + psi(i)
      end do
      call omega_sem_bifactor(r2model, 2, fitted, uniqueness, omega_h, omega_total, score_r2, omega_group, fit, &
         iterations, info, max_iter=800, tolerance=1.0e-8_dp)
      call check(info == 0, 'omegaSem bifactor convergence', failures)
      call check(fit < 1.0e-6_dp, 'omegaSem bifactor reproduction', failures)
      call check(size(fitted, 1) == 6 .and. size(fitted, 2) == 3, 'omegaSem bifactor loading shape', failures)
      call check(omega_h > 0.0_dp .and. omega_total >= omega_h, 'omegaSem omega ordering', failures)
      call check(all(score_r2 > 0.0_dp), 'omegaSem score determinacy', failures)
   end subroutine test_omega_sem

   subroutine test_circadian_reliability(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: angle(10, 48), responses(10, 48, 2), phase, noise
      real(dp), allocatable :: phase_rel(:), fit_rel(:), splits(:, :, :)
      integer :: s, i, j, info
      real(dp), parameter :: pi = acos(-1.0_dp)

      do s = 1, 10
         phase = 1.2_dp + 0.7_dp * real(s - 1, dp)
         do i = 1, 48
            angle(s, i) = 0.5_dp * real(i - 1, dp)
            do j = 1, 2
               noise = 0.05_dp * real(s, dp) * sin(0.37_dp * real(i + 2 * j, dp))
               responses(s, i, j) = cos(2.0_dp * pi * (angle(s, i) - phase - 0.6_dp * real(j - 1, dp)) / 24.0_dp) + noise
            end do
         end do
      end do
      call circadian_reliability_split(angle, responses, phase_rel, fit_rel, splits, odd_even=.true., hours=.true., &
         period=24.0_dp, minimum_observations=16, info=info)
      call check(info == 0, 'circadian reliability info', failures)
      call check(all(phase_rel > 0.8_dp), 'circadian phase reliability', failures)
      call check(size(splits, 1) == 10 .and. size(splits, 2) == 2, 'circadian split shape', failures)
   end subroutine test_circadian_reliability

   subroutine test_cross_validation_bootstrap(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: x(80, 2), y(80, 1), z
      real(dp), allocatable :: trial(:, :), means(:), coefficients(:, :)
      integer, allocatable :: seed(:)
      integer :: i, nseed, info

      do i = 1, 80
         z = real(i, dp)
         x(i, 1) = sin(0.17_dp * z) + 0.015_dp * z
         x(i, 2) = cos(0.11_dp * z) - 0.012_dp * z
         y(i, 1) = 0.75_dp * x(i, 1) - 0.45_dp * x(i, 2) + 0.08_dp * sin(0.37_dp * z)
      end do
      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed = 13579
      call random_seed(put=seed)
      call cross_validation_bootstrap(y, x, 30, trial, means, coefficients, info)
      call check(info == 0, 'crossValidationBoot info', failures)
      call check(means(1) > 0.95_dp, 'crossValidationBoot derivation fit', failures)
      call check(means(2) > 0.90_dp, 'crossValidationBoot validation fit', failures)
      call check(size(coefficients, 1) == 2, 'crossValidationBoot coefficient shape', failures)
   end subroutine test_cross_validation_bootstrap


   subroutine test_factor_solution_correlation(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: l(6, 2), phi0(2, 2), r(6, 6)
      real(dp), allocatable :: rf(:, :), congruence(:, :), first(:, :), second(:, :), v1(:), v2(:)
      integer :: i, info

      l = 0.0_dp
      l(1:3, 1) = [0.82_dp, 0.76_dp, 0.71_dp]
      l(4:6, 2) = [0.79_dp, 0.74_dp, 0.69_dp]
      phi0 = reshape([1.0_dp, 0.30_dp, 0.30_dp, 1.0_dp], [2, 2])
      r = matmul(l, matmul(phi0, transpose(l)))
      do i = 1, 6
         r(i, i) = 1.0_dp
      end do
      call factor_solution_correlation(r, 2, 2, rf, congruence, first, second, v1, v2, info)
      call check(info == 0, 'faCor info', failures)
      call check(abs(rf(1, 1) - 1.0_dp) < 1.0e-10_dp .and. abs(rf(2, 2) - 1.0_dp) < 1.0e-10_dp, &
         'faCor identical solution correlation', failures)
      call check(abs(congruence(1, 1) - 1.0_dp) < 1.0e-10_dp .and. abs(congruence(2, 2) - 1.0_dp) < 1.0e-10_dp, &
         'faCor identical congruence', failures)
      call check(all(v1 > 0.0_dp) .and. all(v2 > 0.0_dp), 'faCor positive score variances', failures)
   end subroutine test_factor_solution_correlation

   subroutine test_organize_factor_loadings(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: loadings(3, 2)
      real(dp), allocatable :: organized(:, :), signs(:)
      integer :: info

      loadings = reshape([-0.8_dp, -0.4_dp, 0.1_dp, 0.2_dp, 0.7_dp, 0.5_dp], [3, 2])
      call organize_factor_loadings(loadings, organized, signs, info, factor_order=[2, 1], item_order=[3, 1, 2])
      call check(info == 0, 'fa.organize info', failures)
      call check(maxval(abs(signs - [1.0_dp, -1.0_dp])) < 1.0e-12_dp, 'fa.organize signs', failures)
      call check(maxval(abs(organized(:, 1) - [0.5_dp, 0.2_dp, 0.7_dp])) < 1.0e-12_dp, &
         'fa.organize factor reorder', failures)
      call check(maxval(abs(organized(:, 2) - [-0.1_dp, 0.8_dp, 0.4_dp])) < 1.0e-12_dp, &
         'fa.organize item reorder', failures)
   end subroutine test_organize_factor_loadings

   subroutine test_parallel_select(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: data(180, 6), z1, z2, t
      integer :: keys(6, 2)
      integer, allocatable :: nf(:), nc(:), status(:), seed(:)
      integer :: i, nseed

      do i = 1, 180
         t = real(i, dp)
         z1 = sin(0.071_dp * t) + 0.30_dp * cos(0.023_dp * t)
         z2 = cos(0.083_dp * t) + 0.25_dp * sin(0.031_dp * t)
         data(i, 1) = z1 + 0.08_dp * sin(0.37_dp * t)
         data(i, 2) = 0.92_dp * z1 + 0.08_dp * cos(0.29_dp * t)
         data(i, 3) = 0.88_dp * z1 + 0.07_dp * sin(0.43_dp * t)
         data(i, 4) = z2 + 0.08_dp * cos(0.41_dp * t)
         data(i, 5) = 0.90_dp * z2 + 0.07_dp * sin(0.33_dp * t)
         data(i, 6) = 0.86_dp * z2 + 0.08_dp * cos(0.27_dp * t)
      end do
      keys = 0
      keys(1:3, 1) = 1
      keys(4:6, 2) = 1
      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed = 24680
      call random_seed(put=seed)
      call parallel_select_by_keys(data, keys, 40, 0.95_dp, nf, nc, status)
      call check(all(status == 0), 'paSelect status', failures)
      call check(all(nf == 1), 'paSelect factor counts', failures)
      call check(all(nc == 1), 'paSelect component counts', failures)
   end subroutine test_parallel_select

   subroutine test_vss_select(failures)
      integer, intent(inout) :: failures !! Running assertion-failure count.
      real(dp) :: data(180, 6), z1, z2, t
      integer :: keys(6, 2)
      integer, allocatable :: v1(:), v2(:), status(:)
      integer :: i

      do i = 1, 180
         t = real(i, dp)
         z1 = sin(0.071_dp * t) + 0.30_dp * cos(0.023_dp * t)
         z2 = cos(0.083_dp * t) + 0.25_dp * sin(0.031_dp * t)
         data(i, 1) = z1 + 0.08_dp * sin(0.37_dp * t)
         data(i, 2) = 0.92_dp * z1 + 0.08_dp * cos(0.29_dp * t)
         data(i, 3) = 0.88_dp * z1 + 0.07_dp * sin(0.43_dp * t)
         data(i, 4) = z2 + 0.08_dp * cos(0.41_dp * t)
         data(i, 5) = 0.90_dp * z2 + 0.07_dp * sin(0.33_dp * t)
         data(i, 6) = 0.86_dp * z2 + 0.08_dp * cos(0.27_dp * t)
      end do
      keys = 0
      keys(1:3, 1) = 1
      keys(4:6, 2) = 1
      call vss_select_by_keys(data, keys, 3, v1, v2, status)
      call check(all(status == 0), 'vssSelect status', failures)
      call check(all(v1 >= 1 .and. v1 <= 2), 'vssSelect complexity-one range', failures)
      call check(all(v2 == 2), 'vssSelect complexity-two counts', failures)
   end subroutine test_vss_select

end program test_computational_extended
