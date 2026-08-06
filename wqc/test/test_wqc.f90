! SPDX-License-Identifier: GPL-3.0-only
program test_wqc
   use wqc, only : apply_quantile_correlation, dp, quantile_correlation, &
      quantile_correlation_analysis, wqc_multi_result, wqc_pair_result
   use waveslim, only : mra, mra_result
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   call test_quantile_correlation_values
   call test_pair_analysis
   call test_multi_series
   call test_invalid_inputs
   write(*, '(a)') 'test_wqc: PASS'

contains

   subroutine test_quantile_correlation_values
      real(dp), parameter :: x(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      real(dp), parameter :: y(5) = [5.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
      real(dp), parameter :: tau(3) = [0.25_dp, 0.5_dp, 0.75_dp]
      real(dp), parameter :: expected(3) = [ &
         0.29211869733608864_dp, 0.0_dp, -0.5842373946721773_dp]
      real(dp), allocatable :: rho(:)
      integer :: stat
      character(len=:), allocatable :: errmsg

      allocate(character(len=1) :: errmsg)
      errmsg = ''
      rho = quantile_correlation(x, y, tau, stat, errmsg)
      call assert_true(stat == 0, 'quantile correlation status: '//errmsg)
      call assert_close_array(rho, expected, 2.0e-14_dp, 'quantile correlation reference values')
   end subroutine test_quantile_correlation_values


   subroutine test_pair_analysis
      integer, parameter :: n = 128
      real(dp) :: x(n), y(n), tau(3)
      type(wqc_pair_result) :: fit1, fit2
      type(mra_result) :: mx, my
      real(dp) :: direct(3)
      integer :: i, stat
      character(len=:), allocatable :: errmsg

      do i = 1, n
         x(i) = sin(0.071_dp * real(i, dp)) + 0.25_dp * cos(0.31_dp * real(i, dp))
         y(i) = 0.65_dp * x(i) + 0.35_dp * sin(0.19_dp * real(i, dp) + 0.4_dp)
      end do
      tau = [0.1_dp, 0.5_dp, 0.9_dp]

      allocate(character(len=1) :: errmsg)
      errmsg = ''
      fit1 = quantile_correlation_analysis(x, y, tau, wf='la8', j_levels=3, &
         n_sim=80, seed=24680, stat=stat, errmsg=errmsg)
      call assert_true(stat == 0, 'pair analysis status: '//errmsg)
      call assert_true(all(shape(fit1%estimated_qc) == [3, 3]), 'estimated QC shape')
      call assert_true(all(shape(fit1%ci_lower) == [3, 3]), 'lower CI shape')
      call assert_true(all(shape(fit1%ci_upper) == [3, 3]), 'upper CI shape')
      call assert_true(all(ieee_is_finite(fit1%estimated_qc)), 'finite estimated QC')
      call assert_true(all(ieee_is_finite(fit1%ci_lower)), 'finite lower CI')
      call assert_true(all(ieee_is_finite(fit1%ci_upper)), 'finite upper CI')
      call assert_true(all(fit1%ci_lower <= fit1%ci_upper), 'ordered confidence intervals')

      mx = mra(x, wf='la8', j_levels=3, method='modwt', boundary='periodic')
      my = mra(y, wf='la8', j_levels=3, method='modwt', boundary='periodic')
      direct = reference_qc(mx%detail(2)%values, my%detail(2)%values, tau)
      call assert_close_array(fit1%estimated_qc(2, :), direct, 3.0e-14_dp, &
         'level-2 QC agrees with independent formula')

      fit2 = quantile_correlation_analysis(x, y, tau, wf='la8', j_levels=3, &
         n_sim=80, seed=24680, stat=stat, errmsg=errmsg)
      call assert_true(stat == 0, 'repeat analysis status: '//errmsg)
      call assert_close_matrix(fit2%estimated_qc, fit1%estimated_qc, 0.0_dp, 'deterministic estimates')
      call assert_close_matrix(fit2%ci_lower, fit1%ci_lower, 0.0_dp, 'deterministic lower CI')
      call assert_close_matrix(fit2%ci_upper, fit1%ci_upper, 0.0_dp, 'deterministic upper CI')
   end subroutine test_pair_analysis


   subroutine test_multi_series
      integer, parameter :: n = 96
      real(dp) :: data(n, 3), tau(2)
      character(len=8) :: names(2)
      type(wqc_multi_result) :: fit
      integer :: i, stat
      character(len=:), allocatable :: errmsg

      do i = 1, n
         data(i, 1) = sin(0.08_dp * real(i, dp))
         data(i, 2) = 0.7_dp * data(i, 1) + 0.2_dp * cos(0.27_dp * real(i, dp))
         data(i, 3) = cos(0.13_dp * real(i, dp) + 0.6_dp)
      end do
      tau = [0.25_dp, 0.75_dp]
      names = ['target_a', 'target_b']

      allocate(character(len=1) :: errmsg)
      errmsg = ''
      fit = apply_quantile_correlation(data, tau, j_levels=2, n_sim=30, seed=13579, &
         series_names=names, stat=stat, errmsg=errmsg)
      call assert_true(stat == 0, 'multi-series status: '//errmsg)
      call assert_true(size(fit%series) == 2, 'multi-series result count')
      call assert_true(fit%series(1)%name == 'target_a', 'first series name')
      call assert_true(fit%series(2)%name == 'target_b', 'second series name')
      call assert_true(all(shape(fit%series(1)%analysis%estimated_qc) == [2, 2]), &
         'multi-series analysis shape')
   end subroutine test_multi_series


   subroutine test_invalid_inputs
      real(dp) :: x(8), y(7), tau(1)
      real(dp), allocatable :: rho(:)
      integer :: stat
      character(len=:), allocatable :: errmsg

      x = 1.0_dp
      y = 2.0_dp
      tau = 0.5_dp
      allocate(character(len=1) :: errmsg)
      errmsg = ''
      rho = quantile_correlation(x, y, tau, stat, errmsg)
      call assert_true(stat /= 0, 'mismatched lengths rejected')
      call assert_true(size(rho) == 0, 'failed result is empty')
   end subroutine test_invalid_inputs


   function reference_qc(x, y, tau) result(rho)
      real(dp), intent(in) :: x(:), y(:), tau(:)
      real(dp) :: rho(size(tau))

      real(dp) :: q, sx, xbar
      real(dp), allocatable :: sorted_y(:)
      integer :: i, j, n

      n = size(x)
      xbar = sum(x) / real(n, dp)
      sx = sqrt(sum((x - xbar)**2) / real(n - 1, dp))
      sorted_y = y
      call insertion_sort(sorted_y)
      do j = 1, size(tau)
         q = type7_sorted(sorted_y, tau(j))
         rho(j) = 0.0_dp
         do i = 1, n
            if (y(i) < q) then
               rho(j) = rho(j) + (x(i) - xbar) * (tau(j) - 1.0_dp)
            else
               rho(j) = rho(j) + (x(i) - xbar) * tau(j)
            end if
         end do
         rho(j) = rho(j) / (real(n, dp) * sqrt(tau(j) * (1.0_dp - tau(j))) * sx)
      end do
   end function reference_qc


   pure function type7_sorted(x, p) result(q)
      real(dp), intent(in) :: x(:), p
      real(dp) :: q
      real(dp) :: fraction, h
      integer :: j

      h = 1.0_dp + real(size(x) - 1, dp) * p
      j = int(floor(h))
      fraction = h - real(j, dp)
      if (j >= size(x)) then
         q = x(size(x))
      else
         q = (1.0_dp - fraction) * x(j) + fraction * x(j + 1)
      end if
   end function type7_sorted


   subroutine insertion_sort(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine insertion_sort


   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write(*, '(a)') 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine assert_true


   subroutine assert_close_array(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label

      call assert_true(size(actual) == size(expected), trim(label)//' size')
      if (size(actual) > 0) then
         call assert_true(maxval(abs(actual - expected)) <= tolerance, label)
      end if
   end subroutine assert_close_array


   subroutine assert_close_matrix(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:, :), expected(:, :), tolerance
      character(len=*), intent(in) :: label

      call assert_true(all(shape(actual) == shape(expected)), trim(label)//' shape')
      if (size(actual) > 0) then
         call assert_true(maxval(abs(actual - expected)) <= tolerance, label)
      end if
   end subroutine assert_close_matrix

end program test_wqc
