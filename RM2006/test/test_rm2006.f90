! SPDX-License-Identifier: GPL-2.0-or-later
program test_rm2006
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use rm2006_kinds, only : dp
   use rm2006_module, only : rm2006, rm2006_covariance, rm2006_scale_weights, &
      rm2006_success, rm2006_bad_parameter, rm2006_nonfinite_data
   implicit none

   integer :: failures

   failures = 0
   call test_reference_values(failures)
   call test_constant_series(failures)
   call test_weights(failures)
   call test_small_sample(failures)
   call test_validation(failures)

   if (failures > 0) then
      print '(a,i0)', 'RM2006 tests failed: ', failures
      error stop 1
   end if
   print '(a)', 'All RM2006 tests passed.'

contains

   subroutine test_reference_values(failures)
      integer, intent(inout) :: failures
      real(dp) :: data(5, 2)
      real(dp) :: expected(2, 2, 6)
      real(dp), allocatable :: actual(:, :, :)
      integer :: status

      data = reshape([ &
          0.01_dp, -0.02_dp, &
          0.03_dp,  0.01_dp, &
         -0.02_dp,  0.04_dp, &
          0.00_dp, -0.01_dp, &
          0.05_dp,  0.02_dp  &
         ], shape(data), order=[2, 1])

      expected(:, :, 1) = reshape([ &
          5.4754524214879947e-4_dp, -5.0891775335741456e-5_dp, &
         -5.0891775335741456e-5_dp,  4.9647564926243637e-4_dp], [2, 2])
      expected(:, :, 2) = reshape([ &
          4.1514016426550965e-4_dp, -9.4673226443340295e-5_dp, &
         -9.4673226443340295e-5_dp,  4.6791399973300220e-4_dp], [2, 2])
      expected(:, :, 3) = reshape([ &
          5.6942100257408989e-4_dp,  2.8281752953193754e-5_dp, &
          2.8281752953193754e-5_dp,  3.5663316495997886e-4_dp], [2, 2])
      expected(:, :, 4) = reshape([ &
          5.1876527272609022e-4_dp, -2.2598810191071250e-4_dp, &
         -2.2598810191071250e-4_dp,  7.4016900959751174e-4_dp], [2, 2])
      expected(:, :, 5) = reshape([ &
          3.6198363676911966e-4_dp, -1.5229636781858105e-4_dp, &
         -1.5229636781858105e-4_dp,  5.3880817644707114e-4_dp], [2, 2])
      expected(:, :, 6) = reshape([ &
          1.0206166024874037e-3_dp,  2.0199121079065137e-4_dp, &
          2.0199121079065137e-4_dp,  4.9647564926243626e-4_dp], [2, 2])

      call rm2006(data, actual, tau0=20.0_dp, tau1=2.0_dp, &
                   kmax=3, rho=1.5_dp, status=status)
      call assert_true(status == rm2006_success, 'reference status', failures)
      call assert_true(all(shape(actual) == [2, 2, 6]), &
                       'reference shape', failures)
      call assert_close(maxval(abs(actual - expected)), 0.0_dp, 5.0e-16_dp, &
                        'reference values', failures)
   end subroutine test_reference_values


   subroutine test_constant_series(failures)
      integer, intent(inout) :: failures
      real(dp) :: data(20, 2)
      real(dp) :: target(2, 2)
      real(dp), allocatable :: actual(:, :, :)
      integer :: t

      data(:, 1) = 0.02_dp
      data(:, 2) = -0.03_dp
      target = reshape([0.0004_dp, -0.0006_dp, -0.0006_dp, 0.0009_dp], [2, 2])
      call rm2006_covariance(data, actual)
      do t = 1, size(actual, 3)
         call assert_close(maxval(abs(actual(:, :, t) - target)), 0.0_dp, &
                           5.0e-17_dp, 'constant series', failures)
      end do
   end subroutine test_constant_series


   subroutine test_weights(failures)
      integer, intent(inout) :: failures
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: weights(:)
      integer :: status

      call rm2006_scale_weights(tau, weights, tau0=20.0_dp, tau1=2.0_dp, &
                                kmax=3, rho=1.5_dp, status=status)
      call assert_true(status == rm2006_success, 'weight status', failures)
      call assert_close(maxval(abs(tau - [2.0_dp, 3.0_dp, 4.5_dp])), &
                        0.0_dp, 1.0e-15_dp, 'time scales', failures)
      call assert_close(sum(weights), 1.0_dp, 2.0e-15_dp, &
                        'normalized weights', failures)
      call assert_close(maxval(abs(weights - [ &
          0.40457555159723768_dp, 0.33333333333333331_dp, &
          0.26209111506942900_dp])), 0.0_dp, 5.0e-15_dp, &
          'scale weights', failures)
   end subroutine test_weights


   subroutine test_small_sample(failures)
      integer, intent(inout) :: failures
      real(dp) :: data(2, 1)
      real(dp), allocatable :: actual(:, :, :)
      integer :: status

      data(:, 1) = [0.01_dp, -0.02_dp]
      call rm2006_covariance(data, actual, kmax=14, status=status)
      call assert_true(status == rm2006_success, 'small sample status', failures)
      call assert_true(all(shape(actual) == [1, 1, 3]), &
                       'small sample shape', failures)
      call assert_true(all(actual >= 0.0_dp), 'small sample covariance', failures)
   end subroutine test_small_sample


   subroutine test_validation(failures)
      integer, intent(inout) :: failures
      real(dp) :: data(2, 1)
      real(dp), allocatable :: actual(:, :, :)
      integer :: status

      data(:, 1) = [0.01_dp, 0.02_dp]
      call rm2006_covariance(data, actual, tau0=1.0_dp, status=status)
      call assert_true(status == rm2006_bad_parameter, &
                       'bad parameter status', failures)

      data(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      call rm2006_covariance(data, actual, status=status)
      call assert_true(status == rm2006_nonfinite_data, &
                       'nonfinite status', failures)
   end subroutine test_validation


   subroutine assert_true(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(name)
         failures = failures + 1
      end if
   end subroutine assert_true


   subroutine assert_close(actual, expected, tolerance, name, failures)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures

      if (abs(actual - expected) > tolerance) then
         print '(a,2(es24.16,1x))', 'FAIL: ' // trim(name) // ': ', &
                                  actual, expected
         failures = failures + 1
      end if
   end subroutine assert_close

end program test_rm2006
