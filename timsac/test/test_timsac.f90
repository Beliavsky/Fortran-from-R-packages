! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program test_timsac
  use timsac, only: dp, autocorrelation, autocorrelation_result, matrix_filter, white_noise
  implicit none

  call test_autocorrelation()
  call test_matrix_filter()
  call test_white_noise()
  print '(a)', 'All tests passed.'

contains

  subroutine test_autocorrelation()
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: x(5)
    type(autocorrelation_result) :: result

    x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    result = autocorrelation(x, 2)

    call assert_close(result%mean, 3.0_dp, tol, 'autocorrelation mean')
    call assert_close(result%covariance(0), 2.0_dp, tol, 'covariance lag 0')
    call assert_close(result%covariance(1), 0.8_dp, tol, 'covariance lag 1')
    call assert_close(result%covariance(2), -0.2_dp, tol, 'covariance lag 2')
    call assert_close(result%correlation(0), 1.0_dp, tol, 'correlation lag 0')
    call assert_close(result%correlation(1), 0.4_dp, tol, 'correlation lag 1')
    call assert_close(result%correlation(2), -0.1_dp, tol, 'correlation lag 2')
  end subroutine test_autocorrelation


  subroutine test_matrix_filter()
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: x(3,1), coefficients(1,1,1), initial(1,1)
    real(dp), allocatable :: y(:,:)

    x(:,1) = [1.0_dp, 2.0_dp, 3.0_dp]
    coefficients(1,1,1) = 0.5_dp
    initial = 0.0_dp

    y = matrix_filter(x, coefficients, recursive=.false., initial=initial)
    call assert_close(y(1,1), 1.0_dp, tol, 'convolution t=1')
    call assert_close(y(2,1), 1.5_dp, tol, 'convolution t=2')
    call assert_close(y(3,1), 2.0_dp, tol, 'convolution t=3')

    y = matrix_filter(x, coefficients, recursive=.true., initial=initial)
    call assert_close(y(1,1), 1.0_dp, tol, 'recursive t=1')
    call assert_close(y(2,1), 2.5_dp, tol, 'recursive t=2')
    call assert_close(y(3,1), 4.25_dp, tol, 'recursive t=3')
  end subroutine test_matrix_filter


  subroutine test_white_noise()
    real(dp) :: covariance(2,2)
    real(dp), allocatable :: noise(:,:)

    covariance = 0.0_dp
    covariance(1,1) = 1.0_dp
    covariance(2,2) = 1.0_dp
    call white_noise(covariance, 100, noise)

    if (any(shape(noise) /= [2, 100])) error stop "white-noise shape"
    if (any(noise /= noise)) error stop "white-noise NaN"
  end subroutine test_white_noise


  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      print '(a)', 'FAILED: ' // label
      print '(a,es24.16)', 'actual   = ', actual
      print '(a,es24.16)', 'expected = ', expected
      error stop 1
    end if
  end subroutine assert_close

end program test_timsac
