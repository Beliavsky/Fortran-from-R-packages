! SPDX-License-Identifier: GPL-3.0-only
module test_support
  use blmodel, only : dp
  implicit none
  private
  public :: assert_true, assert_close, assert_vector_close, assert_matrix_close

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') 'assertion failed: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, atol, rtol, message)
    real(dp), intent(in) :: actual, expected, atol, rtol
    character(len=*), intent(in) :: message
    real(dp) :: tolerance

    tolerance = atol + rtol * abs(expected)
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,3(1x,es24.16))') 'mismatch: ' // trim(message), actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected, atol, rtol, message)
    real(dp), intent(in) :: actual(:), expected(:), atol, rtol
    character(len=*), intent(in) :: message
    integer :: i

    call assert_true(size(actual) == size(expected), trim(message) // ' size')
    do i = 1, size(actual)
      call assert_close(actual(i), expected(i), atol, rtol, trim(message) // ' element')
    end do
  end subroutine assert_vector_close

  subroutine assert_matrix_close(actual, expected, atol, rtol, message)
    real(dp), intent(in) :: actual(:,:), expected(:,:), atol, rtol
    character(len=*), intent(in) :: message
    integer :: i, j

    call assert_true(all(shape(actual) == shape(expected)), trim(message) // ' shape')
    do j = 1, size(actual, 2)
      do i = 1, size(actual, 1)
        call assert_close(actual(i, j), expected(i, j), atol, rtol, trim(message) // ' element')
      end do
    end do
  end subroutine assert_matrix_close

end module test_support
