module rmpfr_matrix
  use rmpfr_types
  implicit none
  private

  public :: mpfr_matmul, mpfr_crossprod, mpfr_tcrossprod

contains

  subroutine mpfr_matmul(a, b, c)
    type(mpfr_real), intent(in) :: a(:, :) !! Left arbitrary-precision matrix with shape (m,k).
    type(mpfr_real), intent(in) :: b(:, :) !! Right arbitrary-precision matrix with shape (k,n).
    type(mpfr_real), allocatable, intent(out) :: c(:, :) !! Allocated product matrix with shape (m,n).
    type(mpfr_real) :: acc
    integer :: i, j, k, p

    if (size(a, 2) /= size(b, 1)) error stop "Rmpfr: nonconformable matrix multiplication"
    p = matrix_precision_pair(a, b)
    allocate(c(size(a, 1), size(b, 2)))
    do j = 1, size(b, 2)
      do i = 1, size(a, 1)
        acc = mpfr_zero(1, p)
        do k = 1, size(a, 2)
          acc = acc + mpfr_copy(a(i, k), p) * mpfr_copy(b(k, j), p)
        end do
        c(i, j) = acc
      end do
    end do
  end subroutine mpfr_matmul

  subroutine mpfr_crossprod(a, c, b)
    type(mpfr_real), intent(in) :: a(:, :) !! Left matrix whose transpose is used, shape (m,k).
    type(mpfr_real), allocatable, intent(out) :: c(:, :) !! Allocated result A^T*B or A^T*A.
    type(mpfr_real), intent(in), optional :: b(:, :) !! Optional right matrix with m rows; defaults to a.
    type(mpfr_real) :: acc
    integer :: i, j, k, p

    if (present(b)) then
      if (size(a, 1) /= size(b, 1)) error stop "Rmpfr: nonconformable crossprod"
      p = matrix_precision_pair(a, b)
      allocate(c(size(a, 2), size(b, 2)))
      do j = 1, size(b, 2)
        do i = 1, size(a, 2)
          acc = mpfr_zero(1, p)
          do k = 1, size(a, 1)
            acc = acc + mpfr_copy(a(k, i), p) * mpfr_copy(b(k, j), p)
          end do
          c(i, j) = acc
        end do
      end do
    else
      p = matrix_precision_one(a)
      allocate(c(size(a, 2), size(a, 2)))
      do j = 1, size(a, 2)
        do i = 1, size(a, 2)
          acc = mpfr_zero(1, p)
          do k = 1, size(a, 1)
            acc = acc + mpfr_copy(a(k, i), p) * mpfr_copy(a(k, j), p)
          end do
          c(i, j) = acc
        end do
      end do
    end if
  end subroutine mpfr_crossprod

  subroutine mpfr_tcrossprod(a, c, b)
    type(mpfr_real), intent(in) :: a(:, :) !! Left matrix, shape (m,k).
    type(mpfr_real), allocatable, intent(out) :: c(:, :) !! Allocated result A*B^T or A*A^T.
    type(mpfr_real), intent(in), optional :: b(:, :) !! Optional matrix with k columns; defaults to a.
    type(mpfr_real) :: acc
    integer :: i, j, k, p

    if (present(b)) then
      if (size(a, 2) /= size(b, 2)) error stop "Rmpfr: nonconformable tcrossprod"
      p = matrix_precision_pair(a, b)
      allocate(c(size(a, 1), size(b, 1)))
      do j = 1, size(b, 1)
        do i = 1, size(a, 1)
          acc = mpfr_zero(1, p)
          do k = 1, size(a, 2)
            acc = acc + mpfr_copy(a(i, k), p) * mpfr_copy(b(j, k), p)
          end do
          c(i, j) = acc
        end do
      end do
    else
      p = matrix_precision_one(a)
      allocate(c(size(a, 1), size(a, 1)))
      do j = 1, size(a, 1)
        do i = 1, size(a, 1)
          acc = mpfr_zero(1, p)
          do k = 1, size(a, 2)
            acc = acc + mpfr_copy(a(i, k), p) * mpfr_copy(a(j, k), p)
          end do
          c(i, j) = acc
        end do
      end do
    end if
  end subroutine mpfr_tcrossprod

  integer function matrix_precision_one(a) result(p)
    type(mpfr_real), intent(in) :: a(:, :) !! Matrix inspected for maximum element precision.
    integer :: i, j

    p = 2
    do j = 1, size(a, 2)
      do i = 1, size(a, 1)
        p = max(p, mpfr_precision(a(i, j)))
      end do
    end do
  end function matrix_precision_one

  integer function matrix_precision_pair(a, b) result(p)
    type(mpfr_real), intent(in) :: a(:, :) !! First matrix inspected for maximum element precision.
    type(mpfr_real), intent(in) :: b(:, :) !! Second matrix inspected for maximum element precision.

    p = max(matrix_precision_one(a), matrix_precision_one(b))
  end function matrix_precision_pair

end module rmpfr_matrix
