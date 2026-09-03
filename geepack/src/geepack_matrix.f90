! Linear-algebra glue for geepack using the shared rfortran-linalg package.
! Upstream geepack license: GPL (>= 3). See LICENSE, NOTICE.md, PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_matrix
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix, solve_system
   implicit none
   private

   public :: inverse_checked, solve_checked, at_inv_b_matrix, at_inv_b_vector
   public :: outer_product, diagonal_matrix, max_abs

contains

   subroutine inverse_checked(a, ainv, status)
      real(dp), intent(in) :: a(:, :) !! Square matrix to invert.
      real(dp), intent(out) :: ainv(:, :) !! Inverse matrix, same shape as a.
      integer, intent(out) :: status !! Zero on success; nonzero from rfortran-linalg on failure.

      call inverse_matrix(a, ainv, status)
   end subroutine inverse_checked

   subroutine solve_checked(a, b, x, status)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand side vector.
      real(dp), intent(out) :: x(:) !! Solution vector.
      integer, intent(out) :: status !! Zero on success; nonzero from rfortran-linalg on failure.

      call solve_system(a, b, x, status)
   end subroutine solve_checked

   subroutine at_inv_b_matrix(a, b, c, result, status)
      real(dp), intent(in) :: a(:, :) !! Left design matrix with row dimension matching b.
      real(dp), intent(in) :: b(:, :) !! Square nonsingular working covariance matrix.
      real(dp), intent(in) :: c(:, :) !! Right design matrix with row dimension matching b.
      real(dp), intent(out) :: result(:, :) !! Product transpose(a) * inverse(b) * c.
      integer, intent(out) :: status !! Zero on success; nonzero if inversion fails.
      real(dp), allocatable :: binv(:, :)

      allocate(binv(size(b, 1), size(b, 2)))
      call inverse_matrix(b, binv, status)
      if (status /= 0) then
         result = 0.0_dp
         return
      end if
      result = matmul(transpose(a), matmul(binv, c))
   end subroutine at_inv_b_matrix

   subroutine at_inv_b_vector(a, b, c, result, status)
      real(dp), intent(in) :: a(:, :) !! Left design matrix with row dimension matching b.
      real(dp), intent(in) :: b(:, :) !! Square nonsingular working covariance matrix.
      real(dp), intent(in) :: c(:) !! Right vector with length matching b.
      real(dp), intent(out) :: result(:) !! Product transpose(a) * inverse(b) * c.
      integer, intent(out) :: status !! Zero on success; nonzero if inversion fails.
      real(dp), allocatable :: binv(:, :)

      allocate(binv(size(b, 1), size(b, 2)))
      call inverse_matrix(b, binv, status)
      if (status /= 0) then
         result = 0.0_dp
         return
      end if
      result = matmul(transpose(a), matmul(binv, c))
   end subroutine at_inv_b_vector

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector.
      real(dp), intent(in) :: y(:) !! Right vector.
      real(dp) :: a(size(x), size(y))
      integer :: i

      do i = 1, size(x)
         a(i, :) = x(i) * y
      end do
   end function outer_product

   pure function diagonal_matrix(x) result(a)
      real(dp), intent(in) :: x(:) !! Diagonal entries.
      real(dp) :: a(size(x), size(x))
      integer :: i

      a = 0.0_dp
      do i = 1, size(x)
         a(i, i) = x(i)
      end do
   end function diagonal_matrix

   pure real(dp) function max_abs(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose largest absolute value is returned.

      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = maxval(abs(x))
      end if
   end function max_abs

end module geepack_matrix
