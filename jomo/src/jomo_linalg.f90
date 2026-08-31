! Small dense linear-algebra kernels used by jomo.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_linalg
   use jomo_kinds, only : dp
   implicit none
   private

   public :: chol_lower
   public :: is_spd
   public :: inverse_spd
   public :: solve_spd
   public :: logdet_spd
   public :: symmetrize
   public :: crossprod
   public :: quadratic_form

   interface solve_spd
      module procedure solve_spd_vector
      module procedure solve_spd_matrix
   end interface solve_spd

contains

   pure subroutine chol_lower(a, l, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix to factor, with shape n by n.
      real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor satisfying a = l*transpose(l) when info is zero.
      integer, intent(out) :: info !! Zero on success; otherwise the first leading principal minor that is not positive definite.
      integer :: i
      integer :: j
      real(dp) :: s

      if (size(a, 1) /= size(a, 2)) error stop "chol_lower: matrix must be square"
      if (any(shape(l) /= shape(a))) error stop "chol_lower: output shape mismatch"

      l = 0.0_dp
      info = 0
      do i = 1, size(a, 1)
         do j = 1, i
            if (j > 1) then
               s = dot_product(l(i, 1:j - 1), l(j, 1:j - 1))
            else
               s = 0.0_dp
            end if
            if (i == j) then
               s = a(i, i) - s
               if (s <= 0.0_dp) then
                  info = i
                  return
               end if
               l(i, j) = sqrt(s)
            else
               l(i, j) = (a(i, j) - s) / l(j, j)
            end if
         end do
      end do
   end subroutine chol_lower

   pure logical function is_spd(a)
      real(dp), intent(in) :: a(:, :) !! Candidate symmetric matrix to test for positive definiteness.
      real(dp), allocatable :: l(:, :)
      integer :: info

      if (size(a, 1) /= size(a, 2)) then
         is_spd = .false.
         return
      end if
      allocate(l(size(a, 1), size(a, 2)))
      call chol_lower(a, l, info)
      is_spd = info == 0
   end function is_spd

   pure subroutine solve_spd_vector(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix, shape n by n.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector of length n.
      real(dp), intent(out) :: x(:) !! Solution vector of length n.
      integer, intent(out) :: info !! Zero on success; positive if the Cholesky factorization fails.
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: y(:)
      integer :: i

      if (size(a, 1) /= size(a, 2)) error stop "solve_spd: matrix must be square"
      if (size(b) /= size(a, 1) .or. size(x) /= size(b)) error stop "solve_spd: shape mismatch"
      allocate(l(size(a, 1), size(a, 2)), y(size(b)))
      call chol_lower(a, l, info)
      if (info /= 0) then
         x = 0.0_dp
         return
      end if

      do i = 1, size(b)
         if (i > 1) then
            y(i) = (b(i) - dot_product(l(i, 1:i - 1), y(1:i - 1))) / l(i, i)
         else
            y(i) = b(i) / l(i, i)
         end if
      end do
      do i = size(b), 1, -1
         if (i < size(b)) then
            x(i) = (y(i) - dot_product(l(i + 1:size(b), i), x(i + 1:size(b)))) / l(i, i)
         else
            x(i) = y(i) / l(i, i)
         end if
      end do
   end subroutine solve_spd_vector

   pure subroutine solve_spd_matrix(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix, shape n by n.
      real(dp), intent(in) :: b(:, :) !! Right-hand-side matrix with n rows.
      real(dp), intent(out) :: x(:, :) !! Solution matrix with the same shape as b.
      integer, intent(out) :: info !! Zero on success; positive if the Cholesky factorization fails.
      integer :: j
      integer :: col_info

      if (size(b, 1) /= size(a, 1)) error stop "solve_spd: shape mismatch"
      if (any(shape(x) /= shape(b))) error stop "solve_spd: output shape mismatch"
      info = 0
      do j = 1, size(b, 2)
         call solve_spd_vector(a, b(:, j), x(:, j), col_info)
         if (col_info /= 0) then
            info = col_info
            return
         end if
      end do
   end subroutine solve_spd_matrix

   pure subroutine inverse_spd(a, ainv, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix to invert.
      real(dp), intent(out) :: ainv(:, :) !! Symmetric inverse of a.
      integer, intent(out) :: info !! Zero on success; positive if the Cholesky factorization fails.
      real(dp), allocatable :: eye(:, :)
      integer :: i

      if (size(a, 1) /= size(a, 2)) error stop "inverse_spd: matrix must be square"
      if (any(shape(ainv) /= shape(a))) error stop "inverse_spd: output shape mismatch"
      allocate(eye(size(a, 1), size(a, 2)))
      eye = 0.0_dp
      do i = 1, size(a, 1)
         eye(i, i) = 1.0_dp
      end do
      call solve_spd_matrix(a, eye, ainv, info)
      if (info == 0) call symmetrize(ainv)
   end subroutine inverse_spd

   pure function logdet_spd(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix whose log determinant is requested.
      real(dp) :: value
      real(dp), allocatable :: l(:, :)
      integer :: i
      integer :: info

      allocate(l(size(a, 1), size(a, 2)))
      call chol_lower(a, l, info)
      if (info /= 0) then
         value = -huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, size(a, 1)
         value = value + 2.0_dp * log(l(i, i))
      end do
   end function logdet_spd

   pure subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:, :) !! Square matrix replaced by the average of itself and its transpose.
      integer :: i
      integer :: j
      real(dp) :: v

      if (size(a, 1) /= size(a, 2)) error stop "symmetrize: matrix must be square"
      do j = 1, size(a, 2)
         do i = j + 1, size(a, 1)
            v = 0.5_dp * (a(i, j) + a(j, i))
            a(i, j) = v
            a(j, i) = v
         end do
      end do
   end subroutine symmetrize

   pure function crossprod(x) result(cp)
      real(dp), intent(in) :: x(:, :) !! Matrix whose column cross-product transpose(x)*x is requested.
      real(dp) :: cp(size(x, 2), size(x, 2))

      cp = matmul(transpose(x), x)
   end function crossprod

   pure function quadratic_form(x, ainv) result(value)
      real(dp), intent(in) :: x(:) !! Vector appearing on both sides of the quadratic form.
      real(dp), intent(in) :: ainv(:, :) !! Square matrix in the quadratic form, generally a precision matrix.
      real(dp) :: value

      if (size(ainv, 1) /= size(x) .or. size(ainv, 2) /= size(x)) error stop "quadratic_form: shape mismatch"
      value = dot_product(x, matmul(ainv, x))
   end function quadratic_form

end module jomo_linalg
