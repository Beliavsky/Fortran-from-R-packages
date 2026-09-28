module optr_linear
   use optr_kinds, only : dp
   use optr_status, only : optr_success, optr_invalid_shape, optr_singular, &
      optr_not_positive_definite
   implicit none
   private

   public :: back_substitution
   public :: cholesky_decompose
   public :: cholesky_solve
   public :: forward_substitution
   public :: gaussian_eliminate
   public :: gaussian_solve
   public :: invert_matrix
   public :: lu_decompose
   public :: lu_solve
   public :: lu_split
   public :: matrix_reorder
   public :: multiply_factor
   public :: zero_small_values

   interface gaussian_solve
      module procedure gaussian_solve_vector
      module procedure gaussian_solve_matrix
   end interface gaussian_solve

   interface lu_solve
      module procedure lu_solve_vector
      module procedure lu_solve_matrix
   end interface lu_solve

   interface cholesky_solve
      module procedure cholesky_solve_vector
      module procedure cholesky_solve_matrix
   end interface cholesky_solve

contains

   pure subroutine matrix_reorder(a, a_ordered, permutation, status, tol)
      real(dp), intent(in) :: a(:, :) !! Matrix whose rows are reordered by scaled partial pivoting; must be square.
      real(dp), allocatable, intent(out) :: a_ordered(:, :) !! Row-reordered copy of `a` with the same square shape.
      integer, allocatable, intent(out) :: permutation(:) !! Row permutation such that `a_ordered = a(permutation,:)`.
      integer, intent(out) :: status !! Zero on success; nonzero when `a` is nonsquare or numerically singular.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-16_dp`.

      real(dp), allocatable :: scale(:)
      real(dp) :: pivot_tol
      real(dp) :: ratio
      real(dp) :: best_ratio
      real(dp), allocatable :: row_tmp(:)
      real(dp) :: scale_tmp
      integer :: i
      integer :: k
      integer :: n
      integer :: p
      integer :: perm_tmp

      status = optr_success
      if (size(a, 1) /= size(a, 2)) then
         status = optr_invalid_shape
         allocate(a_ordered(0, 0), permutation(0))
         return
      end if

      n = size(a, 1)
      allocate(a_ordered(n, n), permutation(n), scale(n), row_tmp(n))
      a_ordered = a
      permutation = [(i, i = 1, n)]
      pivot_tol = 1.0e-16_dp
      if (present(tol)) pivot_tol = max(tol, 0.0_dp)

      do i = 1, n
         scale(i) = maxval(abs(a_ordered(i, :)))
         if (scale(i) <= pivot_tol) then
            status = optr_singular
            return
         end if
      end do

      do k = 1, n - 1
         p = k
         best_ratio = -1.0_dp
         do i = k, n
            ratio = abs(a_ordered(i, k)) / scale(i)
            if (ratio > best_ratio) then
               best_ratio = ratio
               p = i
            end if
         end do

         if (abs(a_ordered(p, k)) <= pivot_tol) then
            status = optr_singular
            return
         end if

         if (p /= k) then
            row_tmp = a_ordered(k, :)
            a_ordered(k, :) = a_ordered(p, :)
            a_ordered(p, :) = row_tmp
            scale_tmp = scale(k)
            scale(k) = scale(p)
            scale(p) = scale_tmp
            perm_tmp = permutation(k)
            permutation(k) = permutation(p)
            permutation(p) = perm_tmp
         end if
      end do

      if (n > 0) then
         if (abs(a_ordered(n, n)) <= pivot_tol) status = optr_singular
      end if
   end subroutine matrix_reorder

   pure function multiply_factor(row_index, a, pivot_index) result(lambda)
      integer, intent(in) :: row_index !! One-based row index whose elimination multiplier is requested.
      real(dp), intent(in) :: a(:, :) !! Matrix containing the current pivot row and target row.
      integer, intent(in) :: pivot_index !! One-based pivot row and column index.
      real(dp) :: lambda

      if (abs(a(row_index, pivot_index)) <= tiny(1.0_dp)) then
         lambda = 0.0_dp
      else
         lambda = a(row_index, pivot_index) / a(pivot_index, pivot_index)
      end if
   end function multiply_factor

   pure subroutine gaussian_eliminate(a, b, u, c, permutation, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:, :) !! Right-hand-side matrix with one row per row of `a`.
      real(dp), allocatable, intent(out) :: u(:, :) !! Upper-triangular matrix after elimination with row pivoting.
      real(dp), allocatable, intent(out) :: c(:, :) !! Transformed right-hand-side matrix corresponding to `u`.
      integer, allocatable, intent(out) :: permutation(:) !! Applied row permutation of the original system.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular pivot.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp) :: factor
      real(dp) :: pivot_tol
      real(dp) :: best
      real(dp), allocatable :: row_tmp(:)
      real(dp), allocatable :: rhs_tmp(:)
      integer :: i
      integer :: k
      integer :: n
      integer :: p
      integer :: perm_tmp

      status = optr_success
      if (size(a, 1) /= size(a, 2) .or. size(b, 1) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(u(0, 0), c(0, 0), permutation(0))
         return
      end if

      n = size(a, 1)
      allocate(u(n, n), c(n, size(b, 2)), permutation(n), row_tmp(n), rhs_tmp(size(b, 2)))
      u = a
      c = b
      permutation = [(i, i = 1, n)]
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = max(tol, 0.0_dp)

      do k = 1, n - 1
         p = k
         best = abs(u(k, k))
         do i = k + 1, n
            if (abs(u(i, k)) > best) then
               best = abs(u(i, k))
               p = i
            end if
         end do
         if (best <= pivot_tol) then
            status = optr_singular
            return
         end if

         if (p /= k) then
            row_tmp = u(k, :)
            u(k, :) = u(p, :)
            u(p, :) = row_tmp
            rhs_tmp = c(k, :)
            c(k, :) = c(p, :)
            c(p, :) = rhs_tmp
            perm_tmp = permutation(k)
            permutation(k) = permutation(p)
            permutation(p) = perm_tmp
         end if

         do i = k + 1, n
            if (abs(u(i, k)) <= tiny(1.0_dp)) cycle
            factor = u(i, k) / u(k, k)
            u(i, k:n) = u(i, k:n) - factor * u(k, k:n)
            c(i, :) = c(i, :) - factor * c(k, :)
         end do
      end do

      if (n > 0) then
         if (abs(u(n, n)) <= pivot_tol) status = optr_singular
      end if
      call zero_small_values(u)
      call zero_small_values(c)
   end subroutine gaussian_eliminate

   pure subroutine gaussian_solve_matrix(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:, :) !! One or more right-hand-side columns.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, one column per right-hand side.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular system.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: c(:, :)
      integer, allocatable :: permutation(:)
      real(dp) :: pivot_tol
      integer :: j
      integer :: n
      real(dp), allocatable :: xj(:)

      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call gaussian_eliminate(a, b, u, c, permutation, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0, 0))
         return
      end if

      n = size(a, 1)
      allocate(x(n, size(b, 2)))
      do j = 1, size(b, 2)
         call back_substitution(u, c(:, j), xj, status, pivot_tol)
         if (status /= optr_success) return
         x(:, j) = xj
      end do
   end subroutine gaussian_solve_matrix

   pure subroutine gaussian_solve_vector(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector with length equal to the order of `a`.
      real(dp), allocatable, intent(out) :: x(:) !! Solution vector.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular system.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: bm(:, :)
      real(dp), allocatable :: xm(:, :)
      real(dp) :: pivot_tol

      if (size(b) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if
      allocate(bm(size(b), 1))
      bm(:, 1) = b
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call gaussian_solve_matrix(a, bm, xm, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0))
         return
      end if
      allocate(x(size(b)))
      x = xm(:, 1)
   end subroutine gaussian_solve_vector

   pure subroutine lu_decompose(a, lu, permutation, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square matrix to factor as a row-permuted `L U` product.
      real(dp), allocatable, intent(out) :: lu(:, :) !! Compact factorization with unit-diagonal `L` below `U`.
      integer, allocatable, intent(out) :: permutation(:) !! Row permutation such that `a(permutation,:) = L U`.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular pivot.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: row_tmp(:)
      real(dp) :: best_ratio
      real(dp) :: pivot_tol
      real(dp) :: ratio
      real(dp) :: scale_tmp
      integer :: i
      integer :: k
      integer :: n
      integer :: p
      integer :: perm_tmp

      status = optr_success
      if (size(a, 1) /= size(a, 2)) then
         status = optr_invalid_shape
         allocate(lu(0, 0), permutation(0))
         return
      end if

      n = size(a, 1)
      allocate(lu(n, n), permutation(n), scale(n), row_tmp(n))
      lu = a
      permutation = [(i, i = 1, n)]
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = max(tol, 0.0_dp)

      do i = 1, n
         scale(i) = maxval(abs(lu(i, :)))
         if (scale(i) <= pivot_tol) then
            status = optr_singular
            return
         end if
      end do

      do k = 1, n - 1
         p = k
         best_ratio = -1.0_dp
         do i = k, n
            ratio = abs(lu(i, k)) / scale(i)
            if (ratio > best_ratio) then
               best_ratio = ratio
               p = i
            end if
         end do

         if (abs(lu(p, k)) <= pivot_tol) then
            status = optr_singular
            return
         end if

         if (p /= k) then
            row_tmp = lu(k, :)
            lu(k, :) = lu(p, :)
            lu(p, :) = row_tmp
            scale_tmp = scale(k)
            scale(k) = scale(p)
            scale(p) = scale_tmp
            perm_tmp = permutation(k)
            permutation(k) = permutation(p)
            permutation(p) = perm_tmp
         end if

         do i = k + 1, n
            lu(i, k) = lu(i, k) / lu(k, k)
            lu(i, k + 1:n) = lu(i, k + 1:n) - lu(i, k) * lu(k, k + 1:n)
         end do
      end do

      if (n > 0) then
         if (abs(lu(n, n)) <= pivot_tol) status = optr_singular
      end if
   end subroutine lu_decompose

   pure subroutine forward_substitution(lu, b, y, status)
      real(dp), intent(in) :: lu(:, :) !! Compact LU matrix whose strict lower triangle stores unit-diagonal `L`.
      real(dp), intent(in) :: b(:) !! Permuted right-hand-side vector.
      real(dp), allocatable, intent(out) :: y(:) !! Solution of `L y = b`.
      integer, intent(out) :: status !! Zero on success; nonzero when dimensions are incompatible.

      integer :: i
      integer :: n

      if (size(lu, 1) /= size(lu, 2) .or. size(b) /= size(lu, 1)) then
         status = optr_invalid_shape
         allocate(y(0))
         return
      end if

      n = size(b)
      allocate(y(n))
      y = b
      do i = 2, n
         y(i) = b(i) - dot_product(lu(i, 1:i - 1), y(1:i - 1))
      end do
      status = optr_success
   end subroutine forward_substitution

   pure subroutine back_substitution(u, c, x, status, tol)
      real(dp), intent(in) :: u(:, :) !! Upper-triangular matrix, or compact LU whose upper triangle is `U`.
      real(dp), intent(in) :: c(:) !! Right-hand-side vector in `U x = c`.
      real(dp), allocatable, intent(out) :: x(:) !! Back-substitution solution vector.
      integer, intent(out) :: status !! Zero on success; nonzero for dimensions or a near-zero diagonal pivot.
      real(dp), optional, intent(in) :: tol !! Absolute diagonal tolerance; defaults to `1.0e-7_dp`.

      real(dp) :: pivot_tol
      integer :: i
      integer :: n

      if (size(u, 1) /= size(u, 2) .or. size(c) /= size(u, 1)) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if

      n = size(c)
      allocate(x(n))
      x = 0.0_dp
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = max(tol, 0.0_dp)
      do i = n, 1, -1
         if (abs(u(i, i)) <= pivot_tol) then
            status = optr_singular
            return
         end if
         if (i < n) then
            x(i) = (c(i) - dot_product(u(i, i + 1:n), x(i + 1:n))) / u(i, i)
         else
            x(i) = c(i) / u(i, i)
         end if
      end do
      status = optr_success
   end subroutine back_substitution

   pure subroutine lu_solve_matrix(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:, :) !! One or more right-hand-side columns.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, one column per right-hand side.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular system.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: lu(:, :)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: xj(:)
      integer, allocatable :: permutation(:)
      real(dp) :: pivot_tol
      integer :: j
      integer :: n

      if (size(b, 1) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0, 0))
         return
      end if

      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call lu_decompose(a, lu, permutation, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0, 0))
         return
      end if

      n = size(a, 1)
      allocate(x(n, size(b, 2)))
      do j = 1, size(b, 2)
         call forward_substitution(lu, b(permutation, j), y, status)
         if (status /= optr_success) return
         call back_substitution(lu, y, xj, status, pivot_tol)
         if (status /= optr_success) return
         x(:, j) = xj
      end do
   end subroutine lu_solve_matrix

   pure subroutine lu_solve_vector(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix in `A x = b`.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), allocatable, intent(out) :: x(:) !! Solution vector.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or a singular system.
      real(dp), optional, intent(in) :: tol !! Absolute pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: bm(:, :)
      real(dp), allocatable :: xm(:, :)
      real(dp) :: pivot_tol

      if (size(b) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if
      allocate(bm(size(b), 1))
      bm(:, 1) = b
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call lu_solve_matrix(a, bm, xm, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0))
         return
      end if
      allocate(x(size(b)))
      x = xm(:, 1)
   end subroutine lu_solve_vector

   pure subroutine lu_split(lu, l, u, status)
      real(dp), intent(in) :: lu(:, :) !! Compact square LU factorization with `L` below and `U` on/above the diagonal.
      real(dp), allocatable, intent(out) :: l(:, :) !! Explicit unit-lower-triangular factor.
      real(dp), allocatable, intent(out) :: u(:, :) !! Explicit upper-triangular factor.
      integer, intent(out) :: status !! Zero on success; nonzero when `lu` is not square.

      integer :: i
      integer :: j
      integer :: n

      if (size(lu, 1) /= size(lu, 2)) then
         status = optr_invalid_shape
         allocate(l(0, 0), u(0, 0))
         return
      end if

      n = size(lu, 1)
      allocate(l(n, n), u(n, n))
      l = 0.0_dp
      u = 0.0_dp
      do i = 1, n
         l(i, i) = 1.0_dp
         do j = 1, i - 1
            l(i, j) = lu(i, j)
         end do
         do j = i, n
            u(i, j) = lu(i, j)
         end do
      end do
      status = optr_success
   end subroutine lu_split

   pure subroutine cholesky_decompose(a, l, status, tol)
      real(dp), intent(in) :: a(:, :) !! Real symmetric matrix to factor as `L transpose(L)`.
      real(dp), allocatable, intent(out) :: l(:, :) !! Lower-triangular Cholesky factor.
      integer, intent(out) :: status !! Zero on success; nonzero if shape is invalid or the matrix is not positive definite.
      real(dp), optional, intent(in) :: tol !! Nonnegative definiteness tolerance; defaults to `1.0e-7_dp`.

      real(dp) :: diagonal
      real(dp) :: pivot_tol
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      if (size(a, 1) /= size(a, 2)) then
         status = optr_invalid_shape
         allocate(l(0, 0))
         return
      end if

      n = size(a, 1)
      allocate(l(n, n))
      l = 0.0_dp
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = max(tol, 0.0_dp)

      do i = 1, n
         diagonal = a(i, i)
         do k = 1, i - 1
            diagonal = diagonal - l(i, k) * l(i, k)
         end do
         if (diagonal <= pivot_tol) then
            status = optr_not_positive_definite
            return
         end if
         l(i, i) = sqrt(diagonal)

         do j = i + 1, n
            l(j, i) = a(j, i)
            do k = 1, i - 1
               l(j, i) = l(j, i) - l(j, k) * l(i, k)
            end do
            l(j, i) = l(j, i) / l(i, i)
         end do
      end do
      status = optr_success
   end subroutine cholesky_decompose

   pure subroutine cholesky_solve_matrix(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! One or more right-hand-side columns.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix.
      integer, intent(out) :: status !! Zero on success; nonzero for shape or positive-definiteness failures.
      real(dp), optional, intent(in) :: tol !! Positive-definiteness tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: y(:)
      real(dp) :: pivot_tol
      integer :: i
      integer :: j
      integer :: n

      if (size(b, 1) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0, 0))
         return
      end if

      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call cholesky_decompose(a, l, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0, 0))
         return
      end if

      n = size(a, 1)
      allocate(x(n, size(b, 2)), y(n))
      do j = 1, size(b, 2)
         y(1) = b(1, j) / l(1, 1)
         do i = 2, n
            y(i) = (b(i, j) - dot_product(l(i, 1:i - 1), y(1:i - 1))) / l(i, i)
         end do
         x(n, j) = y(n) / l(n, n)
         do i = n - 1, 1, -1
            x(i, j) = (y(i) - dot_product(l(i + 1:n, i), x(i + 1:n, j))) / l(i, i)
         end do
      end do
      status = optr_success
   end subroutine cholesky_solve_matrix

   pure subroutine cholesky_solve_vector(a, b, x, status, tol)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), allocatable, intent(out) :: x(:) !! Solution vector.
      integer, intent(out) :: status !! Zero on success; nonzero for shape or positive-definiteness failures.
      real(dp), optional, intent(in) :: tol !! Positive-definiteness tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: bm(:, :)
      real(dp), allocatable :: xm(:, :)
      real(dp) :: pivot_tol

      if (size(b) /= size(a, 1)) then
         status = optr_invalid_shape
         allocate(x(0))
         return
      end if
      allocate(bm(size(b), 1))
      bm(:, 1) = b
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call cholesky_solve_matrix(a, bm, xm, status, pivot_tol)
      if (status /= optr_success) then
         allocate(x(0))
         return
      end if
      allocate(x(size(b)))
      x = xm(:, 1)
   end subroutine cholesky_solve_vector

   pure subroutine invert_matrix(a, a_inverse, status, tol)
      real(dp), intent(in) :: a(:, :) !! Square nonsingular matrix to invert.
      real(dp), allocatable, intent(out) :: a_inverse(:, :) !! Matrix inverse computed through one LU factorization.
      integer, intent(out) :: status !! Zero on success; nonzero for shape errors or singularity.
      real(dp), optional, intent(in) :: tol !! Absolute LU pivot tolerance; defaults to `1.0e-7_dp`.

      real(dp), allocatable :: identity(:, :)
      real(dp) :: pivot_tol
      integer :: i
      integer :: n

      if (size(a, 1) /= size(a, 2)) then
         status = optr_invalid_shape
         allocate(a_inverse(0, 0))
         return
      end if

      n = size(a, 1)
      allocate(identity(n, n))
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      pivot_tol = 1.0e-7_dp
      if (present(tol)) pivot_tol = tol
      call lu_solve_matrix(a, identity, a_inverse, status, pivot_tol)
      if (status == optr_success) call zero_small_values(a_inverse)
   end subroutine invert_matrix

   pure subroutine zero_small_values(a, threshold)
      real(dp), intent(inout) :: a(:, :) !! Matrix whose tiny finite values are replaced by exact zero.
      real(dp), optional, intent(in) :: threshold !! Nonnegative cutoff; defaults to machine epsilon for `dp`.

      real(dp) :: cutoff

      cutoff = epsilon(1.0_dp)
      if (present(threshold)) cutoff = max(threshold, 0.0_dp)
      where (abs(a) < cutoff) a = 0.0_dp
   end subroutine zero_small_values

end module optr_linear
