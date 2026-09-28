! SPDX-License-Identifier: GPL-3.0-or-later
! Direct linear solvers translated from sanic's Eigen-backed solver wrappers.
module sanic_direct
   use sanic_kinds, only: dp
   use sanic_sparse, only: sparse_matrix, sparse_is_symmetric
   use sanic_sparse_direct, only: sparse_chol_solve_core, sparse_lu_solve_core, sparse_qr_solve_core
   implicit none
   private

   public :: solve2_dense, solve2_sparse
   public :: solve_chol_dense, solve_chol_sparse
   public :: solve_lu_dense, solve_lu_sparse
   public :: solve_qr_dense, solve_qr_sparse
   public :: dense_is_symmetric

contains

   pure logical function dense_is_symmetric(a, tol) result(ok)
      !! Test dense symmetry with Eigen-style relative Frobenius norm scaling.
      real(dp), intent(in) :: a(:, :) !! Dense matrix to test; nonsquare input is never symmetric.
      real(dp), intent(in), optional :: tol !! Relative Frobenius tolerance; absent uses machine epsilon.
      real(dp) :: difference_norm, matrix_norm, threshold

      if (size(a, 1) /= size(a, 2)) then
         ok = .false.
         return
      end if
      threshold = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) threshold = tol
      end if
      matrix_norm = sqrt(sum(a*a))
      difference_norm = sqrt(sum((a - transpose(a))*(a - transpose(a))))
      ok = difference_norm <= threshold*matrix_norm
   end function dense_is_symmetric

   pure logical function selector_is_valid(value) result(ok)
      !! Check whether a sanic pivoting/ordering selector is in the public R range 0..2.
      integer, intent(in) :: value !! Integer selector to validate.

      ok = value >= 0 .and. value <= 2
   end function selector_is_valid

   pure subroutine solve2_dense(a, x, info, b, pivot, ordering)
      !! Solve Ax=b using QR for rectangular A, Cholesky for symmetric A, and LU otherwise.
      real(dp), intent(in) :: a(:, :) !! Dense coefficient matrix with m rows and n columns.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix; n rows and one column per RHS.
      integer, intent(out) :: info !! Status code; zero indicates a successful numerical solve.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS; if absent an identity RHS requests an inverse/pseudoinverse.
      integer, intent(in), optional :: pivot !! Pivot selector passed to the chosen direct solver; default is one.
      integer, intent(in), optional :: ordering !! Sparse-ordering compatibility argument; ignored for dense matrices.
      integer :: pivot_value

      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      pivot_value = 1
      if (present(pivot)) pivot_value = pivot

      if (size(a, 1) /= size(a, 2)) then
         if (present(b)) then
            call solve_qr_dense(a, x, info, b=b, pivot=pivot_value)
         else
            call solve_qr_dense(a, x, info, pivot=pivot_value)
         end if
      else if (dense_is_symmetric(a)) then
         if (present(b)) then
            call solve_chol_dense(a, x, info, b=b, pivot=pivot_value)
         else
            call solve_chol_dense(a, x, info, pivot=pivot_value)
         end if
      else
         if (present(b)) then
            call solve_lu_dense(a, x, info, b=b, pivot=pivot_value)
         else
            call solve_lu_dense(a, x, info, pivot=pivot_value)
         end if
      end if
   end subroutine solve2_dense

   pure subroutine solve2_sparse(a, x, info, b, pivot, ordering)
      !! Sparse-input counterpart of solve2 using sparse Cholesky/LU and dense fallback only for sparse QR.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix in package-local coordinate storage.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix; a%ncol rows and one column per RHS.
      integer, intent(out) :: info !! Status code; zero indicates a successful numerical solve.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS; if absent an identity RHS requests an inverse/pseudoinverse.
      integer, intent(in), optional :: pivot !! Pivot selector passed to Cholesky/QR; sparse LU ignores it as upstream does.
      integer, intent(in), optional :: ordering !! Sparse ordering selector; default follows the selected upstream solver.
      integer :: ordering_value, pivot_value

      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      pivot_value = 1
      if (present(pivot)) pivot_value = pivot

      if (a%nrow /= a%ncol) then
         ordering_value = 1
         if (present(ordering)) ordering_value = ordering
         if (present(b)) then
            call solve_qr_sparse(a, x, info, b=b, pivot=pivot_value, ordering=ordering_value)
         else
            call solve_qr_sparse(a, x, info, pivot=pivot_value, ordering=ordering_value)
         end if
      else if (sparse_is_symmetric(a)) then
         ordering_value = 0
         if (present(ordering)) ordering_value = ordering
         if (present(b)) then
            call solve_chol_sparse(a, x, info, b=b, pivot=pivot_value, ordering=ordering_value)
         else
            call solve_chol_sparse(a, x, info, pivot=pivot_value, ordering=ordering_value)
         end if
      else
         ordering_value = 1
         if (present(ordering)) ordering_value = ordering
         if (present(b)) then
            call solve_lu_sparse(a, x, info, b=b, pivot=pivot_value, ordering=ordering_value)
         else
            call solve_lu_sparse(a, x, info, pivot=pivot_value, ordering=ordering_value)
         end if
      end if
   end subroutine solve2_sparse

   pure subroutine solve_chol_dense(a, x, info, b, pivot, ordering)
      !! Solve a symmetric system by LLT (pivot=0) or pivoted LDLT (default).
      real(dp), intent(in) :: a(:, :) !! Square symmetric coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, or inverse when b is absent.
      integer, intent(out) :: info !! Status: zero success, one invalid shape, two failed factorization.
      real(dp), intent(in), optional :: b(:, :) !! RHS with size(a,1) rows; absent means the identity matrix.
      integer, intent(in), optional :: pivot !! Zero selects LLT; other supported values use LDLT as in upstream default.
      integer, intent(in), optional :: ordering !! Sparse-ordering compatibility argument; ignored for dense matrices.
      real(dp), allocatable :: rhs(:, :)
      integer :: mode, n

      n = size(a, 1)
      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (size(a, 2) /= n) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      if (present(ordering)) then
         ! Ordering has no meaning for dense input.
      end if
      call make_rhs(n, b, rhs, info)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if

      mode = 1
      if (present(pivot)) mode = pivot
      if (mode == 0) then
         call llt_solve(a, rhs, x, info)
      else
         call ldlt_solve(a, rhs, x, info)
      end if
   end subroutine solve_chol_dense

   pure subroutine solve_chol_sparse(a, x, info, b, pivot, ordering)
      !! Solve a sparse symmetric system by fill-permitting LLT or LDLT without dense matrix conversion.
      type(sparse_matrix), intent(in) :: a !! Sparse square symmetric coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, or inverse when b is absent.
      integer, intent(out) :: info !! Status: zero success, one invalid input, two failed factorization.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS with a%nrow rows; absent means identity.
      integer, intent(in), optional :: pivot !! Zero selects LLT; one/default selects LDLT; two follows upstream default fallback.
      integer, intent(in), optional :: ordering !! Zero/default uses AMD-like ordering; one falls back to zero; two is natural.
      real(dp), allocatable :: rhs(:, :)
      integer :: mode, ordering_value

      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      call make_rhs(a%nrow, b, rhs, info)
      if (info /= 0 .or. a%nrow /= a%ncol) then
         if (allocated(x)) deallocate(x)
         allocate(x(0, 0))
         if (info == 0) info = 1
         return
      end if
      mode = 1
      if (present(pivot)) mode = pivot
      ordering_value = 0
      if (present(ordering)) ordering_value = ordering
      call sparse_chol_solve_core(a, rhs, x, info, mode == 0, ordering_value)
   end subroutine solve_chol_sparse

   pure subroutine solve_lu_dense(a, x, info, b, pivot, ordering)
      !! Solve a square system by partial-pivot LU or full row/column pivoting.
      real(dp), intent(in) :: a(:, :) !! Square dense coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, or inverse when b is absent.
      integer, intent(out) :: info !! Status: zero success, one invalid shape, two singular factorization.
      real(dp), intent(in), optional :: b(:, :) !! RHS with size(a,1) rows; absent means the identity matrix.
      integer, intent(in), optional :: pivot !! Two requests full pivoting; zero/one use upstream partial-pivot fallback.
      integer, intent(in), optional :: ordering !! Sparse-ordering compatibility argument; ignored for dense matrices.
      real(dp), allocatable :: rhs(:, :)
      integer :: mode, n

      n = size(a, 1)
      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (size(a, 2) /= n) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      if (present(ordering)) then
         ! Ordering has no meaning for dense input.
      end if
      call make_rhs(n, b, rhs, info)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if

      mode = 1
      if (present(pivot)) mode = pivot
      if (mode == 2) then
         call gaussian_solve_full(a, rhs, x, info)
      else
         call gaussian_solve_partial(a, rhs, x, info)
      end if
   end subroutine solve_lu_dense

   pure subroutine solve_lu_sparse(a, x, info, b, pivot, ordering)
      !! Solve a sparse square system by fill-permitting sparse LU with partial row pivoting.
      type(sparse_matrix), intent(in) :: a !! Sparse square coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix, or inverse when b is absent.
      integer, intent(out) :: info !! Status: zero success, one invalid input, two singular factorization.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS with a%nrow rows; absent means identity.
      integer, intent(in), optional :: pivot !! Compatibility argument; sparse upstream ignores dense LU pivot selection.
      integer, intent(in), optional :: ordering !! 0 AMD-like, 1 COLAMD-like, 2 natural ordering.
      real(dp), allocatable :: rhs(:, :)
      integer :: ordering_value

      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      call make_rhs(a%nrow, b, rhs, info)
      if (info /= 0 .or. a%nrow /= a%ncol) then
         if (allocated(x)) deallocate(x)
         allocate(x(0, 0))
         if (info == 0) info = 1
         return
      end if
      ordering_value = 1
      if (present(ordering)) ordering_value = ordering
      call sparse_lu_solve_core(a, rhs, x, info, ordering_value)
   end subroutine solve_lu_sparse

   pure subroutine solve_qr_dense(a, x, info, b, pivot, ordering)
      !! Solve square/rectangular systems by Householder QR with Eigen-like basic solutions when underdetermined.
      real(dp), intent(in) :: a(:, :) !! Dense m-by-n coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Least-squares/basic solution with n rows.
      integer, intent(out) :: info !! Status: zero success, one invalid RHS, two rank deficiency/singularity.
      real(dp), intent(in), optional :: b(:, :) !! RHS with m rows; absent means an m-by-m identity matrix.
      integer, intent(in), optional :: pivot !! Zero disables column pivoting; one/two enable rank-revealing column pivoting.
      integer, intent(in), optional :: ordering !! Sparse-ordering compatibility argument; ignored for dense matrices.
      real(dp), allocatable :: rhs(:, :)
      integer :: m, mode, n

      m = size(a, 1)
      n = size(a, 2)
      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      call make_rhs(m, b, rhs, info)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if

      mode = 1
      if (present(pivot)) mode = pivot
      if (mode == 2) then
         call full_pivot_qr_least_squares(a, rhs, x, info)
      else
         call qr_least_squares(a, rhs, x, info, mode /= 0)
      end if
   end subroutine solve_qr_dense

   pure subroutine solve_qr_sparse(a, x, info, b, pivot, ordering)
      !! Solve a sparse rectangular system with sparse Givens QR and structural ordering.
      type(sparse_matrix), intent(in) :: a !! Sparse m-by-n coefficient matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Least-squares/basic solution in original column order.
      integer, intent(out) :: info !! Status: zero success, one invalid shape, two rank failure for underdetermined input.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS with a%nrow rows; absent means identity.
      integer, intent(in), optional :: pivot !! Compatibility argument; upstream sparse QR ignores dense pivot selection.
      integer, intent(in), optional :: ordering !! 0 dynamic symmetric degree, 1 dynamic column degree, 2 natural order.
      real(dp), allocatable :: rhs(:, :)
      integer :: ordering_value

      info = 0
      if (present(pivot)) then
         if (.not. selector_is_valid(pivot)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      if (present(ordering)) then
         if (.not. selector_is_valid(ordering)) then
            allocate(x(0, 0))
            info = 1
            return
         end if
      end if
      ordering_value = 1
      if (present(ordering)) ordering_value = ordering
      call make_rhs(a%nrow, b, rhs, info)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if
      call sparse_qr_solve_core(a, rhs, x, info, ordering_value)
   end subroutine solve_qr_sparse

   pure subroutine make_rhs(nrow, b, rhs, info)
      !! Form a validated RHS, using the identity matrix when b is absent.
      integer, intent(in) :: nrow !! Required number of RHS rows.
      real(dp), intent(in), optional :: b(:, :) !! Optional user RHS with nrow rows.
      real(dp), allocatable, intent(out) :: rhs(:, :) !! Validated RHS or identity matrix.
      integer, intent(out) :: info !! Status: zero success, one for row mismatch.
      integer :: i

      info = 0
      if (present(b)) then
         if (size(b, 1) /= nrow) then
            allocate(rhs(0, 0))
            info = 1
            return
         end if
         rhs = b
      else
         allocate(rhs(nrow, nrow))
         rhs = 0.0_dp
         do i = 1, nrow
            rhs(i, i) = 1.0_dp
         end do
      end if
   end subroutine make_rhs

   pure subroutine gaussian_solve_partial(a, b, x, info)
      !! Gaussian elimination with partial row pivoting for one or more RHS columns.
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with the same row count as a.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix.
      integer, intent(out) :: info !! Status: zero success, two if matrix is numerically singular.
      real(dp), allocatable :: u(:, :), rhs(:, :), tmp_row(:)
      real(dp) :: factor, scale, threshold
      integer :: i, k, n, p

      n = size(a, 1)
      u = a
      rhs = b
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)

      do k = 1, n - 1
         p = k - 1 + maxloc(abs(u(k:n, k)), dim=1)
         if (abs(u(p, k)) <= threshold) then
            info = 2
            allocate(x(0, 0))
            return
         end if
         if (p /= k) then
            tmp_row = u(k, :)
            u(k, :) = u(p, :)
            u(p, :) = tmp_row
            tmp_row = rhs(k, :)
            rhs(k, :) = rhs(p, :)
            rhs(p, :) = tmp_row
         end if
         do i = k + 1, n
            factor = u(i, k)/u(k, k)
            u(i, k:n) = u(i, k:n) - factor*u(k, k:n)
            rhs(i, :) = rhs(i, :) - factor*rhs(k, :)
         end do
      end do

      if (n > 0) then
         scale = abs(u(n, n))
      else
         scale = 1.0_dp
      end if
      if (n > 0 .and. scale <= threshold) then
         info = 2
         allocate(x(0, 0))
         return
      end if

      allocate(x(n, size(b, 2)))
      x = rhs
      do i = n, 1, -1
         if (abs(u(i, i)) <= threshold) then
            info = 2
            deallocate(x)
            allocate(x(0, 0))
            return
         end if
         if (i < n) x(i, :) = x(i, :) - matmul(u(i, i + 1:n), x(i + 1:n, :))
         x(i, :) = x(i, :)/u(i, i)
      end do
      info = 0
   end subroutine gaussian_solve_partial

   pure subroutine gaussian_solve_full(a, b, x, info)
      !! Gaussian elimination with full row/column pivoting for one or more RHS columns.
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with the same row count as a.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix in original variable order.
      integer, intent(out) :: info !! Status: zero success, two if matrix is numerically singular.
      real(dp), allocatable :: u(:, :), rhs(:, :), z(:, :), tmp_row(:), tmp_col(:)
      real(dp) :: factor, threshold
      integer, allocatable :: perm(:)
      integer :: i, j, k, loc(2), n, p, q, tmp

      n = size(a, 1)
      u = a
      rhs = b
      allocate(perm(n))
      perm = [(i, i=1, n)]
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)

      do k = 1, n - 1
         loc = maxloc(abs(u(k:n, k:n)))
         p = k + loc(1) - 1
         q = k + loc(2) - 1
         if (abs(u(p, q)) <= threshold) then
            info = 2
            allocate(x(0, 0))
            return
         end if
         if (p /= k) then
            tmp_row = u(k, :)
            u(k, :) = u(p, :)
            u(p, :) = tmp_row
            tmp_row = rhs(k, :)
            rhs(k, :) = rhs(p, :)
            rhs(p, :) = tmp_row
         end if
         if (q /= k) then
            tmp_col = u(:, k)
            u(:, k) = u(:, q)
            u(:, q) = tmp_col
            tmp = perm(k)
            perm(k) = perm(q)
            perm(q) = tmp
         end if
         do i = k + 1, n
            factor = u(i, k)/u(k, k)
            u(i, k:n) = u(i, k:n) - factor*u(k, k:n)
            rhs(i, :) = rhs(i, :) - factor*rhs(k, :)
         end do
      end do

      if (n > 0 .and. abs(u(n, n)) <= threshold) then
         info = 2
         allocate(x(0, 0))
         return
      end if
      allocate(z(n, size(b, 2)))
      z = rhs
      do i = n, 1, -1
         if (i < n) z(i, :) = z(i, :) - matmul(u(i, i + 1:n), z(i + 1:n, :))
         z(i, :) = z(i, :)/u(i, i)
      end do
      allocate(x(n, size(b, 2)))
      x = 0.0_dp
      do i = 1, n
         j = perm(i)
         x(j, :) = z(i, :)
      end do
      info = 0
   end subroutine gaussian_solve_full

   pure subroutine llt_solve(a, b, x, info)
      !! Cholesky LLT factorization and solve for a positive-definite matrix.
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with size(a,1) rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix.
      integer, intent(out) :: info !! Status: zero success, two if matrix is not numerically positive definite.
      real(dp), allocatable :: l(:, :), y(:, :)
      real(dp) :: diagonal, threshold
      integer :: i, j, n

      n = size(a, 1)
      allocate(l(n, n))
      l = 0.0_dp
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)
      do i = 1, n
         diagonal = a(i, i)
         if (i > 1) diagonal = diagonal - dot_product(l(i, 1:i - 1), l(i, 1:i - 1))
         if (diagonal <= threshold) then
            info = 2
            allocate(x(0, 0))
            return
         end if
         l(i, i) = sqrt(diagonal)
         do j = i + 1, n
            diagonal = a(j, i)
            if (i > 1) diagonal = diagonal - dot_product(l(j, 1:i - 1), l(i, 1:i - 1))
            l(j, i) = diagonal/l(i, i)
         end do
      end do

      allocate(y(n, size(b, 2)))
      y = b
      do i = 1, n
         if (i > 1) y(i, :) = y(i, :) - matmul(l(i, 1:i - 1), y(1:i - 1, :))
         y(i, :) = y(i, :)/l(i, i)
      end do
      allocate(x(n, size(b, 2)))
      x = y
      do i = n, 1, -1
         if (i < n) x(i, :) = x(i, :) - matmul(l(i + 1:n, i), x(i + 1:n, :))
         x(i, :) = x(i, :)/l(i, i)
      end do
      info = 0
   end subroutine llt_solve

   pure subroutine ldlt_solve(a, b, x, info)
      !! Symmetric-pivoted LDLT factorization and solve for a nonsingular self-adjoint matrix.
      real(dp), intent(in) :: a(:, :) !! Symmetric coefficient matrix; positive or indefinite matrices are supported.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with size(a,1) rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix in the original variable ordering.
      integer, intent(out) :: info !! Status: zero success, two if no usable one-by-one D pivot remains.
      integer, allocatable :: perm(:)
      real(dp), allocatable :: d(:), l(:, :), rhs_perm(:, :), temporary(:), work(:, :), y(:, :), z(:, :)
      real(dp) :: pivot_value, threshold
      integer :: i, j, k, n, pivot_index, temporary_index

      n = size(a, 1)
      allocate(work(n, n), l(n, n), d(n), perm(n))
      work = a
      l = 0.0_dp
      perm = [(i, i=1, n)]
      do i = 1, n
         l(i, i) = 1.0_dp
      end do
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(1, n), dp)

      do k = 1, n
         pivot_index = k
         pivot_value = abs(work(k, k))
         do i = k + 1, n
            if (abs(work(i, i)) > pivot_value) then
               pivot_value = abs(work(i, i))
               pivot_index = i
            end if
         end do
         if (pivot_value <= threshold) then
            info = 2
            allocate(x(0, 0))
            return
         end if
         if (pivot_index /= k) then
            temporary = work(k, :)
            work(k, :) = work(pivot_index, :)
            work(pivot_index, :) = temporary
            temporary = work(:, k)
            work(:, k) = work(:, pivot_index)
            work(:, pivot_index) = temporary
            if (k > 1) then
               temporary = l(k, 1:k - 1)
               l(k, 1:k - 1) = l(pivot_index, 1:k - 1)
               l(pivot_index, 1:k - 1) = temporary
            end if
            temporary_index = perm(k)
            perm(k) = perm(pivot_index)
            perm(pivot_index) = temporary_index
         end if

         d(k) = work(k, k)
         do i = k + 1, n
            l(i, k) = work(i, k)/d(k)
         end do
         do j = k + 1, n
            do i = j, n
               work(i, j) = work(i, j) - l(i, k)*d(k)*l(j, k)
               work(j, i) = work(i, j)
            end do
         end do
      end do

      allocate(rhs_perm(n, size(b, 2)))
      do i = 1, n
         rhs_perm(i, :) = b(perm(i), :)
      end do
      y = rhs_perm
      do i = 1, n
         if (i > 1) y(i, :) = y(i, :) - matmul(l(i, 1:i - 1), y(1:i - 1, :))
      end do
      allocate(z(n, size(b, 2)))
      do i = 1, n
         z(i, :) = y(i, :)/d(i)
      end do
      y = z
      do i = n, 1, -1
         if (i < n) y(i, :) = y(i, :) - matmul(l(i + 1:n, i), y(i + 1:n, :))
      end do
      allocate(x(n, size(b, 2)))
      x = 0.0_dp
      do i = 1, n
         x(perm(i), :) = y(i, :)
      end do
      info = 0
   end subroutine ldlt_solve

   pure subroutine full_pivot_qr_least_squares(a, b, x, info)
      !! Householder solve with complete row/column pivoting on each trailing block.
      real(dp), intent(in) :: a(:, :) !! Dense m-by-n coefficient matrix; tall, square, or wide.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with m rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Rank-revealing basic least-squares solution with n rows.
      integer, intent(out) :: info !! Status: zero success, one invalid shape.
      real(dp), allocatable :: r(:, :), rhs_work(:, :), temporary(:), vector(:), z(:, :)
      integer, allocatable :: perm(:)
      real(dp) :: alpha, beta, denominator, largest, norm_value, projection, threshold
      integer :: i, j, k, kmax, m, n, pivot_col, pivot_row, rank_value, temporary_index

      m = size(a, 1)
      n = size(a, 2)
      if (m < 1 .or. size(b, 1) /= m .or. n < 1) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      r = a
      rhs_work = b
      allocate(perm(n))
      perm = [(i, i=1, n)]
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(m, n), dp)
      kmax = min(m, n)
      rank_value = kmax

      do k = 1, kmax
         largest = -1.0_dp
         pivot_row = k
         pivot_col = k
         do j = k, n
            do i = k, m
               if (abs(r(i, j)) > largest) then
                  largest = abs(r(i, j))
                  pivot_row = i
                  pivot_col = j
               end if
            end do
         end do
         if (largest <= threshold) then
            rank_value = k - 1
            exit
         end if
         if (pivot_row /= k) then
            temporary = r(k, :)
            r(k, :) = r(pivot_row, :)
            r(pivot_row, :) = temporary
            temporary = rhs_work(k, :)
            rhs_work(k, :) = rhs_work(pivot_row, :)
            rhs_work(pivot_row, :) = temporary
         end if
         if (pivot_col /= k) then
            temporary = r(:, k)
            r(:, k) = r(:, pivot_col)
            r(:, pivot_col) = temporary
            temporary_index = perm(k)
            perm(k) = perm(pivot_col)
            perm(pivot_col) = temporary_index
         end if
         vector = r(k:m, k)
         norm_value = sqrt(max(0.0_dp, dot_product(vector, vector)))
         if (norm_value <= threshold) then
            rank_value = k - 1
            exit
         end if
         alpha = -sign(norm_value, vector(1))
         vector(1) = vector(1) - alpha
         denominator = dot_product(vector, vector)
         if (denominator <= tiny(1.0_dp)) then
            rank_value = k - 1
            exit
         end if
         beta = 2.0_dp/denominator
         do j = k, n
            projection = beta*dot_product(vector, r(k:m, j))
            r(k:m, j) = r(k:m, j) - vector*projection
         end do
         do j = 1, size(rhs_work, 2)
            projection = beta*dot_product(vector, rhs_work(k:m, j))
            rhs_work(k:m, j) = rhs_work(k:m, j) - vector*projection
         end do
         r(k, k) = alpha
         if (k < m) r(k + 1:m, k) = 0.0_dp
      end do

      allocate(z(n, size(b, 2)))
      z = 0.0_dp
      if (rank_value > 0) then
         z(1:rank_value, :) = rhs_work(1:rank_value, :)
         do i = rank_value, 1, -1
            if (i < rank_value) then
               z(i, :) = z(i, :) - matmul(r(i, i + 1:rank_value), z(i + 1:rank_value, :))
            end if
            z(i, :) = z(i, :)/r(i, i)
         end do
      end if
      allocate(x(n, size(b, 2)))
      x = 0.0_dp
      do i = 1, n
         x(perm(i), :) = z(i, :)
      end do
      info = 0
   end subroutine full_pivot_qr_least_squares

   pure subroutine qr_least_squares(a, b, x, info, pivoting)
      !! Householder solve with optional rank-revealing column pivoting and Eigen-style basic solutions.
      real(dp), intent(in) :: a(:, :) !! Dense m-by-n coefficient matrix; tall, square, or wide.
      real(dp), intent(in) :: b(:, :) !! RHS matrix with m rows.
      real(dp), allocatable, intent(out) :: x(:, :) !! Least-squares solution with n rows.
      integer, intent(out) :: info !! Status: zero success, two for rank deficiency when pivoting is disabled.
      logical, intent(in) :: pivoting !! If true, choose columns by residual norm and allow numerical rank deficiency.
      real(dp), allocatable :: beta(:), r(:, :), reflectors(:, :), rhs_work(:, :), z(:, :)
      integer, allocatable :: perm(:)
      integer :: i, k, n, rank_value

      n = size(a, 2)
      call householder_qr_factor(a, r, reflectors, beta, perm, rank_value, info, pivoting)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if
      rhs_work = b
      do k = 1, rank_value
         call apply_householder_left(reflectors(:, k), beta(k), k, rhs_work)
      end do

      allocate(z(n, size(b, 2)))
      z = 0.0_dp
      if (rank_value > 0) then
         z(1:rank_value, :) = rhs_work(1:rank_value, :)
         do i = rank_value, 1, -1
            if (i < rank_value) then
               z(i, :) = z(i, :) - matmul(r(i, i + 1:rank_value), z(i + 1:rank_value, :))
            end if
            z(i, :) = z(i, :)/r(i, i)
         end do
      end if
      allocate(x(n, size(b, 2)))
      x = 0.0_dp
      do i = 1, n
         x(perm(i), :) = z(i, :)
      end do
   end subroutine qr_least_squares



   pure subroutine householder_qr_factor(a, r, reflectors, beta, perm, rank_value, info, pivoting)
      !! Factor a real matrix with Householder reflectors and optional residual-norm column pivoting.
      real(dp), intent(in) :: a(:, :) !! Dense m-by-n matrix; at most min(m,n) reflectors are generated.
      real(dp), allocatable, intent(out) :: r(:, :) !! Upper-trapezoidal factor in the permuted column order.
      real(dp), allocatable, intent(out) :: reflectors(:, :) !! Householder vectors stored in their full row coordinates.
      real(dp), allocatable, intent(out) :: beta(:) !! Householder scale factors, one for each matrix column.
      integer, allocatable, intent(out) :: perm(:) !! New-index to original-index column permutation.
      integer, intent(out) :: rank_value !! Numerical rank detected from diagonal pivots.
      integer, intent(out) :: info !! Status: zero success, one invalid shape, two unpivoted rank deficiency.
      logical, intent(in) :: pivoting !! Whether to use residual column norms for rank-revealing pivoting.
      real(dp), allocatable :: temporary_column(:), vector(:)
      real(dp) :: alpha, denominator, norm_value, projection, threshold
      integer :: i, j, k, kmax, m, n, p, temporary_index

      m = size(a, 1)
      n = size(a, 2)
      if (m < 1 .or. n < 1) then
         allocate(r(0, 0), reflectors(0, 0), beta(0), perm(0))
         rank_value = 0
         info = 1
         return
      end if
      r = a
      kmax = min(m, n)
      allocate(reflectors(m, kmax), beta(kmax), perm(n))
      reflectors = 0.0_dp
      beta = 0.0_dp
      perm = [(i, i=1, n)]
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(max(m, n), dp)
      rank_value = kmax

      do k = 1, kmax
         if (pivoting) then
            p = k
            norm_value = -1.0_dp
            do j = k, n
               denominator = dot_product(r(k:m, j), r(k:m, j))
               if (denominator > norm_value) then
                  norm_value = denominator
                  p = j
               end if
            end do
            if (p /= k) then
               temporary_column = r(:, k)
               r(:, k) = r(:, p)
               r(:, p) = temporary_column
               temporary_index = perm(k)
               perm(k) = perm(p)
               perm(p) = temporary_index
            end if
         end if

         vector = r(k:m, k)
         norm_value = sqrt(max(0.0_dp, dot_product(vector, vector)))
         if (norm_value <= threshold) then
            if (.not. pivoting) then
               rank_value = k - 1
               info = 2
               return
            end if
            rank_value = k - 1
            exit
         end if
         alpha = -sign(norm_value, vector(1))
         vector(1) = vector(1) - alpha
         denominator = dot_product(vector, vector)
         if (denominator <= tiny(1.0_dp)) then
            rank_value = k - 1
            if (.not. pivoting) then
               info = 2
               return
            end if
            exit
         end if
         beta(k) = 2.0_dp/denominator
         reflectors(k:m, k) = vector
         do j = k, n
            projection = beta(k)*dot_product(vector, r(k:m, j))
            r(k:m, j) = r(k:m, j) - vector*projection
         end do
         r(k, k) = alpha
         if (k < m) r(k + 1:m, k) = 0.0_dp
      end do
      info = 0
   end subroutine householder_qr_factor

   pure subroutine apply_householder_left(reflector, beta, first_row, matrix)
      !! Premultiply a dense matrix by one stored Householder reflector.
      real(dp), intent(in) :: reflector(:) !! Full-length reflector vector with zeros before first_row.
      real(dp), intent(in) :: beta !! Householder scale factor; zero leaves the matrix unchanged.
      integer, intent(in) :: first_row !! First row on which the reflector has nonzero support.
      real(dp), intent(inout) :: matrix(:, :) !! Dense matrix updated in place as H*matrix.
      real(dp), allocatable :: projection(:)

      if (abs(beta) <= tiny(1.0_dp)) return
      projection = beta*matmul(reflector(first_row:), matrix(first_row:, :))
      matrix(first_row:, :) = matrix(first_row:, :) - spread(reflector(first_row:), 2, size(matrix, 2))* &
                               spread(projection, 1, size(matrix, 1) - first_row + 1)
   end subroutine apply_householder_left

end module sanic_direct
