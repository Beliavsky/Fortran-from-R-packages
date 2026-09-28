! SPDX-License-Identifier: GPL-3.0-or-later
! Iterative Krylov solvers translated from sanic's Eigen-backed wrappers.
module sanic_iterative
   use sanic_kinds, only: dp
   use sanic_sparse, only: sparse_matrix, sparse_matvec, sparse_transpose_matvec
   use sanic_sparse_direct, only: sparse_column_permutation
   implicit none
   private

   integer, parameter, public :: sanic_bicgstab = 1
   integer, parameter, public :: sanic_lscg = 2
   integer, parameter, public :: sanic_cg = 3

   integer, parameter :: pc_identity = 0
   integer, parameter :: pc_diagonal = 1
   integer, parameter :: pc_ilut = 2
   integer, parameter :: pc_ichol = 3

   type :: preconditioner_factor
      integer :: kind = pc_identity
      integer :: n = 0
      real(dp), allocatable :: diag_inv(:)
      integer, allocatable :: row_ptr(:)
      integer, allocatable :: col(:)
      real(dp), allocatable :: value(:)
      integer, allocatable :: t_row_ptr(:)
      integer, allocatable :: t_col(:)
      real(dp), allocatable :: t_value(:)
      integer, allocatable :: perm(:)
      integer, allocatable :: inverse_perm(:)
      real(dp), allocatable :: scale(:)
   end type preconditioner_factor

   type :: factor_row
      integer, allocatable :: col(:)
      real(dp), allocatable :: value(:)
   end type factor_row

   type :: factor_column
      integer, allocatable :: row(:)
      real(dp), allocatable :: value(:)
   end type factor_column

   public :: solve_cg_sparse

contains

   pure subroutine solve_cg_sparse(a, x, info, b, x0, method, iter, tol, precond, &
                                   iterations, error, converged)
      !! Solve sparse systems with BiCGSTAB, least-squares CG, or classic CG.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix; CG/BiCGSTAB require square shape.
      real(dp), allocatable, intent(out) :: x(:, :) !! Solution matrix with a%ncol rows and one column per RHS.
      integer, intent(out) :: info !! Status: zero success, one invalid input/setup, two nonconvergence/breakdown.
      real(dp), intent(in), optional :: b(:, :) !! Dense RHS with a%nrow rows; absent means an identity RHS.
      real(dp), intent(in), optional :: x0(:, :) !! Initial guess with a%ncol rows and the same number of RHS columns.
      integer, intent(in), optional :: method !! Solver selector: sanic_bicgstab, sanic_lscg, or sanic_cg; default BiCGSTAB.
      integer, intent(in), optional :: iter !! Maximum iterations per RHS; zero/absent defaults to 2*a%ncol.
      real(dp), intent(in), optional :: tol !! Relative residual tolerance; zero/absent uses machine epsilon.
      integer, intent(in), optional :: precond !! 0 identity, 1 diagonal/Jacobi, 2 incomplete factorization when supported.
      integer, intent(out), optional :: iterations !! Maximum iteration count used across all RHS columns.
      real(dp), intent(out), optional :: error !! Maximum final relative residual estimate across RHS columns.
      logical, intent(out), optional :: converged !! True only if every RHS converged.
      real(dp), allocatable :: guess(:, :), rhs(:, :), solution(:)
      type(preconditioner_factor) :: factor
      real(dp) :: err_col, tolerance
      integer :: i, it_col, max_iter, max_iterations_used, nrhs, preconditioner, solver, status
      real(dp) :: max_error
      logical :: all_converged, conv_col

      solver = sanic_bicgstab
      if (present(method)) solver = method
      preconditioner = 1
      if (present(precond)) preconditioner = precond
      info = 0
      if (present(iter)) then
         if (iter < 0) then
            allocate(x(0, 0))
            info = 1
            call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
            return
         end if
      end if
      if (present(tol)) then
         if (tol < 0.0_dp) then
            allocate(x(0, 0))
            info = 1
            call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
            return
         end if
      end if
      max_iter = max(1, 2*a%ncol)
      if (present(iter)) then
         if (iter > 0) max_iter = iter
      end if
      tolerance = epsilon(1.0_dp)
      if (present(tol)) then
         if (tol > 0.0_dp) tolerance = tol
      end if

      if (a%nrow < 1 .or. a%ncol < 1 .or. solver < sanic_bicgstab .or. solver > sanic_cg .or. &
          preconditioner < 0 .or. preconditioner > 2) then
         allocate(x(0, 0))
         info = 1
         call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
         return
      end if
      if (solver /= sanic_lscg .and. a%nrow /= a%ncol) then
         allocate(x(0, 0))
         info = 1
         call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
         return
      end if

      if (present(b)) then
         if (size(b, 1) /= a%nrow) then
            allocate(x(0, 0))
            info = 1
            call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
            return
         end if
         rhs = b
      else
         allocate(rhs(a%nrow, a%nrow))
         rhs = 0.0_dp
         do i = 1, a%nrow
            rhs(i, i) = 1.0_dp
         end do
      end if
      nrhs = size(rhs, 2)

      allocate(guess(a%ncol, nrhs))
      guess = 0.0_dp
      if (present(x0)) then
         if (size(x0, 1) /= a%ncol .or. size(x0, 2) /= nrhs) then
            allocate(x(0, 0))
            info = 1
            call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
            return
         end if
         guess = x0
      end if

      call build_preconditioner(a, solver, preconditioner, factor, status)
      if (status /= 0) then
         allocate(x(0, 0))
         info = 1
         call assign_diagnostics(0, huge(1.0_dp), .false., iterations, error, converged)
         return
      end if

      allocate(x(a%ncol, nrhs))
      x = 0.0_dp
      all_converged = .true.
      err_col = 0.0_dp
      it_col = 0
      max_iterations_used = 0
      max_error = 0.0_dp
      do i = 1, nrhs
         select case (solver)
         case (sanic_bicgstab)
            call bicgstab_vector(a, rhs(:, i), guess(:, i), solution, max_iter, tolerance, &
                                 factor, status, conv_col, it_col, err_col)
         case (sanic_lscg)
            call lscg_vector(a, rhs(:, i), guess(:, i), solution, max_iter, tolerance, &
                             factor, status, conv_col, it_col, err_col)
         case (sanic_cg)
            call cg_vector(a, rhs(:, i), guess(:, i), solution, max_iter, tolerance, &
                           factor, status, conv_col, it_col, err_col)
         end select
         if (status == 1) then
            info = 1
            x(:, i) = guess(:, i)
         else
            x(:, i) = solution
            if (.not. conv_col .and. info == 0) info = 2
         end if
         all_converged = all_converged .and. conv_col
         max_iterations_used = max(max_iterations_used, it_col)
         max_error = max(max_error, err_col)
      end do
      if (present(iterations)) iterations = max_iterations_used
      if (present(error)) error = max_error
      if (present(converged)) converged = all_converged
   end subroutine solve_cg_sparse

   pure subroutine cg_vector(a, b, x0, x, max_iter, tol, factor, info, converged, iterations, error)
      !! Classic preconditioned conjugate-gradient solve for one SPD RHS.
      type(sparse_matrix), intent(in) :: a !! Sparse symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:) !! RHS vector of length a%nrow.
      real(dp), intent(in) :: x0(:) !! Initial solution guess of length a%ncol.
      real(dp), allocatable, intent(out) :: x(:) !! Final solution estimate.
      integer, intent(in) :: max_iter !! Positive maximum iteration count.
      real(dp), intent(in) :: tol !! Positive relative-residual convergence tolerance.
      type(preconditioner_factor), intent(in) :: factor !! Prepared identity, Jacobi, or incomplete-Cholesky preconditioner.
      integer, intent(out) :: info !! Status: zero normal completion, one invalid dimensions/apply failure, two breakdown.
      logical, intent(out) :: converged !! True if the residual criterion was met.
      integer, intent(out) :: iterations !! Number of iterations executed.
      real(dp), intent(out) :: error !! Final relative residual norm.
      real(dp), allocatable :: ax(:), p(:), r(:), z(:), ap(:)
      real(dp) :: alpha, beta, bnorm, denom, rho, rho_new
      integer :: k, status

      if (size(b) /= a%nrow .or. size(x0) /= a%ncol .or. a%nrow /= a%ncol) then
         allocate(x(0))
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      x = x0
      call sparse_matvec(a, x, ax, status)
      if (status /= 0) then
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      r = b - ax
      call apply_preconditioner(factor, r, z, status)
      if (status /= 0) then
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      p = z
      rho = dot_product(r, z)
      bnorm = max(sqrt(dot_product(b, b)), tiny(1.0_dp))
      error = sqrt(dot_product(r, r))/bnorm
      converged = error <= tol
      iterations = 0
      info = 0
      if (converged) return

      do k = 1, max_iter
         call sparse_matvec(a, p, ap, status)
         denom = dot_product(p, ap)
         if (abs(denom) <= tiny(1.0_dp) .or. abs(rho) <= tiny(1.0_dp)) then
            info = 2
            exit
         end if
         alpha = rho/denom
         x = x + alpha*p
         r = r - alpha*ap
         error = sqrt(dot_product(r, r))/bnorm
         iterations = k
         if (error <= tol) then
            converged = .true.
            return
         end if
         call apply_preconditioner(factor, r, z, status)
         if (status /= 0) then
            info = 1
            converged = .false.
            return
         end if
         rho_new = dot_product(r, z)
         beta = rho_new/rho
         p = z + beta*p
         rho = rho_new
      end do
      converged = .false.
      if (info == 0) info = 2
   end subroutine cg_vector

   pure subroutine bicgstab_vector(a, b, x0, x, max_iter, tol, factor, info, converged, iterations, error)
      !! Preconditioned stabilized bi-conjugate-gradient solve for one RHS.
      type(sparse_matrix), intent(in) :: a !! Sparse square coefficient matrix.
      real(dp), intent(in) :: b(:) !! RHS vector of length a%nrow.
      real(dp), intent(in) :: x0(:) !! Initial solution guess of length a%ncol.
      real(dp), allocatable, intent(out) :: x(:) !! Final solution estimate.
      integer, intent(in) :: max_iter !! Positive maximum iteration count.
      real(dp), intent(in) :: tol !! Positive relative-residual convergence tolerance.
      type(preconditioner_factor), intent(in) :: factor !! Prepared identity, Jacobi, or dual-threshold ILUT preconditioner.
      integer, intent(out) :: info !! Status: zero normal completion, one invalid dimensions/apply failure, two breakdown.
      logical, intent(out) :: converged !! True if the residual criterion was met.
      integer, intent(out) :: iterations !! Number of iterations executed.
      real(dp), intent(out) :: error !! Final relative residual norm.
      real(dp), allocatable :: ax(:), p(:), phat(:), r(:), rhat(:), s(:), shat(:), t(:), v(:)
      real(dp) :: alpha, beta, bnorm, denom, omega, rho, rho_old
      integer :: k, status

      if (size(b) /= a%nrow .or. size(x0) /= a%ncol .or. a%nrow /= a%ncol) then
         allocate(x(0))
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      x = x0
      call sparse_matvec(a, x, ax, status)
      if (status /= 0) then
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      r = b - ax
      rhat = r
      allocate(p(a%ncol), v(a%nrow))
      p = 0.0_dp
      v = 0.0_dp
      rho_old = 1.0_dp
      alpha = 1.0_dp
      omega = 1.0_dp
      bnorm = max(sqrt(dot_product(b, b)), tiny(1.0_dp))
      error = sqrt(dot_product(r, r))/bnorm
      converged = error <= tol
      iterations = 0
      info = 0
      if (converged) return

      do k = 1, max_iter
         rho = dot_product(rhat, r)
         if (abs(rho) <= tiny(1.0_dp) .or. abs(omega) <= tiny(1.0_dp)) then
            info = 2
            exit
         end if
         if (k == 1) then
            p = r
         else
            beta = (rho/rho_old)*(alpha/omega)
            p = r + beta*(p - omega*v)
         end if
         call apply_preconditioner(factor, p, phat, status)
         if (status /= 0) then
            info = 1
            converged = .false.
            return
         end if
         call sparse_matvec(a, phat, v, status)
         denom = dot_product(rhat, v)
         if (abs(denom) <= tiny(1.0_dp)) then
            info = 2
            exit
         end if
         alpha = rho/denom
         s = r - alpha*v
         error = sqrt(dot_product(s, s))/bnorm
         if (error <= tol) then
            x = x + alpha*phat
            converged = .true.
            iterations = k
            return
         end if
         call apply_preconditioner(factor, s, shat, status)
         if (status /= 0) then
            info = 1
            converged = .false.
            return
         end if
         call sparse_matvec(a, shat, t, status)
         denom = dot_product(t, t)
         if (denom <= tiny(1.0_dp)) then
            info = 2
            exit
         end if
         omega = dot_product(t, s)/denom
         x = x + alpha*phat + omega*shat
         r = s - omega*t
         error = sqrt(dot_product(r, r))/bnorm
         iterations = k
         if (error <= tol) then
            converged = .true.
            return
         end if
         rho_old = rho
      end do
      converged = .false.
      if (info == 0) info = 2
   end subroutine bicgstab_vector

   pure subroutine lscg_vector(a, b, x0, x, max_iter, tol, factor, info, converged, iterations, error)
      !! Least-squares conjugate-gradient solve using A and transpose(A) products without forming normal equations.
      type(sparse_matrix), intent(in) :: a !! Sparse rectangular coefficient matrix.
      real(dp), intent(in) :: b(:) !! RHS vector of length a%nrow.
      real(dp), intent(in) :: x0(:) !! Initial solution guess of length a%ncol.
      real(dp), allocatable, intent(out) :: x(:) !! Final least-squares solution estimate.
      integer, intent(in) :: max_iter !! Positive maximum iteration count.
      real(dp), intent(in) :: tol !! Positive relative normal-residual convergence tolerance.
      type(preconditioner_factor), intent(in) :: factor !! Prepared identity or least-squares diagonal preconditioner.
      integer, intent(out) :: info !! Status: zero normal completion, one invalid dimensions/apply failure, two breakdown.
      logical, intent(out) :: converged !! True if the normal-residual criterion was met.
      integer, intent(out) :: iterations !! Number of iterations executed.
      real(dp), intent(out) :: error !! Final normal-residual norm relative to transpose(A)*b.
      real(dp), allocatable :: ax(:), atb(:), atr(:), p(:), q(:), r(:), z(:), z_new(:)
      real(dp) :: alpha, beta, denom, reference, rho, rho_new
      integer :: k, status

      if (size(b) /= a%nrow .or. size(x0) /= a%ncol) then
         allocate(x(0))
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      x = x0
      call sparse_matvec(a, x, ax, status)
      if (status /= 0) then
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      r = b - ax
      call sparse_transpose_matvec(a, r, atr, status)
      call sparse_transpose_matvec(a, b, atb, status)
      call apply_preconditioner(factor, atr, z, status)
      if (status /= 0) then
         info = 1
         converged = .false.
         iterations = 0
         error = huge(1.0_dp)
         return
      end if
      p = z
      rho = dot_product(atr, z)
      reference = max(sqrt(dot_product(atb, atb)), tiny(1.0_dp))
      error = sqrt(dot_product(atr, atr))/reference
      converged = error <= tol
      iterations = 0
      info = 0
      if (converged) return

      do k = 1, max_iter
         call sparse_matvec(a, p, q, status)
         denom = dot_product(q, q)
         if (denom <= tiny(1.0_dp) .or. abs(rho) <= tiny(1.0_dp)) then
            info = 2
            exit
         end if
         alpha = rho/denom
         x = x + alpha*p
         r = r - alpha*q
         call sparse_transpose_matvec(a, r, atr, status)
         error = sqrt(dot_product(atr, atr))/reference
         iterations = k
         if (error <= tol) then
            converged = .true.
            return
         end if
         call apply_preconditioner(factor, atr, z_new, status)
         if (status /= 0) then
            info = 1
            converged = .false.
            return
         end if
         rho_new = dot_product(atr, z_new)
         beta = rho_new/rho
         p = z_new + beta*p
         z = z_new
         rho = rho_new
      end do
      converged = .false.
      if (info == 0) info = 2
   end subroutine lscg_vector

   pure subroutine build_preconditioner(a, solver, precond, factor, info)
      !! Prepare the preconditioner corresponding to Eigen's solver-specific choices.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix to precondition.
      integer, intent(in) :: solver !! Solver selector controlling which incomplete factorization is available.
      integer, intent(in) :: precond !! Requested preconditioner: zero identity, one diagonal, two incomplete when supported.
      type(preconditioner_factor), intent(out) :: factor !! Prepared preconditioner data reused for every RHS.
      integer, intent(out) :: info !! Status: zero success, one when an incomplete-factorization setup fails.

      factor%n = a%ncol
      info = 0
      if (precond == 0) then
         factor%kind = pc_identity
         return
      end if

      if (precond == 1 .or. solver == sanic_lscg) then
         call build_diagonal_preconditioner(a, solver == sanic_lscg, factor)
         return
      end if

      select case (solver)
      case (sanic_bicgstab)
         call build_ilut_preconditioner(a, factor, info)
      case (sanic_cg)
         call build_ichol_preconditioner(a, factor, info)
      case default
         call build_diagonal_preconditioner(a, solver == sanic_lscg, factor)
      end select
   end subroutine build_preconditioner

   pure subroutine build_diagonal_preconditioner(a, normal_equations, factor)
      !! Form Eigen-like diagonal/Jacobi inverse scaling for A or the normal equations.
      type(sparse_matrix), intent(in) :: a !! Sparse coefficient matrix.
      logical, intent(in) :: normal_equations !! True for least-squares diagonal scaling based on diag(A^T A).
      type(preconditioner_factor), intent(out) :: factor !! Diagonal factor containing reciprocal safe diagonal entries.
      real(dp), allocatable :: diagonal(:)
      real(dp) :: threshold
      integer :: k, n

      n = a%ncol
      factor%n = n
      factor%kind = pc_diagonal
      allocate(factor%diag_inv(n), diagonal(n))
      factor%diag_inv = 1.0_dp
      diagonal = 0.0_dp
      if (allocated(a%value)) then
         if (normal_equations) then
            do k = 1, size(a%value)
               if (a%col(k) >= 1 .and. a%col(k) <= n) then
                  diagonal(a%col(k)) = diagonal(a%col(k)) + a%value(k)*a%value(k)
               end if
            end do
         else
            do k = 1, size(a%value)
               if (a%row(k) == a%col(k) .and. a%col(k) >= 1 .and. a%col(k) <= n) then
                  diagonal(a%col(k)) = diagonal(a%col(k)) + a%value(k)
               end if
            end do
         end if
      end if
      threshold = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(diagonal)))
      do k = 1, n
         if (abs(diagonal(k)) > threshold) factor%diag_inv(k) = 1.0_dp/diagonal(k)
      end do
   end subroutine build_diagonal_preconditioner

   pure subroutine build_ilut_preconditioner(a, factor, info)
      !! Build Eigen-style dual-threshold ILUT with symmetric fill-reducing permutation.
      type(sparse_matrix), intent(in) :: a !! Sparse square coefficient matrix to factor approximately.
      type(preconditioner_factor), intent(out) :: factor !! Combined unit-lower and upper ILUT factors in CSR storage.
      integer, intent(out) :: info !! Status: zero success, one for invalid shape, zero row, or unusable pivot.
      type(factor_row), allocatable :: rows(:)
      real(dp), allocatable :: work(:)
      integer, allocatable :: keep(:)
      real(dp), parameter :: droptol = 1.0e-12_dp
      integer, parameter :: fillfactor = 10
      real(dp) :: diagonal, fact, rownorm
      integer :: fill_in, i, j, k, keep_l, keep_u, n, nnz_total, old_row, p

      factor%n = a%ncol
      factor%kind = pc_ilut
      if (a%nrow /= a%ncol .or. a%ncol < 1) then
         info = 1
         return
      end if
      n = a%ncol
      call sparse_column_permutation(a, 0, factor%perm, factor%inverse_perm)
      nnz_total = 0
      if (allocated(a%value)) then
         do k = 1, size(a%value)
            if (a%row(k) < 1 .or. a%row(k) > n .or. a%col(k) < 1 .or. a%col(k) > n) cycle
            nnz_total = nnz_total + 1
         end do
      end if
      fill_in = min(n, (nnz_total*fillfactor)/n + 1)
      keep_l = fill_in/2
      keep_u = keep_l
      allocate(rows(n), work(n))

      do i = 1, n
         work = 0.0_dp
         old_row = factor%perm(i)
         if (allocated(a%value)) then
            do k = 1, size(a%value)
               if (a%row(k) /= old_row) cycle
               if (a%col(k) < 1 .or. a%col(k) > n) cycle
               j = factor%inverse_perm(a%col(k))
               work(j) = work(j) + a%value(k)
            end do
         end if
         rownorm = sqrt(dot_product(work, work))
         if (rownorm <= tiny(1.0_dp)) then
            info = 1
            return
         end if

         do k = 1, i - 1
            if (abs(work(k)) <= tiny(1.0_dp)) cycle
            diagonal = factor_row_get(rows(k), k)
            if (abs(diagonal) <= tiny(1.0_dp)) then
               info = 1
               return
            end if
            fact = work(k)/diagonal
            if (abs(fact) <= droptol) then
               work(k) = 0.0_dp
               cycle
            end if
            work(k) = fact
            if (allocated(rows(k)%col)) then
               do p = 1, size(rows(k)%col)
                  j = rows(k)%col(p)
                  if (j <= k) cycle
                  work(j) = work(j) - fact*rows(k)%value(p)
               end do
            end if
         end do

         call select_largest_entries(work, 1, i - 1, keep_l, 0.0_dp, keep)
         do p = 1, size(keep)
            call factor_row_set(rows(i), keep(p), work(keep(p)))
         end do
         diagonal = work(i)
         if (abs(diagonal) <= tiny(1.0_dp)) diagonal = sqrt(droptol)*rownorm
         if (abs(diagonal) <= tiny(1.0_dp)) then
            info = 1
            return
         end if
         call factor_row_set(rows(i), i, diagonal)
         call select_largest_entries(work, i + 1, n, keep_u, droptol*rownorm, keep)
         do p = 1, size(keep)
            call factor_row_set(rows(i), keep(p), work(keep(p)))
         end do
      end do

      call flatten_factor_rows(rows, factor%row_ptr, factor%col, factor%value)
      info = 0
   end subroutine build_ilut_preconditioner

   pure subroutine build_ichol_preconditioner(a, factor, info)
      !! Build Eigen-style scaled, permuted, shifted limited-memory incomplete Cholesky factors.
      type(sparse_matrix), intent(in) :: a !! Sparse symmetric coefficient matrix; the lower self-adjoint view is used.
      type(preconditioner_factor), intent(out) :: factor !! Sparse lower factor, transpose, scaling, and permutation data.
      integer, intent(out) :: info !! Status: zero success, one after ten unsuccessful shifted factorization attempts.
      type(sparse_matrix) :: lower
      type(factor_column), allocatable :: columns(:)
      integer, allocatable :: base_col(:), base_row_ptr(:), keep(:)
      real(dp), allocatable :: base_value(:), diagonal_base(:), diagonal_work(:), work(:)
      real(dp), parameter :: initial_shift = 1.0e-3_dp
      real(dp) :: diagonal, ljk, mindiag, shift
      integer :: attempt, i, j, k, n, original_keep, p, status
      logical :: success

      factor%n = a%ncol
      factor%kind = pc_ichol
      if (a%nrow /= a%ncol .or. a%ncol < 1) then
         info = 1
         return
      end if
      n = a%ncol
      call sparse_column_permutation(a, 0, factor%perm, factor%inverse_perm)
      call build_permuted_lower(a, factor%inverse_perm, lower)
      call csr_from_sparse(lower, base_row_ptr, base_col, base_value, status)
      if (status /= 0) then
         info = 1
         return
      end if
      allocate(factor%scale(n), diagonal_base(n), diagonal_work(n), work(n))
      factor%scale = 0.0_dp
      do i = 1, n
         do k = base_row_ptr(i), base_row_ptr(i + 1) - 1
            factor%scale(i) = factor%scale(i) + base_value(k)*base_value(k)
            if (base_col(k) /= i) then
               factor%scale(base_col(k)) = factor%scale(base_col(k)) + base_value(k)*base_value(k)
            end if
         end do
      end do
      do i = 1, n
         factor%scale(i) = sqrt(sqrt(factor%scale(i)))
         if (factor%scale(i) > tiny(1.0_dp)) then
            factor%scale(i) = 1.0_dp/factor%scale(i)
         else
            factor%scale(i) = 1.0_dp
         end if
      end do
      do i = 1, n
         do k = base_row_ptr(i), base_row_ptr(i + 1) - 1
            base_value(k) = base_value(k)*factor%scale(i)*factor%scale(base_col(k))
         end do
      end do

      mindiag = huge(1.0_dp)
      do i = 1, n
         k = csr_find_index(base_row_ptr, base_col, i, i)
         if (k == 0) then
            info = 1
            return
         end if
         diagonal_base(i) = base_value(k)
         mindiag = min(mindiag, diagonal_base(i))
      end do
      shift = 0.0_dp
      if (mindiag <= 0.0_dp) shift = initial_shift - mindiag

      do attempt = 1, 10
         if (allocated(columns)) deallocate(columns)
         allocate(columns(n))
         diagonal_work = diagonal_base + shift
         success = .true.
         do j = 1, n
            diagonal = diagonal_work(j)
            if (diagonal <= tiny(1.0_dp)) then
               success = .false.
               exit
            end if
            diagonal = sqrt(diagonal)
            call factor_column_set(columns(j), j, diagonal)
            work = 0.0_dp
            original_keep = 0
            do i = j + 1, n
               k = csr_find_index(base_row_ptr, base_col, i, j)
               if (k > 0) then
                  work(i) = base_value(k)
                  original_keep = original_keep + 1
               end if
            end do

            do k = 1, j - 1
               ljk = factor_column_get(columns(k), j)
               if (abs(ljk) <= tiny(1.0_dp)) cycle
               if (.not. allocated(columns(k)%row)) cycle
               do p = 1, size(columns(k)%row)
                  i = columns(k)%row(p)
                  if (i <= j) cycle
                  work(i) = work(i) - columns(k)%value(p)*ljk
               end do
            end do
            do i = j + 1, n
               if (abs(work(i)) <= tiny(1.0_dp)) cycle
               work(i) = work(i)/diagonal
               diagonal_work(i) = diagonal_work(i) - work(i)*work(i)
            end do
            call select_largest_entries(work, j + 1, n, original_keep, 0.0_dp, keep)
            do p = 1, size(keep)
               call factor_column_set(columns(j), keep(p), work(keep(p)))
            end do
         end do
         if (success) then
            call flatten_factor_columns(columns, factor%row_ptr, factor%col, factor%value)
            call csr_transpose(n, factor%row_ptr, factor%col, factor%value, &
                               factor%t_row_ptr, factor%t_col, factor%t_value)
            info = 0
            return
         end if
         shift = max(initial_shift, 2.0_dp*shift)
      end do
      info = 1
   end subroutine build_ichol_preconditioner

   pure subroutine build_permuted_lower(a, inverse_perm, lower)
      !! Form the permuted lower self-adjoint view used by incomplete Cholesky.
      type(sparse_matrix), intent(in) :: a !! Original sparse matrix; only its lower-triangular stored entries are used.
      integer, intent(in) :: inverse_perm(:) !! Original-index to permuted-index mapping.
      type(sparse_matrix), intent(out) :: lower !! Permuted lower-triangular coordinate matrix.
      integer, allocatable :: col_work(:), row_work(:)
      real(dp), allocatable :: value_work(:)
      integer :: count_value, i, j, k, n, new_i, new_j

      n = a%ncol
      lower%nrow = n
      lower%ncol = n
      if (.not. allocated(a%value)) then
         allocate(lower%row(0), lower%col(0), lower%value(0))
         return
      end if
      allocate(row_work(size(a%value)), col_work(size(a%value)), value_work(size(a%value)))
      count_value = 0
      do k = 1, size(a%value)
         i = a%row(k)
         j = a%col(k)
         if (i < 1 .or. i > n .or. j < 1 .or. j > n .or. i < j) cycle
         new_i = inverse_perm(i)
         new_j = inverse_perm(j)
         if (new_i < new_j) then
            i = new_i
            new_i = new_j
            new_j = i
         end if
         count_value = count_value + 1
         row_work(count_value) = new_i
         col_work(count_value) = new_j
         value_work(count_value) = a%value(k)
      end do
      allocate(lower%row(count_value), lower%col(count_value), lower%value(count_value))
      if (count_value > 0) then
         lower%row = row_work(1:count_value)
         lower%col = col_work(1:count_value)
         lower%value = value_work(1:count_value)
      end if
   end subroutine build_permuted_lower

   pure real(dp) function factor_column_get(column, target_row) result(value)
      !! Read one entry from a sorted factor column, returning zero when absent.
      type(factor_column), intent(in) :: column !! Sorted factor column to query.
      integer, intent(in) :: target_row !! One-based row index to query.
      integer :: i

      value = 0.0_dp
      if (.not. allocated(column%row)) return
      do i = 1, size(column%row)
         if (column%row(i) == target_row) then
            value = column%value(i)
            return
         end if
         if (column%row(i) > target_row) return
      end do
   end function factor_column_get

   pure subroutine factor_column_set(column, target_row, new_value)
      !! Insert one factor-column entry while preserving sorted row order.
      type(factor_column), intent(inout) :: column !! Factor column to modify.
      integer, intent(in) :: target_row !! One-based row index to insert or replace.
      real(dp), intent(in) :: new_value !! Numerical value stored at target_row.
      integer, allocatable :: new_row(:)
      real(dp), allocatable :: new_values(:)
      integer :: insert_at, n

      if (.not. allocated(column%row)) then
         allocate(column%row(1), column%value(1))
         column%row(1) = target_row
         column%value(1) = new_value
         return
      end if
      n = size(column%row)
      insert_at = 1
      do while (insert_at <= n)
         if (column%row(insert_at) >= target_row) exit
         insert_at = insert_at + 1
      end do
      if (insert_at <= n) then
         if (column%row(insert_at) == target_row) then
            column%value(insert_at) = new_value
            return
         end if
      end if
      allocate(new_row(n + 1), new_values(n + 1))
      if (insert_at > 1) then
         new_row(1:insert_at - 1) = column%row(1:insert_at - 1)
         new_values(1:insert_at - 1) = column%value(1:insert_at - 1)
      end if
      new_row(insert_at) = target_row
      new_values(insert_at) = new_value
      if (insert_at <= n) then
         new_row(insert_at + 1:n + 1) = column%row(insert_at:n)
         new_values(insert_at + 1:n + 1) = column%value(insert_at:n)
      end if
      call move_alloc(new_row, column%row)
      call move_alloc(new_values, column%value)
   end subroutine factor_column_set

   pure subroutine flatten_factor_columns(columns, row_ptr, col, value)
      !! Convert lower factor columns to row-sorted CSR arrays for triangular solves.
      type(factor_column), intent(in) :: columns(:) !! Sorted lower factor columns.
      integer, allocatable, intent(out) :: row_ptr(:) !! One-based CSR row pointers.
      integer, allocatable, intent(out) :: col(:) !! Sorted lower-factor column indices.
      real(dp), allocatable, intent(out) :: value(:) !! Numerical lower-factor entries.
      type(sparse_matrix) :: lower
      integer :: j, offset_value, p, total

      total = 0
      do j = 1, size(columns)
         if (allocated(columns(j)%row)) total = total + size(columns(j)%row)
      end do
      lower%nrow = size(columns)
      lower%ncol = size(columns)
      allocate(lower%row(total), lower%col(total), lower%value(total))
      offset_value = 0
      do j = 1, size(columns)
         if (.not. allocated(columns(j)%row)) cycle
         do p = 1, size(columns(j)%row)
            offset_value = offset_value + 1
            lower%row(offset_value) = columns(j)%row(p)
            lower%col(offset_value) = j
            lower%value(offset_value) = columns(j)%value(p)
         end do
      end do
      call csr_from_sparse(lower, row_ptr, col, value, p)
   end subroutine flatten_factor_columns

   pure subroutine select_largest_entries(values, first_index, last_index, max_keep, drop_threshold, indices)
      !! Select up to max_keep entries of greatest magnitude above a drop threshold, returning sorted indices.
      real(dp), intent(in) :: values(:) !! Dense row workspace containing candidate values.
      integer, intent(in) :: first_index !! First candidate index, inclusive.
      integer, intent(in) :: last_index !! Last candidate index, inclusive.
      integer, intent(in) :: max_keep !! Maximum number of entries retained.
      real(dp), intent(in) :: drop_threshold !! Entries at or below this absolute magnitude are dropped.
      integer, allocatable, intent(out) :: indices(:) !! Retained indices sorted in ascending column order.
      integer, allocatable :: candidates(:)
      integer :: count_value, current, i, j, keep_count

      if (first_index > last_index .or. max_keep <= 0) then
         allocate(indices(0))
         return
      end if
      allocate(candidates(last_index - first_index + 1))
      count_value = 0
      do i = first_index, last_index
         if (abs(values(i)) <= drop_threshold) cycle
         count_value = count_value + 1
         candidates(count_value) = i
      end do
      do i = 2, count_value
         current = candidates(i)
         j = i - 1
         do while (j >= 1)
            if (abs(values(candidates(j))) >= abs(values(current))) exit
            candidates(j + 1) = candidates(j)
            j = j - 1
         end do
         candidates(j + 1) = current
      end do
      keep_count = min(max_keep, count_value)
      allocate(indices(keep_count))
      if (keep_count > 0) indices = candidates(1:keep_count)
      call sort_integer_ascending(indices)
   end subroutine select_largest_entries

   pure subroutine sort_integer_ascending(values)
      !! Sort a short integer index list in ascending order using insertion sort.
      integer, intent(inout) :: values(:) !! Integer values to sort in place.
      integer :: current, i, j

      do i = 2, size(values)
         current = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= current) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = current
      end do
   end subroutine sort_integer_ascending

   pure real(dp) function factor_row_get(row, target_col) result(value)
      !! Read one entry from a sorted factor row, returning zero when absent.
      type(factor_row), intent(in) :: row !! Sorted factor row to query.
      integer, intent(in) :: target_col !! One-based column index to query.
      integer :: i

      value = 0.0_dp
      if (.not. allocated(row%col)) return
      do i = 1, size(row%col)
         if (row%col(i) == target_col) then
            value = row%value(i)
            return
         end if
         if (row%col(i) > target_col) return
      end do
   end function factor_row_get

   pure subroutine factor_row_set(row, target_col, new_value)
      !! Insert one factor entry while preserving sorted column order.
      type(factor_row), intent(inout) :: row !! Factor row to modify.
      integer, intent(in) :: target_col !! One-based column index to insert or replace.
      real(dp), intent(in) :: new_value !! Numerical value stored at target_col.
      integer, allocatable :: new_col(:)
      real(dp), allocatable :: new_values(:)
      integer :: insert_at, n

      if (.not. allocated(row%col)) then
         allocate(row%col(1), row%value(1))
         row%col(1) = target_col
         row%value(1) = new_value
         return
      end if
      n = size(row%col)
      insert_at = 1
      do while (insert_at <= n)
         if (row%col(insert_at) >= target_col) exit
         insert_at = insert_at + 1
      end do
      if (insert_at <= n) then
         if (row%col(insert_at) == target_col) then
            row%value(insert_at) = new_value
            return
         end if
      end if
      allocate(new_col(n + 1), new_values(n + 1))
      if (insert_at > 1) then
         new_col(1:insert_at - 1) = row%col(1:insert_at - 1)
         new_values(1:insert_at - 1) = row%value(1:insert_at - 1)
      end if
      new_col(insert_at) = target_col
      new_values(insert_at) = new_value
      if (insert_at <= n) then
         new_col(insert_at + 1:n + 1) = row%col(insert_at:n)
         new_values(insert_at + 1:n + 1) = row%value(insert_at:n)
      end if
      call move_alloc(new_col, row%col)
      call move_alloc(new_values, row%value)
   end subroutine factor_row_set

   pure subroutine flatten_factor_rows(rows, row_ptr, col, value)
      !! Flatten dynamic factor rows into one-based CSR arrays.
      type(factor_row), intent(in) :: rows(:) !! Sorted factor rows to flatten.
      integer, allocatable, intent(out) :: row_ptr(:) !! One-based CSR row pointers.
      integer, allocatable, intent(out) :: col(:) !! Flattened sorted column indices.
      real(dp), allocatable, intent(out) :: value(:) !! Flattened numerical factor entries.
      integer :: i, offset_value, total

      allocate(row_ptr(size(rows) + 1))
      row_ptr(1) = 1
      do i = 1, size(rows)
         if (allocated(rows(i)%col)) then
            row_ptr(i + 1) = row_ptr(i) + size(rows(i)%col)
         else
            row_ptr(i + 1) = row_ptr(i)
         end if
      end do
      total = row_ptr(size(rows) + 1) - 1
      allocate(col(total), value(total))
      offset_value = 0
      do i = 1, size(rows)
         if (.not. allocated(rows(i)%col)) cycle
         col(offset_value + 1:offset_value + size(rows(i)%col)) = rows(i)%col
         value(offset_value + 1:offset_value + size(rows(i)%value)) = rows(i)%value
         offset_value = offset_value + size(rows(i)%col)
      end do
   end subroutine flatten_factor_rows

   pure subroutine apply_preconditioner(factor, rhs, solution, info)
      !! Apply a prepared identity, Jacobi, ILUT, or incomplete-Cholesky preconditioner.
      type(preconditioner_factor), intent(in) :: factor !! Prepared factorization or diagonal scaling.
      real(dp), intent(in) :: rhs(:) !! Vector to which M^{-1} is applied.
      real(dp), allocatable, intent(out) :: solution(:) !! Preconditioned vector.
      integer, intent(out) :: info !! Status: zero success, one for incompatible dimensions or invalid pivots.
      real(dp), allocatable :: work(:)
      real(dp) :: diagonal
      integer :: diag_pos, i, j, p

      if (size(rhs) /= factor%n) then
         allocate(solution(0))
         info = 1
         return
      end if
      select case (factor%kind)
      case (pc_identity)
         solution = rhs
         info = 0
      case (pc_diagonal)
         if (.not. allocated(factor%diag_inv) .or. size(factor%diag_inv) /= factor%n) then
            allocate(solution(0))
            info = 1
            return
         end if
         solution = factor%diag_inv*rhs
         info = 0
      case (pc_ilut)
         allocate(work(factor%n), solution(factor%n))
         do i = 1, factor%n
            work(i) = rhs(factor%perm(i))
         end do
         do i = 1, factor%n
            do p = factor%row_ptr(i), factor%row_ptr(i + 1) - 1
               j = factor%col(p)
               if (j >= i) exit
               work(i) = work(i) - factor%value(p)*work(j)
            end do
         end do
         solution = work
         do i = factor%n, 1, -1
            diag_pos = csr_find_index(factor%row_ptr, factor%col, i, i)
            if (diag_pos == 0) then
               info = 1
               return
            end if
            diagonal = factor%value(diag_pos)
            if (abs(diagonal) <= tiny(1.0_dp)) then
               info = 1
               return
            end if
            do p = diag_pos + 1, factor%row_ptr(i + 1) - 1
               j = factor%col(p)
               solution(i) = solution(i) - factor%value(p)*solution(j)
            end do
            solution(i) = solution(i)/diagonal
         end do
         work = solution
         do i = 1, factor%n
            solution(factor%perm(i)) = work(i)
         end do
         info = 0
      case (pc_ichol)
         allocate(work(factor%n), solution(factor%n))
         do i = 1, factor%n
            work(i) = factor%scale(i)*rhs(factor%perm(i))
         end do
         do i = 1, factor%n
            diag_pos = csr_find_index(factor%row_ptr, factor%col, i, i)
            if (diag_pos == 0) then
               info = 1
               return
            end if
            diagonal = factor%value(diag_pos)
            if (diagonal <= tiny(1.0_dp)) then
               info = 1
               return
            end if
            do p = factor%row_ptr(i), diag_pos - 1
               j = factor%col(p)
               work(i) = work(i) - factor%value(p)*work(j)
            end do
            work(i) = work(i)/diagonal
         end do
         solution = work
         do i = factor%n, 1, -1
            diag_pos = csr_find_index(factor%row_ptr, factor%col, i, i)
            diagonal = factor%value(diag_pos)
            do p = factor%t_row_ptr(i), factor%t_row_ptr(i + 1) - 1
               j = factor%t_col(p)
               if (j > i) solution(i) = solution(i) - factor%t_value(p)*solution(j)
            end do
            solution(i) = solution(i)/diagonal
         end do
         work = solution
         do i = 1, factor%n
            solution(factor%perm(i)) = factor%scale(i)*work(i)
         end do
         info = 0
      case default
         allocate(solution(0))
         info = 1
      end select
   end subroutine apply_preconditioner

   pure subroutine csr_from_sparse(a, row_ptr, col, value, info)
      !! Convert coordinate sparse storage to row-sorted CSR while summing duplicate coordinates.
      type(sparse_matrix), intent(in) :: a !! Coordinate sparse matrix to canonicalize.
      integer, allocatable, intent(out) :: row_ptr(:) !! One-based CSR row pointers of length a%nrow+1.
      integer, allocatable, intent(out) :: col(:) !! One-based sorted column indices for canonical nonzeros.
      real(dp), allocatable, intent(out) :: value(:) !! Values matching col, with duplicate coordinates summed.
      integer, intent(out) :: info !! Status: zero success, one for inconsistent coordinate-array lengths.
      integer, allocatable :: col_tmp(:), counts(:), next_slot(:), unique_counts(:)
      real(dp), allocatable :: value_tmp(:)
      integer :: current_col, i, j, k, nvalid, p, q, start_pos, stop_pos, total_unique
      real(dp) :: current_value, value_swap
      integer :: col_swap

      allocate(row_ptr(a%nrow + 1))
      row_ptr = 1
      if (.not. allocated(a%row) .or. .not. allocated(a%col) .or. .not. allocated(a%value)) then
         allocate(col(0), value(0))
         info = 0
         return
      end if
      if (size(a%row) /= size(a%col) .or. size(a%row) /= size(a%value)) then
         allocate(col(0), value(0))
         info = 1
         return
      end if

      allocate(counts(a%nrow))
      counts = 0
      nvalid = 0
      do k = 1, size(a%value)
         if (a%row(k) >= 1 .and. a%row(k) <= a%nrow .and. &
             a%col(k) >= 1 .and. a%col(k) <= a%ncol) then
            counts(a%row(k)) = counts(a%row(k)) + 1
            nvalid = nvalid + 1
         end if
      end do
      row_ptr(1) = 1
      do i = 1, a%nrow
         row_ptr(i + 1) = row_ptr(i) + counts(i)
      end do
      allocate(col_tmp(nvalid), value_tmp(nvalid), next_slot(a%nrow))
      next_slot = row_ptr(1:a%nrow)
      do k = 1, size(a%value)
         i = a%row(k)
         j = a%col(k)
         if (i >= 1 .and. i <= a%nrow .and. j >= 1 .and. j <= a%ncol) then
            p = next_slot(i)
            col_tmp(p) = j
            value_tmp(p) = a%value(k)
            next_slot(i) = p + 1
         end if
      end do

      do i = 1, a%nrow
         start_pos = row_ptr(i)
         stop_pos = row_ptr(i + 1) - 1
         do p = start_pos + 1, stop_pos
            col_swap = col_tmp(p)
            value_swap = value_tmp(p)
            q = p - 1
            do while (q >= start_pos)
               if (col_tmp(q) <= col_swap) exit
               col_tmp(q + 1) = col_tmp(q)
               value_tmp(q + 1) = value_tmp(q)
               q = q - 1
            end do
            col_tmp(q + 1) = col_swap
            value_tmp(q + 1) = value_swap
         end do
      end do

      allocate(unique_counts(a%nrow))
      unique_counts = 0
      do i = 1, a%nrow
         start_pos = row_ptr(i)
         stop_pos = row_ptr(i + 1) - 1
         if (stop_pos >= start_pos) then
            unique_counts(i) = 1
            do p = start_pos + 1, stop_pos
               if (col_tmp(p) /= col_tmp(p - 1)) unique_counts(i) = unique_counts(i) + 1
            end do
         end if
      end do
      row_ptr(1) = 1
      do i = 1, a%nrow
         row_ptr(i + 1) = row_ptr(i) + unique_counts(i)
      end do
      total_unique = row_ptr(a%nrow + 1) - 1
      allocate(col(total_unique), value(total_unique))
      q = 0
      do i = 1, a%nrow
         start_pos = sum(counts(1:max(0, i - 1))) + 1
         stop_pos = start_pos + counts(i) - 1
         if (stop_pos < start_pos) cycle
         current_col = col_tmp(start_pos)
         current_value = value_tmp(start_pos)
         do p = start_pos + 1, stop_pos
            if (col_tmp(p) == current_col) then
               current_value = current_value + value_tmp(p)
            else
               q = q + 1
               col(q) = current_col
               value(q) = current_value
               current_col = col_tmp(p)
               current_value = value_tmp(p)
            end if
         end do
         q = q + 1
         col(q) = current_col
         value(q) = current_value
      end do
      info = 0
   end subroutine csr_from_sparse

   pure integer function csr_find_index(row_ptr, col, row, target_col) result(position)
      !! Locate a sorted CSR entry by row and column, returning zero when the structural entry is absent.
      integer, intent(in) :: row_ptr(:) !! One-based CSR row pointers.
      integer, intent(in) :: col(:) !! Sorted one-based column indices.
      integer, intent(in) :: row !! One-based row to search.
      integer, intent(in) :: target_col !! One-based column to locate in the selected row.
      integer :: left, middle, right

      position = 0
      if (row < 1 .or. row + 1 > size(row_ptr)) return
      left = row_ptr(row)
      right = row_ptr(row + 1) - 1
      do while (left <= right)
         middle = (left + right)/2
         if (col(middle) == target_col) then
            position = middle
            return
         else if (col(middle) < target_col) then
            left = middle + 1
         else
            right = middle - 1
         end if
      end do
   end function csr_find_index

   pure subroutine csr_transpose(n, row_ptr, col, value, t_row_ptr, t_col, t_value)
      !! Transpose a square CSR matrix into another row-sorted CSR structure.
      integer, intent(in) :: n !! Number of rows and columns in the square matrix.
      integer, intent(in) :: row_ptr(:) !! One-based CSR row pointers for the source matrix.
      integer, intent(in) :: col(:) !! One-based source column indices.
      real(dp), intent(in) :: value(:) !! Source values corresponding to col.
      integer, allocatable, intent(out) :: t_row_ptr(:) !! One-based CSR row pointers for the transpose.
      integer, allocatable, intent(out) :: t_col(:) !! Transposed column indices.
      real(dp), allocatable, intent(out) :: t_value(:) !! Transposed values.
      integer, allocatable :: counts(:), next_slot(:)
      integer :: i, j, p, q

      allocate(counts(n), t_row_ptr(n + 1))
      counts = 0
      do p = 1, size(col)
         if (col(p) >= 1 .and. col(p) <= n) counts(col(p)) = counts(col(p)) + 1
      end do
      t_row_ptr(1) = 1
      do i = 1, n
         t_row_ptr(i + 1) = t_row_ptr(i) + counts(i)
      end do
      allocate(t_col(size(col)), t_value(size(value)), next_slot(n))
      next_slot = t_row_ptr(1:n)
      do i = 1, n
         do p = row_ptr(i), row_ptr(i + 1) - 1
            j = col(p)
            q = next_slot(j)
            t_col(q) = i
            t_value(q) = value(p)
            next_slot(j) = q + 1
         end do
      end do
   end subroutine csr_transpose

   pure subroutine assign_diagnostics(iteration_value, error_value, convergence_value, iterations, error, converged)
      !! Assign optional diagnostic outputs in invalid-input early returns.
      integer, intent(in) :: iteration_value !! Iteration count to return.
      real(dp), intent(in) :: error_value !! Error estimate to return.
      logical, intent(in) :: convergence_value !! Convergence flag to return.
      integer, intent(out), optional :: iterations !! Optional caller-visible iteration count.
      real(dp), intent(out), optional :: error !! Optional caller-visible error estimate.
      logical, intent(out), optional :: converged !! Optional caller-visible convergence flag.

      if (present(iterations)) iterations = iteration_value
      if (present(error)) error = error_value
      if (present(converged)) converged = convergence_value
   end subroutine assign_diagnostics

end module sanic_iterative
