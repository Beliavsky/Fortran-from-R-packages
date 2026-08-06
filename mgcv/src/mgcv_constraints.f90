! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_constraints
   use mgcv_kinds, only : dp
   use mgcv_linalg, only : spd_solve, solve_linear
   implicit none
   private
   public :: pcls_fit, monotonicity_constraints, convexity_constraints
   public :: tri_cholesky, band_cholesky

contains

   subroutine pcls_fit(x, y, weights, penalties, lambda, beta, status, &
                       a_ineq, b_ineq, a_eq, b_eq, max_iter, tolerance)
      real(dp), intent(in) :: x(:, :), y(:), weights(:)
      real(dp), intent(in), optional :: penalties(:, :, :), lambda(:)
      real(dp), allocatable, intent(out) :: beta(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: a_ineq(:, :), b_ineq(:), a_eq(:, :), b_eq(:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: h(:, :), rhs(:, :), sol(:, :), wx(:, :), pmat(:, :), old(:), grad(:)
      real(dp) :: step, violation, norm2, tol
      integer :: i, j, p, niter, it

      status = 0; p = size(x, 2)
      if (size(x, 1) /= size(y) .or. size(weights) /= size(y) .or. any(weights < 0.0_dp)) then
         allocate(beta(0)); status = 1; return
      end if
      allocate(wx(size(x, 1), p)); wx = x * spread(weights, 2, p)
      allocate(h(p, p), rhs(p, 1), pmat(p, p)); pmat = 0.0_dp
      h = matmul(transpose(x), wx)
      rhs(:, 1) = matmul(transpose(x), weights * y)
      if (present(penalties)) then
         if (.not. present(lambda)) then; allocate(beta(0)); status = 2; return; end if
         if (size(lambda) /= size(penalties, 3)) then; allocate(beta(0)); status = 3; return; end if
         do j = 1, size(lambda)
            pmat = pmat + lambda(j) * penalties(:, :, j)
         end do
      end if
      h = h + pmat
      do i = 1, p; h(i, i) = h(i, i) + 1.0e-10_dp; end do
      call spd_solve(h, rhs, sol, status, 1.0e-12_dp)
      if (status /= 0) then; allocate(beta(0)); return; end if
      allocate(beta(p)); beta = sol(:, 1)
      call project_equalities(beta, a_eq, b_eq, status)
      if (status /= 0) return
      call project_inequalities(beta, a_ineq, b_ineq)

      niter = 2000; if (present(max_iter)) niter = max_iter
      tol = 1.0e-8_dp; if (present(tolerance)) tol = tolerance
      step = 1.0_dp / max(1.0_dp, maxval(sum(abs(h), dim=2)))
      allocate(old(p), grad(p))
      do it = 1, niter
         old = beta
         grad = matmul(h, beta) - rhs(:, 1)
         beta = beta - step * grad
         call project_equalities(beta, a_eq, b_eq, status)
         if (status /= 0) return
         call project_inequalities(beta, a_ineq, b_ineq)
         if (maxval(abs(beta - old)) <= tol * (1.0_dp + maxval(abs(beta)))) exit
      end do
      if (it > niter) status = 4
      if (present(a_ineq) .and. present(b_ineq)) then
         do i = 1, size(a_ineq, 1)
            violation = b_ineq(i) - dot_product(a_ineq(i, :), beta)
            if (violation > 10.0_dp * tol) status = 5
         end do
      end if
      if (present(a_eq) .and. present(b_eq)) then
         norm2 = maxval(abs(matmul(a_eq, beta) - b_eq))
         if (norm2 > 10.0_dp * tol) status = 6
      end if
   end subroutine pcls_fit

   subroutine project_equalities(beta, a_eq, b_eq, status)
      real(dp), intent(inout) :: beta(:)
      real(dp), intent(in), optional :: a_eq(:, :), b_eq(:)
      integer, intent(out) :: status
      real(dp), allocatable :: gram(:, :), rhs(:, :), alpha(:, :)
      status = 0
      if (.not. present(a_eq)) return
      if (.not. present(b_eq) .or. size(a_eq, 1) /= size(b_eq) .or. size(a_eq, 2) /= size(beta)) then
         status = 1; return
      end if
      if (size(a_eq, 1) == 0) return
      gram = matmul(a_eq, transpose(a_eq))
      allocate(rhs(size(b_eq), 1)); rhs(:, 1) = b_eq - matmul(a_eq, beta)
      call solve_linear(gram, rhs, alpha, status)
      if (status == 0) beta = beta + matmul(transpose(a_eq), alpha(:, 1))
   end subroutine project_equalities

   subroutine project_inequalities(beta, a_ineq, b_ineq)
      real(dp), intent(inout) :: beta(:)
      real(dp), intent(in), optional :: a_ineq(:, :), b_ineq(:)
      real(dp) :: violation, norm2
      integer :: i, sweep
      if (.not. present(a_ineq) .or. .not. present(b_ineq)) return
      if (size(a_ineq, 1) /= size(b_ineq) .or. size(a_ineq, 2) /= size(beta)) return
      do sweep = 1, 20
         do i = 1, size(a_ineq, 1)
            violation = b_ineq(i) - dot_product(a_ineq(i, :), beta)
            if (violation > 0.0_dp) then
               norm2 = dot_product(a_ineq(i, :), a_ineq(i, :))
               if (norm2 > 0.0_dp) beta = beta + violation * a_ineq(i, :) / norm2
            end if
         end do
      end do
   end subroutine project_inequalities

   subroutine monotonicity_constraints(n, a, b, increasing)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: a(:, :), b(:)
      logical, intent(in), optional :: increasing
      logical :: up
      integer :: i
      up = .true.; if (present(increasing)) up = increasing
      if (n < 2) then; allocate(a(0, n), b(0)); return; end if
      allocate(a(n - 1, n), b(n - 1)); a = 0.0_dp; b = 0.0_dp
      do i = 1, n - 1
         a(i, i) = merge(-1.0_dp, 1.0_dp, up)
         a(i, i + 1) = merge(1.0_dp, -1.0_dp, up)
      end do
   end subroutine monotonicity_constraints

   subroutine convexity_constraints(n, a, b, convex)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: a(:, :), b(:)
      logical, intent(in), optional :: convex
      logical :: cvx
      integer :: i
      cvx = .true.; if (present(convex)) cvx = convex
      if (n < 3) then; allocate(a(0, n), b(0)); return; end if
      allocate(a(n - 2, n), b(n - 2)); a = 0.0_dp; b = 0.0_dp
      do i = 1, n - 2
         if (cvx) then
            a(i, i:i + 2) = [1.0_dp, -2.0_dp, 1.0_dp]
         else
            a(i, i:i + 2) = [-1.0_dp, 2.0_dp, -1.0_dp]
         end if
      end do
   end subroutine convexity_constraints

   subroutine tri_cholesky(main_diag, off_diag, r_diag, r_off, status)
      real(dp), intent(in) :: main_diag(:), off_diag(:)
      real(dp), allocatable, intent(out) :: r_diag(:), r_off(:)
      integer, intent(out) :: status
      integer :: i, n
      n = size(main_diag)
      if (n < 2 .or. size(off_diag) /= n - 1) then
         allocate(r_diag(0), r_off(0)); status = 1; return
      end if
      allocate(r_diag(n), r_off(n - 1)); status = 0
      if (main_diag(1) <= 0.0_dp) then; status = 2; return; end if
      r_diag(1) = sqrt(main_diag(1))
      do i = 1, n - 1
         r_off(i) = off_diag(i) / r_diag(i)
         if (main_diag(i + 1) - r_off(i)**2 <= 0.0_dp) then; status = 2; return; end if
         r_diag(i + 1) = sqrt(main_diag(i + 1) - r_off(i)**2)
      end do
   end subroutine tri_cholesky

   subroutine band_cholesky(band, rband, status)
      real(dp), intent(in) :: band(:, :)
      real(dp), allocatable, intent(out) :: rband(:, :)
      integer, intent(out) :: status
      real(dp) :: s
      integer :: n, k, j, i, l, low
      k = size(band, 1); n = size(band, 2)
      if (k < 1 .or. n < 1) then; allocate(rband(0, 0)); status = 1; return; end if
      allocate(rband(k, n)); rband = 0.0_dp; status = 0
      do j = 1, n
         s = band(1, j)
         low = max(1, j - k + 1)
         do l = low, j - 1
            s = s - rband(j - l + 1, l)**2
         end do
         if (s <= 0.0_dp) then; status = 2; return; end if
         rband(1, j) = sqrt(s)
         do i = j + 1, min(n, j + k - 1)
            s = band(i - j + 1, j)
            do l = max(1, i - k + 1, j - k + 1), j - 1
               s = s - rband(j - l + 1, l) * rband(i - l + 1, l)
            end do
            rband(i - j + 1, j) = s / rband(1, j)
         end do
      end do
   end subroutine band_cholesky

end module mgcv_constraints
