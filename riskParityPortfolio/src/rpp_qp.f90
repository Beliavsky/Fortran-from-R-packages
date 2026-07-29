! SPDX-License-Identifier: GPL-3.0-only
! Self-contained quadratic-programming and projection routines for the Fortran port.
module rpp_qp
   use rpp_kinds, only: dp
   use rpp_linalg, only: solve_linear, vector_norm2, identity_matrix
   use rpp_types, only: RPP_OK, RPP_INFEASIBLE, RPP_LINEAR_SOLVE_FAILED, RPP_MAX_ITER
   implicit none
   private
   public :: solve_equality_qp, solve_qp_active_set
   public :: project_equality, project_feasible, project_linear_constraints
   public :: project_budget_box, is_feasible
contains
   subroutine solve_equality_qp(qmat, qvec, amat, bvec, x, lambda, info)
      real(dp), intent(in) :: qmat(:, :), qvec(:), amat(:, :), bvec(:)
      real(dp), intent(out) :: x(:)
      real(dp), intent(out), optional :: lambda(:)
      integer, intent(out) :: info
      real(dp), allocatable :: kkt(:, :), rhs(:), sol(:)
      integer :: n, m, linfo
      n = size(qvec)
      m = size(amat, 1)
      if (size(qmat, 1) /= n .or. size(qmat, 2) /= n .or. size(amat, 2) /= n .or. &
          size(bvec) /= m .or. size(x) /= n) then
         info = RPP_LINEAR_SOLVE_FAILED
         return
      end if
      if (m == 0) then
         call solve_linear(qmat, -qvec, x, linfo)
         if (linfo == 0) then
            info = RPP_OK
         else
            info = RPP_LINEAR_SOLVE_FAILED
         end if
         return
      end if
      allocate(kkt(n + m, n + m), rhs(n + m), sol(n + m))
      kkt = 0.0_dp
      kkt(1:n, 1:n) = qmat
      kkt(1:n, n + 1:n + m) = transpose(amat)
      kkt(n + 1:n + m, 1:n) = amat
      rhs(1:n) = -qvec
      rhs(n + 1:n + m) = bvec
      call solve_linear(kkt, rhs, sol, linfo)
      if (linfo /= 0) then
         info = RPP_LINEAR_SOLVE_FAILED
         return
      end if
      x = sol(1:n)
      if (present(lambda)) then
         if (size(lambda) == m) lambda = sol(n + 1:n + m)
      end if
      info = RPP_OK
   end subroutine solve_equality_qp

   subroutine project_equality(w, cmat, cvec, projected, info)
      real(dp), intent(in) :: w(:), cmat(:, :), cvec(:)
      real(dp), intent(out) :: projected(:)
      integer, intent(out) :: info
      real(dp), allocatable :: gram(:, :), rhs(:), y(:)
      integer :: m, linfo
      m = size(cmat, 1)
      if (m == 0) then
         projected = w
         info = RPP_OK
         return
      end if
      allocate(gram(m, m), rhs(m), y(m))
      gram = matmul(cmat, transpose(cmat))
      rhs = matmul(cmat, w) - cvec
      call solve_linear(gram, rhs, y, linfo)
      if (linfo /= 0) then
         info = RPP_LINEAR_SOLVE_FAILED
         return
      end if
      projected = w - matmul(transpose(cmat), y)
      info = RPP_OK
   end subroutine project_equality

   subroutine project_feasible(w0, cmat, cvec, dmat, dvec, w, info, tol, maxiter)
      real(dp), intent(in) :: w0(:), cmat(:, :), cvec(:), dmat(:, :), dvec(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: direction(:), projected_direction(:), trial(:), zero_rhs(:)
      real(dp) :: eps, violation, denom, max_violation
      integer :: it, j, mit, pinfo
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      mit = 20000
      if (present(maxiter)) mit = maxiter
      call project_equality(w0, cmat, cvec, w, pinfo)
      if (pinfo /= RPP_OK) then
         info = pinfo
         return
      end if
      allocate(direction(size(w)), projected_direction(size(w)), trial(size(w)), zero_rhs(size(cmat, 1)))
      zero_rhs = 0.0_dp
      do it = 1, mit
         max_violation = 0.0_dp
         do j = 1, size(dmat, 1)
            violation = dot_product(dmat(j, :), w) - dvec(j)
            max_violation = max(max_violation, violation)
            if (violation > eps) then
               direction = dmat(j, :)
               if (size(cmat, 1) > 0) then
                  call project_equality(direction, cmat, zero_rhs, projected_direction, pinfo)
                  if (pinfo /= RPP_OK) then
                     info = pinfo
                     return
                  end if
               else
                  projected_direction = direction
               end if
               denom = dot_product(direction, projected_direction)
               if (denom <= tiny(1.0_dp)) then
                  info = RPP_INFEASIBLE
                  return
               end if
               trial = w - violation * projected_direction / denom
               call project_equality(trial, cmat, cvec, w, pinfo)
               if (pinfo /= RPP_OK) then
                  info = pinfo
                  return
               end if
            end if
         end do
         if (max_violation <= eps .and. is_feasible(w, cmat, cvec, dmat, dvec, 10.0_dp * eps)) then
            info = RPP_OK
            return
         end if
      end do
      info = RPP_INFEASIBLE
   end subroutine project_feasible


   subroutine project_linear_constraints(w0, cmat, cvec, dmat, dvec, w, info, tol, maxiter)
      real(dp), intent(in) :: w0(:), cmat(:, :), cvec(:), dmat(:, :), dvec(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: feasible(:), qmat(:, :)
      real(dp) :: eps
      integer :: mit, qiter
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      mit = 2000
      if (present(maxiter)) mit = maxiter
      allocate(feasible(size(w0)), qmat(size(w0), size(w0)))
      call project_feasible(w0, cmat, cvec, dmat, dvec, feasible, info, eps, 20 * mit)
      if (info /= RPP_OK) return
      qmat = identity_matrix(size(w0))
      call solve_qp_active_set(qmat, -w0, cmat, cvec, dmat, dvec, feasible, w, info, qiter, eps, mit)
   end subroutine project_linear_constraints

   subroutine project_budget_box(w0, lower, upper, w, info, tol)
      real(dp), intent(in) :: w0(:), lower(:), upper(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp) :: lo, hi, mid, value, eps
      integer :: it
      eps = 1.0e-12_dp
      if (present(tol)) eps = tol
      if (sum(lower) > 1.0_dp + eps .or. sum(upper) < 1.0_dp - eps .or. any(lower > upper)) then
         info = RPP_INFEASIBLE
         return
      end if
      lo = minval(w0 - upper)
      hi = maxval(w0 - lower)
      do it = 1, 200
         mid = 0.5_dp * (lo + hi)
         w = min(max(w0 - mid, lower), upper)
         value = sum(w) - 1.0_dp
         if (abs(value) <= eps) exit
         if (value > 0.0_dp) then
            lo = mid
         else
            hi = mid
         end if
      end do
      w = min(max(w0 - 0.5_dp * (lo + hi), lower), upper)
      info = RPP_OK
   end subroutine project_budget_box

   pure logical function is_feasible(w, cmat, cvec, dmat, dvec, tol) result(ok)
      real(dp), intent(in) :: w(:), cmat(:, :), cvec(:), dmat(:, :), dvec(:), tol
      ok = .true.
      if (size(cmat, 1) > 0) ok = maxval(abs(matmul(cmat, w) - cvec)) <= tol
      if (ok .and. size(dmat, 1) > 0) ok = maxval(matmul(dmat, w) - dvec) <= tol
   end function is_feasible

   subroutine solve_qp_active_set(qmat, qvec, cmat, cvec, dmat, dvec, x0, x, info, &
                                  iterations, tol, maxiter)
      real(dp), intent(in) :: qmat(:, :), qvec(:), cmat(:, :), cvec(:)
      real(dp), intent(in) :: dmat(:, :), dvec(:), x0(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      logical, allocatable :: active(:)
      real(dp), allocatable :: aw(:, :), bw(:), grad(:), kkt(:, :), rhs(:), sol(:), p(:), lambda(:)
      real(dp) :: eps, alpha, candidate, min_mu
      integer :: n, meq, mineq, mw, mit, it, j, blocker, remove_j, linfo, pos

      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      mit = 1000
      if (present(maxiter)) mit = maxiter
      n = size(qvec)
      meq = size(cmat, 1)
      mineq = size(dmat, 1)
      x = x0
      if (.not. is_feasible(x, cmat, cvec, dmat, dvec, 100.0_dp * eps)) then
         call project_feasible(x0, cmat, cvec, dmat, dvec, x, linfo, eps, 20000)
         if (linfo /= RPP_OK) then
            info = linfo
            iterations = 0
            return
         end if
      end if
      allocate(active(mineq))
      active = .false.
      do j = 1, mineq
         active(j) = abs(dot_product(dmat(j, :), x) - dvec(j)) <= 100.0_dp * eps
      end do

      do it = 1, mit
         mw = meq + count(active)
         allocate(aw(mw, n), bw(mw), grad(n), kkt(n + mw, n + mw), rhs(n + mw), &
                  sol(n + mw), p(n), lambda(mw))
         if (meq > 0) then
            aw(1:meq, :) = cmat
            bw(1:meq) = 0.0_dp
         end if
         pos = meq
         do j = 1, mineq
            if (active(j)) then
               pos = pos + 1
               aw(pos, :) = dmat(j, :)
               bw(pos) = 0.0_dp
            end if
         end do
         grad = matmul(qmat, x) + qvec
         if (mw == 0) then
            call solve_linear(qmat, -grad, p, linfo)
            lambda = 0.0_dp
         else
            kkt = 0.0_dp
            kkt(1:n, 1:n) = qmat
            kkt(1:n, n + 1:n + mw) = transpose(aw)
            kkt(n + 1:n + mw, 1:n) = aw
            rhs(1:n) = -grad
            rhs(n + 1:n + mw) = 0.0_dp
            call solve_linear(kkt, rhs, sol, linfo)
            if (linfo == 0) then
               p = sol(1:n)
               lambda = sol(n + 1:n + mw)
            end if
         end if
         if (linfo /= 0) then
            info = RPP_LINEAR_SOLVE_FAILED
            iterations = it
            return
         end if

         if (vector_norm2(p) <= eps * max(1.0_dp, vector_norm2(x))) then
            min_mu = 0.0_dp
            remove_j = 0
            pos = meq
            do j = 1, mineq
               if (active(j)) then
                  pos = pos + 1
                  if (lambda(pos) < min_mu) then
                     min_mu = lambda(pos)
                     remove_j = j
                  end if
               end if
            end do
            if (remove_j == 0 .or. min_mu >= -sqrt(eps)) then
               info = RPP_OK
               iterations = it
               return
            end if
            active(remove_j) = .false.
         else
            alpha = 1.0_dp
            blocker = 0
            do j = 1, mineq
               if (.not. active(j)) then
                  candidate = dot_product(dmat(j, :), p)
                  if (candidate > eps) then
                     candidate = (dvec(j) - dot_product(dmat(j, :), x)) / candidate
                     if (candidate < alpha) then
                        alpha = max(0.0_dp, candidate)
                        blocker = j
                     end if
                  end if
               end if
            end do
            x = x + alpha * p
            if (blocker > 0 .and. alpha < 1.0_dp - eps) active(blocker) = .true.
         end if
         deallocate(aw, bw, grad, kkt, rhs, sol, p, lambda)
      end do
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine solve_qp_active_set
end module rpp_qp
