! SPDX-License-Identifier: GPL-3.0-only
! Fast risk-parity solvers derived from riskParityPortfolio 0.2.2.9000.
module rpp_solvers
   use rpp_kinds, only: dp
   use rpp_core, only: relative_risk_contributions, objective_spinu, objective_roncalli, portfolio_variance
   use rpp_linalg, only: solve_spd, max_abs
   use rpp_types, only: RPP_OK, RPP_INVALID_INPUT, RPP_LINEAR_SOLVE_FAILED, RPP_MAX_ITER
   implicit none
   private
   public :: risk_parity_ccd_spinu, risk_parity_ccd_roncalli
   public :: risk_parity_ccd_choi, risk_parity_newton
   public :: active_risk_parity_ccd
contains
   pure real(dp) function positive_quadratic_root(aux, a, cterm) result(root)
      real(dp), intent(in) :: aux, a, cterm
      real(dp) :: disc
      disc = sqrt(aux * aux + 4.0_dp * a * cterm)
      if (aux >= 0.0_dp) then
         root = (aux + disc) / (2.0_dp * a)
      else
         root = 2.0_dp * cterm / (disc - aux)
      end if
   end function positive_quadratic_root

   subroutine risk_parity_ccd_spinu(sigma, b, w, info, iterations, tol, maxiter)
      real(dp), intent(in) :: sigma(:, :), b(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: x(:), sx(:), rrc(:)
      real(dp) :: eps, aux, xnew, diff, s
      integer :: n, i, k, mit
      n = size(b)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      mit = 1000
      if (present(maxiter)) mit = maxiter
      if (size(w) /= n .or. size(sigma, 1) /= n .or. size(sigma, 2) /= n .or. &
          any(b <= 0.0_dp) .or. any([(sigma(i,i), i=1,n)] <= 0.0_dp)) then
         info = RPP_INVALID_INPUT
         iterations = 0
         return
      end if
      allocate(x(n), sx(n), rrc(n))
      s = sum(sigma)
      if (s <= 0.0_dp) s = sum([(sigma(i, i), i=1,n)])
      x = sqrt(sum(b) / s)
      sx = matmul(sigma, x)
      do k = 1, mit
         do i = 1, n
            aux = x(i) * sigma(i, i) - sx(i)
            xnew = positive_quadratic_root(aux, sigma(i, i), b(i))
            diff = xnew - x(i)
            sx = sx + sigma(:, i) * diff
            x(i) = xnew
         end do
         w = x / sum(x)
         rrc = relative_risk_contributions(sigma, w)
         if (max_abs(rrc - b / sum(b)) < eps) then
            info = RPP_OK
            iterations = k
            return
         end if
      end do
      w = x / sum(x)
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine risk_parity_ccd_spinu

   subroutine risk_parity_ccd_roncalli(sigma, b, w, info, iterations, tol, maxiter)
      real(dp), intent(in) :: sigma(:, :), b(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: x(:), sx(:), rrc(:)
      real(dp) :: eps, aux, xnew, diff, scale, sig, cross
      integer :: n, i, k, mit
      n = size(b)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      mit = 1000
      if (present(maxiter)) mit = maxiter
      if (size(w) /= n .or. size(sigma, 1) /= n .or. size(sigma, 2) /= n .or. any(b <= 0.0_dp)) then
         info = RPP_INVALID_INPUT
         iterations = 0
         return
      end if
      allocate(x(n), sx(n), rrc(n))
      scale = sum(sigma)
      if (scale <= 0.0_dp) scale = sum([(sigma(i, i), i=1,n)])
      x = sqrt(sum(b) / scale)
      sx = matmul(sigma, x)
      sig = sqrt(max(tiny(1.0_dp), dot_product(x, sx)))
      do k = 1, mit
         do i = 1, n
            aux = x(i) * sigma(i, i) - sx(i)
            xnew = positive_quadratic_root(aux, sigma(i, i), b(i) * sig)
            diff = xnew - x(i)
            cross = dot_product(sigma(i, :), x)
            sx = sx + sigma(:, i) * diff
            sig = sqrt(max(tiny(1.0_dp), sig * sig + 2.0_dp * diff * cross + sigma(i, i) * diff * diff))
            x(i) = xnew
         end do
         w = x / sum(x)
         rrc = relative_risk_contributions(sigma, w)
         if (max_abs(rrc - b / sum(b)) < eps) then
            info = RPP_OK
            iterations = k
            return
         end if
      end do
      w = x / sum(x)
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine risk_parity_ccd_roncalli

   subroutine active_risk_parity_ccd(sigma, b, mu, tradeoff, risk_free, w, info, iterations, tol, maxiter)
      real(dp), intent(in) :: sigma(:, :), b(:), mu(:), tradeoff, risk_free
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: x(:), sx(:), pi(:), rrc(:)
      real(dp) :: eps, aux, xnew, diff, scale, sig, cross
      integer :: n, i, k, mit
      n = size(b)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      mit = 1000
      if (present(maxiter)) mit = maxiter
      if (size(w) /= n .or. size(mu) /= n .or. tradeoff <= 0.0_dp) then
         info = RPP_INVALID_INPUT
         iterations = 0
         return
      end if
      allocate(x(n), sx(n), pi(n), rrc(n))
      scale = sum(sigma)
      if (scale <= 0.0_dp) scale = sum([(sigma(i, i), i=1,n)])
      x = sqrt(sum(b) / scale)
      sx = matmul(sigma, x)
      pi = mu - risk_free
      sig = sqrt(max(tiny(1.0_dp), dot_product(x, sx)))
      do k = 1, mit
         do i = 1, n
            aux = tradeoff * (x(i) * sigma(i, i) - sx(i)) + pi(i) * sig
            xnew = (aux + sqrt(aux * aux + 4.0_dp * tradeoff * sigma(i, i) * b(i) * sig)) / &
                   (2.0_dp * tradeoff * sigma(i, i))
            diff = xnew - x(i)
            cross = dot_product(sigma(i, :), x)
            sx = sx + sigma(:, i) * diff
            sig = sqrt(max(tiny(1.0_dp), sig * sig + 2.0_dp * diff * cross + sigma(i, i) * diff * diff))
            x(i) = xnew
         end do
         w = x / sum(x)
         rrc = relative_risk_contributions(sigma, w)
         if (max_abs(rrc - b / sum(b)) < eps) then
            info = RPP_OK
            iterations = k
            return
         end if
      end do
      w = x / sum(x)
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine active_risk_parity_ccd

   subroutine risk_parity_ccd_choi(sigma, b, w, info, iterations, tol, maxiter)
      real(dp), intent(in) :: sigma(:, :), b(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: vol(:), invvol(:), corr(:, :), adj(:, :), x(:), rc(:)
      real(dp) :: eps, ai, scale
      integer :: n, i, k, mit
      n = size(b)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      mit = 1000
      if (present(maxiter)) mit = maxiter
      if (size(w) /= n .or. any([(sigma(i,i), i=1,n)] <= 0.0_dp)) then
         info = RPP_INVALID_INPUT
         iterations = 0
         return
      end if
      allocate(vol(n), invvol(n), corr(n, n), adj(n, n), x(n), rc(n))
      vol = sqrt([(sigma(i, i), i=1,n)])
      invvol = 1.0_dp / vol
      do i = 1, n
         corr(:, i) = sigma(:, i) * invvol * invvol(i)
      end do
      adj = corr
      do i = 1, n
         adj(i, i) = 0.0_dp
      end do
      scale = sum(corr)
      if (scale <= 0.0_dp) scale = real(n, dp)
      x = 1.0_dp / sqrt(scale)
      do k = 1, mit
         do i = 1, n
            ai = 0.5_dp * dot_product(adj(:, i), x)
            x(i) = sqrt(ai * ai + b(i)) - ai
         end do
         x = x / sqrt(dot_product(x, matmul(corr, x)))
         rc = x * matmul(corr, x)
         if (max_abs(rc - b) < eps) then
            w = x / vol
            w = w / sum(w)
            info = RPP_OK
            iterations = k
            return
         end if
      end do
      w = x / vol
      w = w / sum(w)
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine risk_parity_ccd_choi

   subroutine risk_parity_newton(sigma, b, w, info, iterations, tol, maxiter)
      real(dp), intent(in) :: sigma(:, :), b(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: x(:), grad(:), hess(:, :), direction(:), rrc(:), trial(:)
      real(dp) :: eps, scale, decrement, dx, step, f0, f1
      integer :: n, i, k, mit, linfo
      n = size(b)
      eps = 1.0e-8_dp
      if (present(tol)) eps = tol
      mit = 100
      if (present(maxiter)) mit = maxiter
      allocate(x(n), grad(n), hess(n, n), direction(n), rrc(n), trial(n))
      scale = sum(sigma)
      if (scale <= 0.0_dp) scale = sum([(sigma(i, i), i=1,n)])
      x = sqrt(sum(b) / scale)
      do k = 1, mit
         grad = matmul(sigma, x) - b / x
         hess = sigma
         do i = 1, n
            hess(i, i) = hess(i, i) + b(i) / (x(i) * x(i))
         end do
         call solve_spd(hess, grad, direction, linfo)
         if (linfo /= 0) then
            info = RPP_LINEAR_SOLVE_FAILED
            iterations = k
            return
         end if
         decrement = sqrt(max(0.0_dp, dot_product(grad, direction)))
         dx = maxval(direction / x)
         step = 1.0_dp / (1.0_dp + max(0.0_dp, dx))
         f0 = objective_spinu(sigma, x, b)
         do
            trial = x - step * direction
            if (all(trial > 0.0_dp)) then
               f1 = objective_spinu(sigma, trial, b)
               if (f1 <= f0 - 1.0e-4_dp * step * dot_product(grad, direction)) exit
            end if
            step = 0.5_dp * step
            if (step < 1.0e-14_dp) exit
         end do
         x = trial
         w = x / sum(x)
         rrc = relative_risk_contributions(sigma, w)
         if (max_abs(rrc - b / sum(b)) <= 2.0_dp * eps .or. decrement < eps) then
            info = RPP_OK
            iterations = k
            return
         end if
      end do
      w = x / sum(x)
      info = RPP_MAX_ITER
      iterations = mit
   end subroutine risk_parity_newton
end module rpp_solvers
