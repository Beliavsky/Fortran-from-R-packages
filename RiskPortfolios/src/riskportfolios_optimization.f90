! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_optimization
   use riskportfolios_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: OBJECTIVE_MEAN_VARIANCE = 1
   integer, parameter, public :: OBJECTIVE_MINIMUM_VARIANCE = 2
   integer, parameter, public :: OBJECTIVE_ERC = 3
   integer, parameter, public :: OBJECTIVE_MAXIMUM_DIVERSIFICATION = 4
   integer, parameter, public :: OBJECTIVE_RISK_EFFICIENT = 5
   integer, parameter, public :: OBJECTIVE_MAXIMUM_DECORRELATION = 6

   public :: projected_gradient, project_feasible

contains

   subroutine projected_gradient(objective, sigma, vector, gamma, x0, lower, &
      upper, use_gross, gross_limit, max_iterations, tolerance, x, info)
      integer, intent(in) :: objective
      real(dp), intent(in) :: sigma(:, :), vector(:), gamma
      real(dp), intent(in) :: x0(:), lower(:), upper(:)
      logical, intent(in) :: use_gross
      real(dp), intent(in) :: gross_limit
      integer, intent(in) :: max_iterations
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp) :: g(size(x0)), g_new(size(x0)), trial(size(x0))
      real(dp) :: s(size(x0)), y(size(x0))
      real(dp) :: f, f_new, alpha, step, denom, projected_norm
      integer :: iter, ls, project_info

      call project_feasible(x0, lower, upper, use_gross, gross_limit, x, project_info)
      if (project_info /= 0) then
         info = project_info
         return
      end if
      call evaluate_objective(objective, sigma, vector, gamma, x, f, g)
      step = 1.0_dp
      info = 1

      do iter = 1, max_iterations
         call project_feasible(x - g, lower, upper, use_gross, gross_limit, &
            trial, project_info)
         projected_norm = sqrt(sum((trial - x) ** 2))
         if (projected_norm <= tolerance * (1.0_dp + sqrt(sum(x * x)))) then
            info = 0
            exit
         end if

         alpha = step
         do ls = 1, 40
            call project_feasible(x - alpha * g, lower, upper, use_gross, &
               gross_limit, trial, project_info)
            call evaluate_objective(objective, sigma, vector, gamma, trial, &
               f_new, g_new)
            if (f_new <= f - 1.0e-4_dp * dot_product(g, x - trial)) exit
            alpha = 0.5_dp * alpha
         end do
         if (ls > 40) then
            info = 2
            exit
         end if

         s = trial - x
         y = g_new - g
         denom = dot_product(s, y)
         if (denom > 1.0e-20_dp) then
            step = max(1.0e-8_dp, min(1.0e3_dp, dot_product(s, s) / denom))
         else
            step = 1.0_dp
         end if
         x = trial
         f = f_new
         g = g_new
         if (sqrt(sum(s * s)) <= tolerance * &
            (1.0_dp + sqrt(sum(x * x)))) then
            info = 0
            exit
         end if
      end do
   end subroutine projected_gradient

   subroutine evaluate_objective(objective, sigma, vector, gamma, w, value, gradient)
      integer, intent(in) :: objective
      real(dp), intent(in) :: sigma(:, :), vector(:), gamma, w(:)
      real(dp), intent(out) :: value, gradient(:)
      real(dp) :: sigma_w(size(w)), residual(size(w)), prc(size(w))
      real(dp) :: vol(size(w)), variance, portfolio_sd, numerator
      real(dp) :: weighted_residual
      integer :: i

      sigma_w = matmul(sigma, w)
      select case (objective)
      case (OBJECTIVE_MEAN_VARIANCE)
         value = -dot_product(vector, w) + 0.5_dp * gamma * dot_product(w, sigma_w)
         gradient = -vector + gamma * sigma_w

      case (OBJECTIVE_MINIMUM_VARIANCE, OBJECTIVE_MAXIMUM_DECORRELATION)
         value = dot_product(w, sigma_w)
         gradient = 2.0_dp * sigma_w

      case (OBJECTIVE_ERC)
         variance = max(dot_product(w, sigma_w), tiny(1.0_dp))
         prc = w * sigma_w / variance
         residual = prc - 1.0_dp / real(size(w), dp)
         value = sum(residual * residual)
         weighted_residual = dot_product(w * residual, sigma_w)
         gradient = 2.0_dp * (variance * (matmul(sigma, w * residual) + &
            residual * sigma_w) - 2.0_dp * sigma_w * weighted_residual) / variance ** 2

      case (OBJECTIVE_MAXIMUM_DIVERSIFICATION)
         do i = 1, size(w)
            vol(i) = sqrt(max(sigma(i, i), tiny(1.0_dp)))
         end do
         portfolio_sd = sqrt(max(dot_product(w, sigma_w), tiny(1.0_dp)))
         numerator = dot_product(w, vol)
         value = -numerator / portfolio_sd
         gradient = -(portfolio_sd * vol - numerator * sigma_w / portfolio_sd) / &
            portfolio_sd ** 2

      case (OBJECTIVE_RISK_EFFICIENT)
         portfolio_sd = sqrt(max(dot_product(w, sigma_w), tiny(1.0_dp)))
         numerator = dot_product(w, vector)
         value = -numerator / portfolio_sd
         gradient = -(portfolio_sd * vector - numerator * sigma_w / portfolio_sd) / &
            portfolio_sd ** 2

      case default
         value = huge(1.0_dp)
         gradient = 0.0_dp
      end select
   end subroutine evaluate_objective

   subroutine project_feasible(y, lower, upper, use_gross, gross_limit, x, info)
      real(dp), intent(in) :: y(:), lower(:), upper(:)
      logical, intent(in) :: use_gross
      real(dp), intent(in) :: gross_limit
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp) :: p(size(y)), q(size(y)), a(size(y)), b(size(y))
      real(dp) :: previous(size(y))
      integer :: iter, stat

      if (size(y) /= size(lower) .or. size(y) /= size(upper)) then
         info = -1
         x = 0.0_dp
         return
      end if
      if (sum(lower) > 1.0_dp + 1.0e-12_dp .or. &
          sum(upper) < 1.0_dp - 1.0e-12_dp) then
         info = -2
         x = 0.0_dp
         return
      end if
      if (.not. use_gross) then
         call project_box_sum(y, lower, upper, x, stat)
         info = stat
         return
      end if
      if (gross_limit < 1.0_dp - 1.0e-12_dp) then
         info = -3
         x = 0.0_dp
         return
      end if

      x = y
      p = 0.0_dp
      q = 0.0_dp
      do iter = 1, 250
         previous = x
         call project_box_sum(x + p, lower, upper, a, stat)
         if (stat /= 0) then
            info = stat
            return
         end if
         p = x + p - a
         call project_l1_ball(a + q, gross_limit, b)
         q = a + q - b
         x = b
         if (sqrt(sum((x - previous) ** 2)) <= 1.0e-13_dp .and. &
             abs(sum(x) - 1.0_dp) <= 1.0e-11_dp .and. &
             sum(abs(x)) <= gross_limit + 1.0e-11_dp) exit
      end do
      call project_box_sum(x, lower, upper, a, stat)
      x = a
      if (sum(abs(x)) > gross_limit + 1.0e-8_dp) then
         info = 3
      else
         info = 0
      end if
   end subroutine project_feasible

   subroutine project_box_sum(y, lower, upper, x, info)
      real(dp), intent(in) :: y(:), lower(:), upper(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp) :: tau_low, tau_high, tau_mid, total
      integer :: iter

      if (any(lower > upper)) then
         info = -4
         x = 0.0_dp
         return
      end if
      tau_low = minval(y - upper)
      tau_high = maxval(y - lower)
      do iter = 1, 200
         tau_mid = 0.5_dp * (tau_low + tau_high)
         x = min(upper, max(lower, y - tau_mid))
         total = sum(x)
         if (total > 1.0_dp) then
            tau_low = tau_mid
         else
            tau_high = tau_mid
         end if
         if (abs(total - 1.0_dp) <= 1.0e-14_dp) then
            info = 0
            return
         end if
      end do
      x = min(upper, max(lower, y - 0.5_dp * (tau_low + tau_high)))
      info = 0
   end subroutine project_box_sum

   subroutine project_l1_ball(y, radius, x)
      real(dp), intent(in) :: y(:), radius
      real(dp), intent(out) :: x(:)
      real(dp) :: u(size(y)), cssv, theta
      integer :: i

      if (sum(abs(y)) <= radius) then
         x = y
         return
      end if
      u = abs(y)
      call sort_descending(u)
      cssv = 0.0_dp
      theta = 0.0_dp
      do i = 1, size(u)
         cssv = cssv + u(i)
         if (u(i) - (cssv - radius) / real(i, dp) > 0.0_dp) then
            theta = (cssv - radius) / real(i, dp)
         end if
      end do
      x = sign(max(abs(y) - theta, 0.0_dp), y)
   end subroutine project_l1_ball

   subroutine sort_descending(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) >= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_descending

end module riskportfolios_optimization
