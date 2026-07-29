! SPDX-License-Identifier: GPL-3.0-only
! Derived from riskParityPortfolio 0.2.2.9000, Copyright Ze Vinicius and Daniel P. Palomar.
module rpp_core
   use rpp_kinds, only: dp
   implicit none
   private
   public :: portfolio_variance, risk_contributions, relative_risk_contributions
   public :: objective_spinu, objective_roncalli, diagonal_risk_parity
   public :: validate_covariance, validate_budgets
contains
   pure real(dp) function portfolio_variance(sigma, w) result(v)
      real(dp), intent(in) :: sigma(:, :), w(:)
      v = dot_product(w, matmul(sigma, w))
   end function portfolio_variance

   pure function risk_contributions(sigma, w) result(rc)
      real(dp), intent(in) :: sigma(:, :), w(:)
      real(dp) :: rc(size(w))
      rc = w * matmul(sigma, w)
   end function risk_contributions

   pure function relative_risk_contributions(sigma, w) result(rrc)
      real(dp), intent(in) :: sigma(:, :), w(:)
      real(dp) :: rrc(size(w)), total
      rrc = risk_contributions(sigma, w)
      total = sum(rrc)
      if (abs(total) > tiny(1.0_dp)) then
         rrc = rrc / total
      else
         rrc = 0.0_dp
      end if
   end function relative_risk_contributions

   pure real(dp) function objective_spinu(sigma, x, b) result(v)
      real(dp), intent(in) :: sigma(:, :), x(:), b(:)
      if (any(x <= 0.0_dp)) then
         v = huge(1.0_dp)
      else
         v = 0.5_dp * portfolio_variance(sigma, x) - sum(b * log(x))
      end if
   end function objective_spinu

   pure real(dp) function objective_roncalli(sigma, x, b) result(v)
      real(dp), intent(in) :: sigma(:, :), x(:), b(:)
      if (any(x <= 0.0_dp)) then
         v = huge(1.0_dp)
      else
         v = sqrt(max(0.0_dp, portfolio_variance(sigma, x))) - sum(b * log(x))
      end if
   end function objective_roncalli

   pure function diagonal_risk_parity(sigma, b) result(w)
      real(dp), intent(in) :: sigma(:, :), b(:)
      real(dp) :: w(size(b))
      integer :: i
      do i = 1, size(b)
         w(i) = sqrt(max(0.0_dp, b(i))) / sqrt(max(tiny(1.0_dp), sigma(i, i)))
      end do
      if (sum(w) > 0.0_dp) w = w / sum(w)
   end function diagonal_risk_parity

   pure logical function validate_covariance(sigma, tol) result(ok)
      real(dp), intent(in) :: sigma(:, :)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      ok = size(sigma, 1) == size(sigma, 2)
      if (.not. ok) return
      ok = maxval(abs(sigma - transpose(sigma))) <= eps * max(1.0_dp, maxval(abs(sigma)))
      if (.not. ok) return
      do i = 1, size(sigma, 1)
         if (sigma(i, i) <= 0.0_dp) then
            ok = .false.
            return
         end if
      end do
   end function validate_covariance

   pure logical function validate_budgets(b) result(ok)
      real(dp), intent(in) :: b(:)
      ok = size(b) > 0 .and. all(b > 0.0_dp) .and. sum(b) > 0.0_dp
   end function validate_budgets
end module rpp_core
