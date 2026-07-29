! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_risk
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use gcpm_kinds, only: dp
   use gcpm_types, only: credit_portfolio, gcpm_model, loss_distribution, &
      risk_measures, model_analytical_crp, model_simulation
   implicit none
   private

   public :: distribution_index
   public :: value_at_risk
   public :: expected_shortfall
   public :: economic_capital
   public :: calculate_risk_measures
   public :: search_tau_for_target_es
   public :: analytical_risk_contributions
   public :: simulation_risk_contributions

contains

   integer function distribution_index(distribution, alpha) result(index_value)
      type(loss_distribution), intent(in) :: distribution
      real(dp), intent(in) :: alpha
      integer :: i

      index_value = 0
      if (.not. allocated(distribution%cdf)) return
      do i = 1, size(distribution%cdf)
         if (distribution%cdf(i) >= alpha - 1.0e-12_dp) then
            index_value = i
            return
         end if
      end do
   end function distribution_index

   function value_at_risk(distribution, alpha, status) result(var_value)
      type(loss_distribution), intent(in) :: distribution
      real(dp), intent(in) :: alpha
      integer, intent(out), optional :: status
      real(dp) :: var_value
      integer :: idx

      if (present(status)) status = 0
      var_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         if (present(status)) status = 1
         return
      end if
      idx = distribution_index(distribution, alpha)
      if (idx == 0) then
         if (present(status)) status = 2
         return
      end if
      var_value = distribution%loss(idx)
   end function value_at_risk

   function expected_shortfall(distribution, portfolio, alpha, status) result(es_value)
      type(loss_distribution), intent(in) :: distribution
      type(credit_portfolio), intent(in) :: portfolio
      real(dp), intent(in) :: alpha
      integer, intent(out), optional :: status
      real(dp) :: es_value
      real(dp) :: tail_probability, numerator
      integer :: idx

      if (present(status)) status = 0
      es_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         if (present(status)) status = 1
         return
      end if
      idx = distribution_index(distribution, alpha)
      if (idx == 0) then
         if (present(status)) status = 2
         return
      end if
      if (idx == 1) then
         es_value = portfolio%analytical_expected_loss
         return
      end if

      tail_probability = 1.0_dp - distribution%cdf(idx - 1)
      if (tail_probability <= tiny(1.0_dp)) then
         if (present(status)) status = 3
         return
      end if

      if (distribution%model_kind == model_analytical_crp) then
         numerator = portfolio%analytical_expected_loss - &
            sum(distribution%pdf(1:idx - 1) * distribution%loss(1:idx - 1))
      else
         numerator = sum(distribution%pdf(idx:) * distribution%loss(idx:))
      end if
      es_value = numerator / tail_probability
   end function expected_shortfall

   function economic_capital(distribution, portfolio, alpha, status) result(ec_value)
      type(loss_distribution), intent(in) :: distribution
      type(credit_portfolio), intent(in) :: portfolio
      real(dp), intent(in) :: alpha
      integer, intent(out), optional :: status
      real(dp) :: ec_value
      integer :: local_status

      ec_value = value_at_risk(distribution, alpha, local_status) - &
         portfolio%analytical_expected_loss
      if (present(status)) status = local_status
   end function economic_capital

   subroutine calculate_risk_measures(distribution, portfolio, alpha, measures, status, message)
      type(loss_distribution), intent(in) :: distribution
      type(credit_portfolio), intent(in) :: portfolio
      real(dp), intent(in) :: alpha(:)
      type(risk_measures), intent(out) :: measures
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      integer :: i, local_status

      status = 0
      message = ''
      allocate(measures%alpha(size(alpha)), measures%var(size(alpha)), &
               measures%economic_capital(size(alpha)), &
               measures%expected_shortfall(size(alpha)))
      measures%alpha = alpha

      do i = 1, size(alpha)
         measures%var(i) = value_at_risk(distribution, alpha(i), local_status)
         if (local_status /= 0) then
            status = local_status
            message = 'one or more requested VaR levels are unavailable'
         end if
         measures%economic_capital(i) = measures%var(i) - &
            portfolio%analytical_expected_loss
         measures%expected_shortfall(i) = expected_shortfall(distribution, portfolio, &
                                                              alpha(i), local_status)
         if (local_status /= 0) then
            status = local_status
            message = 'one or more requested expected-shortfall levels are unavailable'
         end if
      end do
   end subroutine calculate_risk_measures

   function search_tau_for_target_es(distribution, portfolio, target_loss, status) result(tau)
      type(loss_distribution), intent(in) :: distribution
      type(credit_portfolio), intent(in) :: portfolio
      real(dp), intent(in) :: target_loss
      integer, intent(out), optional :: status
      real(dp) :: tau
      real(dp) :: candidate_alpha, candidate_es, best_error, current_error
      integer :: i, local_status, best_i, upper

      if (present(status)) status = 0
      tau = ieee_value(0.0_dp, ieee_quiet_nan)
      if (target_loss > distribution%loss(size(distribution%loss))) then
         if (present(status)) status = 1
         return
      end if

      upper = size(distribution%cdf)
      if (distribution%cdf(upper) >= 1.0_dp - 1.0e-12_dp) upper = max(1, upper - 1)
      best_error = huge(1.0_dp)
      best_i = 1
      do i = 1, upper
         candidate_alpha = min(1.0_dp - 1.0e-12_dp, &
                               max(1.0e-12_dp, distribution%cdf(i)))
         candidate_es = expected_shortfall(distribution, portfolio, candidate_alpha, local_status)
         if (local_status /= 0) cycle
         current_error = abs(candidate_es - target_loss)
         if (current_error < best_error) then
            best_error = current_error
            best_i = i
         end if
      end do
      tau = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, distribution%cdf(best_i)))
   end function search_tau_for_target_es

   subroutine analytical_risk_contributions(portfolio, model, distribution, alpha, &
                                             contributions, status, message)
      type(credit_portfolio), intent(in) :: portfolio
      type(gcpm_model), intent(in) :: model
      type(loss_distribution), intent(in) :: distribution
      real(dp), intent(in) :: alpha
      real(dp), allocatable, intent(out) :: contributions(:,:)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message

      integer :: i, h, j, r, ec_start, var_start, es_start
      integer :: k_index
      real(dp) :: var_value, es_value, tau, denominator
      real(dp), allocatable :: distorted_pdf(:,:), distorted_cdf(:,:), full_weight(:)
      real(dp), allocatable :: coefficient(:)
      integer :: local_status, n_support, n_factors

      status = 0
      message = ''
      allocate(contributions(portfolio%n_counterparties, 3))
      contributions = ieee_value(0.0_dp, ieee_quiet_nan)
      if (model%model_kind /= model_analytical_crp .or. &
          distribution%model_kind /= model_analytical_crp) then
         status = 1
         message = 'analytical contributions require an analytical CreditRisk+ model and distribution'
         return
      end if
      if (.not. allocated(distribution%recursion_a) .or. &
          .not. allocated(distribution%recursion_b)) then
         status = 2
         message = 'analytical recursion coefficients are unavailable'
         return
      end if

      var_value = value_at_risk(distribution, alpha, local_status)
      if (local_status /= 0) then
         status = 3
         message = 'requested alpha is unavailable'
         return
      end if
      es_value = expected_shortfall(distribution, portfolio, alpha, local_status)
      if (local_status /= 0 .or. es_value <= 0.0_dp) then
         status = 4
         message = 'expected shortfall is unavailable'
         return
      end if

      tau = search_tau_for_target_es(distribution, portfolio, &
                                     var_value - portfolio%analytical_expected_loss, local_status)
      if (local_status /= 0) then
         status = 5
         message = 'economic-capital contribution threshold is unavailable'
         return
      end if
      ec_start = distribution_index(distribution, tau)

      tau = search_tau_for_target_es(distribution, portfolio, var_value, local_status)
      if (local_status /= 0) then
         status = 6
         message = 'VaR contribution threshold is unavailable'
         return
      end if
      var_start = distribution_index(distribution, tau)
      es_start = distribution_index(distribution, alpha)

      n_support = size(distribution%pdf)
      n_factors = portfolio%n_sectors + 1
      allocate(distorted_pdf(n_factors, n_support), distorted_cdf(n_factors, n_support))
      allocate(coefficient(n_support), full_weight(n_factors))
      distorted_pdf = 0.0_dp
      distorted_pdf(1, :) = distribution%pdf

      do h = 1, portfolio%n_sectors
         distorted_pdf(h + 1, 1) = exp(distribution%recursion_a(1) + &
                                         distribution%recursion_b(h, 1))
         coefficient = distribution%recursion_a + distribution%recursion_b(h, :)
         do j = 2, n_support
            distorted_pdf(h + 1, j) = 0.0_dp
            do r = 1, j - 1
               distorted_pdf(h + 1, j) = distorted_pdf(h + 1, j) + &
                  real(r, dp) / real(j - 1, dp) * coefficient(r + 1) * &
                  distorted_pdf(h + 1, j - r)
            end do
         end do
      end do

      do h = 1, n_factors
         distorted_cdf(h, 1) = distorted_pdf(h, 1)
         do j = 2, n_support
            distorted_cdf(h, j) = distorted_cdf(h, j - 1) + distorted_pdf(h, j)
         end do
      end do

      do i = 1, portfolio%n_counterparties
         full_weight(1) = portfolio%idiosyncratic(i)
         full_weight(2:) = portfolio%weight(i, :)

         k_index = ec_start - portfolio%loss_multiple(i)
         denominator = 1.0_dp
         if (ec_start > 1) denominator = 1.0_dp - distribution%cdf(ec_start - 1)
         contributions(i, 1) = tail_contribution(k_index, denominator, full_weight, &
            distorted_cdf, portfolio%discretized_loss(i), portfolio%discretized_pd(i))

         k_index = var_start - portfolio%loss_multiple(i)
         denominator = 1.0_dp
         if (var_start > 1) denominator = 1.0_dp - distribution%cdf(var_start - 1)
         contributions(i, 2) = tail_contribution(k_index, denominator, full_weight, &
            distorted_cdf, portfolio%discretized_loss(i), portfolio%discretized_pd(i))

         k_index = es_start - portfolio%loss_multiple(i)
         denominator = 1.0_dp
         if (es_start > 1) denominator = 1.0_dp - distribution%cdf(es_start - 1)
         contributions(i, 3) = tail_contribution(k_index, denominator, full_weight, &
            distorted_cdf, portfolio%discretized_loss(i), portfolio%discretized_pd(i))
      end do
   end subroutine analytical_risk_contributions

   function tail_contribution(k_index, denominator, weight, cdf, loss, pd) result(value)
      integer, intent(in) :: k_index
      real(dp), intent(in) :: denominator
      real(dp), intent(in) :: weight(:)
      real(dp), intent(in) :: cdf(:,:)
      real(dp), intent(in) :: loss, pd
      real(dp) :: value
      real(dp) :: weighted_tail

      if (denominator <= tiny(1.0_dp)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (k_index <= 1) then
         value = loss * pd / denominator
      else if (k_index - 1 > size(cdf, 2)) then
         value = 0.0_dp
      else
         weighted_tail = sum(weight * (1.0_dp - cdf(:, k_index - 1)))
         value = loss * pd * weighted_tail / denominator
      end if
   end function tail_contribution

   subroutine simulation_risk_contributions(portfolio, distribution, alpha, contributions, &
                                             status, message)
      type(credit_portfolio), intent(in) :: portfolio
      type(loss_distribution), intent(in) :: distribution
      real(dp), intent(in) :: alpha
      real(dp), allocatable, intent(out) :: contributions(:,:)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message

      real(dp) :: var_value, tau, threshold
      integer :: local_status

      status = 0
      message = ''
      allocate(contributions(portfolio%n_counterparties, 3))
      contributions = ieee_value(0.0_dp, ieee_quiet_nan)
      if (distribution%model_kind /= model_simulation) then
         status = 1
         message = 'simulation contributions require a simulated distribution'
         return
      end if
      if (.not. distribution%contributions_available) then
         status = 2
         message = 'stored counterparty loss scenarios are unavailable'
         return
      end if

      var_value = value_at_risk(distribution, alpha, local_status)
      if (local_status /= 0) then
         status = 3
         message = 'requested alpha is unavailable'
         return
      end if

      tau = search_tau_for_target_es(distribution, portfolio, &
                                     var_value - portfolio%analytical_expected_loss, local_status)
      if (local_status == 0) then
         threshold = value_at_risk(distribution, tau, local_status)
         call average_stored_tail(distribution, threshold, contributions(:, 1), local_status)
      end if

      tau = search_tau_for_target_es(distribution, portfolio, var_value, local_status)
      if (local_status == 0) then
         threshold = value_at_risk(distribution, tau, local_status)
         call average_stored_tail(distribution, threshold, contributions(:, 2), local_status)
      end if

      call average_stored_tail(distribution, var_value, contributions(:, 3), local_status)
      if (local_status /= 0) then
         status = 4
         message = 'stored scenarios do not extend far enough down the loss distribution'
      end if
   end subroutine simulation_risk_contributions

   subroutine average_stored_tail(distribution, threshold, average_loss, status)
      type(loss_distribution), intent(in) :: distribution
      real(dp), intent(in) :: threshold
      real(dp), intent(out) :: average_loss(:)
      integer, intent(out) :: status
      real(dp) :: denominator, scenario_total
      integer :: j, scenario

      average_loss = 0.0_dp
      denominator = 0.0_dp
      status = 0
      do j = 1, size(distribution%stored_scenario)
         scenario_total = sum(distribution%stored_counterparty_losses(:, j))
         if (scenario_total >= threshold - 1.0e-10_dp) then
            scenario = distribution%stored_scenario(j)
            average_loss = average_loss + distribution%stored_counterparty_losses(:, j) * &
                           distribution%likelihood(scenario)
            denominator = denominator + distribution%likelihood(scenario)
         end if
      end do
      if (denominator <= tiny(1.0_dp)) then
         average_loss = ieee_value(0.0_dp, ieee_quiet_nan)
         status = 1
      else
         average_loss = average_loss / denominator
      end if
   end subroutine average_stored_tail

end module gcpm_risk
