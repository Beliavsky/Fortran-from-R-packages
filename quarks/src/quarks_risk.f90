module quarks_risk
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use quarks_kinds, only : dp
   use quarks_rng, only : rng_state, random_index
   use quarks_stats, only : sample_variance, quantile_type7, sort_real_with_index, &
      nan_value
   use quarks_types, only : risk_result, quarks_ok, quarks_invalid_input, &
      quarks_empty_tail, quarks_fit_failed, volatility_ewma, volatility_garch
   use rugarch, only : garch_fit_result, fit_garch11, dist_norm, &
      forecast_volatility
   implicit none
   private

   public :: ewma, hs, vwhs, fhs

contains

   function ewma(x, lambda, status) result(variance)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: lambda
      integer, intent(out), optional :: status
      real(dp), allocatable :: variance(:)
      real(dp) :: decay
      integer :: i, n
      n = size(x)
      allocate(variance(n))
      decay = 0.94_dp
      if (present(lambda)) decay = lambda
      if (present(status)) status = quarks_ok
      if (n <= 1 .or. decay < 0.0_dp .or. decay >= 1.0_dp .or. &
          any(.not. ieee_is_finite(x))) then
         variance = nan_value()
         if (present(status)) status = quarks_invalid_input
         return
      end if
      variance(1) = sample_variance(x)
      do i = 2, n
         variance(i) = decay * variance(i - 1) + (1.0_dp - decay) * x(i - 1)**2
      end do
   end function ewma

   function hs(x, p, method, lambda) result(result)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: p, lambda
      integer, intent(in), optional :: method
      type(risk_result) :: result
      real(dp) :: probability, decay, total_tail_weight
      real(dp), allocatable :: loss(:), sorted_loss(:), weights(:), tail_weights(:)
      integer, allocatable :: indices(:)
      integer :: selected_method, i, n, tail_count

      probability = 0.975_dp
      decay = 0.98_dp
      selected_method = 1
      if (present(p)) probability = p
      if (present(lambda)) decay = lambda
      if (present(method)) selected_method = method
      result%p = probability
      if (size(x) <= 1 .or. probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
          any(.not. ieee_is_finite(x))) then
         result%status = quarks_invalid_input
         result%message = 'invalid return series or probability'
         result%var = nan_value()
         result%es = nan_value()
         return
      end if
      n = size(x)
      allocate(loss(n))
      loss = -x
      call sort_real_with_index(loss, sorted_loss, indices)
      select case (selected_method)
      case (1)
         result%var = quantile_type7(loss, probability)
         tail_count = count(loss > result%var)
         if (tail_count == 0) then
            result%es = nan_value()
            result%status = quarks_empty_tail
            result%message = 'no losses strictly above VaR'
         else
            result%es = sum(loss, mask=loss > result%var) / real(tail_count, dp)
         end if
      case (2)
         if (decay < 0.0_dp .or. decay >= 1.0_dp) then
            result%status = quarks_invalid_input
            result%message = 'invalid age-decay factor'
            result%var = nan_value()
            result%es = nan_value()
            return
         end if
         allocate(weights(n))
         if (decay <= 0.0_dp) then
            weights = 0.0_dp
            weights(n) = 1.0_dp
         else
            do i = 1, n
               weights(i) = decay**real(n - i, dp)
            end do
            weights = weights / sum(weights)
         end if
         result%var = upstream_age_quantile(sorted_loss, indices, weights, probability)
         tail_count = count(sorted_loss > result%var)
         if (tail_count == 0) then
            result%es = nan_value()
            result%status = quarks_empty_tail
            result%message = 'no losses strictly above VaR'
         else
            allocate(tail_weights(tail_count))
            tail_weights = pack(weights(indices), sorted_loss > result%var)
            total_tail_weight = sum(tail_weights)
            if (total_tail_weight <= 0.0_dp) then
               result%es = nan_value()
               result%status = quarks_empty_tail
               result%message = 'zero weight in the ES tail'
            else
               tail_weights = tail_weights / total_tail_weight
               result%es = sum(pack(sorted_loss, sorted_loss > result%var) * &
                  tail_weights)
            end if
         end if
      case default
         result%status = quarks_invalid_input
         result%message = 'method must be 1 (plain) or 2 (age)'
         result%var = nan_value()
         result%es = nan_value()
      end select
   end function hs

   function upstream_age_quantile(sorted_loss, indices, weights, p) result(value)
      real(dp), intent(in) :: sorted_loss(:), weights(:), p
      integer, intent(in) :: indices(:)
      real(dp) :: value, cumulative, previous, fraction
      integer :: i
      cumulative = weights(indices(1))
      if (p <= cumulative) then
         value = sorted_loss(1)
         return
      end if
      previous = cumulative
      do i = 2, size(sorted_loss)
         cumulative = cumulative + weights(indices(i))
         if (cumulative > p) then
            fraction = (p - previous) / max(cumulative - previous, tiny(1.0_dp))
            value = sorted_loss(i - 1) + fraction * &
               (sorted_loss(i) - sorted_loss(i - 1))
            return
         end if
         previous = cumulative
      end do
      value = sorted_loss(size(sorted_loss))
   end function upstream_age_quantile

   function vwhs(x, p, model, lambda, garch_max_iterations) result(result)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: p, lambda
      integer, intent(in), optional :: model, garch_max_iterations
      type(risk_result) :: result
      type(garch_fit_result) :: fit
      real(dp), allocatable :: conditional_variance(:), sigma(:), standardized(:), loss(:)
      real(dp) :: probability, decay
      integer :: selected_model, max_iterations, n, tail_count, ewma_status

      probability = 0.975_dp
      decay = 0.94_dp
      selected_model = volatility_ewma
      max_iterations = 1200
      if (present(p)) probability = p
      if (present(lambda)) decay = lambda
      if (present(model)) selected_model = model
      if (present(garch_max_iterations)) max_iterations = garch_max_iterations
      result%p = probability
      result%volatility_model = selected_model
      n = size(x)
      if (n <= 1 .or. probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
          any(.not. ieee_is_finite(x))) then
         result%status = quarks_invalid_input
         result%message = 'invalid return series or probability'
         result%var = nan_value()
         result%es = nan_value()
         return
      end if
      allocate(sigma(n), standardized(n), loss(n), conditional_variance(n))
      if (selected_model == volatility_garch) then
         fit = fit_garch11(x, cond_dist=dist_norm, fit_mean=.true., &
            max_iterations=max_iterations)
         if (fit%status == 0 .and. allocated(fit%sigma) .and. &
             all(ieee_is_finite(fit%sigma)) .and. all(fit%sigma > 0.0_dp)) then
            sigma = fit%sigma
         else
            conditional_variance = ewma(x, decay, ewma_status)
            sigma = sqrt(conditional_variance)
            result%used_fallback = .true.
            result%status = quarks_fit_failed
            result%message = 'GARCH fit failed; EWMA fallback used'
         end if
      else if (selected_model == volatility_ewma) then
         conditional_variance = ewma(x, decay, ewma_status)
         if (ewma_status /= quarks_ok) then
            result%status = quarks_invalid_input
            result%message = 'invalid EWMA inputs'
            result%var = nan_value()
            result%es = nan_value()
            return
         end if
         sigma = sqrt(conditional_variance)
      else
         result%status = quarks_invalid_input
         result%message = 'unknown volatility model'
         result%var = nan_value()
         result%es = nan_value()
         return
      end if
      standardized = x / max(sigma, sqrt(tiny(1.0_dp)))
      loss = -(standardized * sigma(n))
      result%var = quantile_type7(loss, probability)
      tail_count = count(loss > result%var)
      if (tail_count == 0) then
         result%es = nan_value()
         result%status = quarks_empty_tail
         result%message = 'no losses strictly above VaR'
      else
         result%es = sum(loss, mask=loss > result%var) / real(tail_count, dp)
      end if
   end function vwhs

   function fhs(x, p, model, lambda, nboot, rng, garch_max_iterations) result(result)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: p, lambda
      integer, intent(in), optional :: model, nboot, garch_max_iterations
      type(rng_state), intent(inout), optional :: rng
      type(risk_result) :: result
      type(garch_fit_result) :: fit
      real(dp), allocatable :: conditional_variance(:), sigma(:), standardized(:)
      real(dp), allocatable :: bootstrap_loss(:), forecast(:)
      real(dp) :: probability, decay, one_ahead_sigma, u
      integer :: selected_model, bootstrap_size, max_iterations
      integer :: n, i, index, tail_count, ewma_status

      probability = 0.975_dp
      decay = 0.94_dp
      selected_model = volatility_ewma
      bootstrap_size = 10000
      max_iterations = 1200
      if (present(p)) probability = p
      if (present(lambda)) decay = lambda
      if (present(model)) selected_model = model
      if (present(nboot)) bootstrap_size = nboot
      if (present(garch_max_iterations)) max_iterations = garch_max_iterations
      result%p = probability
      result%volatility_model = selected_model
      n = size(x)
      if (n <= 1 .or. bootstrap_size <= 0 .or. probability <= 0.0_dp .or. &
          probability >= 1.0_dp .or. any(.not. ieee_is_finite(x))) then
         result%status = quarks_invalid_input
         result%message = 'invalid filtered-historical-simulation inputs'
         result%var = nan_value()
         result%es = nan_value()
         return
      end if
      allocate(sigma(n), standardized(n), bootstrap_loss(bootstrap_size), conditional_variance(n))
      if (selected_model == volatility_garch) then
         fit = fit_garch11(x, cond_dist=dist_norm, fit_mean=.true., &
            max_iterations=max_iterations)
         if (fit%status == 0 .and. allocated(fit%sigma) .and. &
             allocated(fit%residuals)) then
            sigma = fit%sigma
            allocate(forecast(1))
            call forecast_volatility(fit%spec, fit%residuals, fit%sigma, 1, forecast)
            one_ahead_sigma = forecast(1)
         else
            conditional_variance = ewma(x, decay, ewma_status)
            sigma = sqrt(conditional_variance)
            one_ahead_sigma = sqrt(decay * conditional_variance(n) + &
               (1.0_dp - decay) * x(n)**2)
            result%used_fallback = .true.
            result%status = quarks_fit_failed
            result%message = 'GARCH fit failed; EWMA fallback used'
         end if
      else if (selected_model == volatility_ewma) then
         conditional_variance = ewma(x, decay, ewma_status)
         if (ewma_status /= quarks_ok) then
            result%status = quarks_invalid_input
            result%message = 'invalid EWMA inputs'
            result%var = nan_value()
            result%es = nan_value()
            return
         end if
         sigma = sqrt(conditional_variance)
         one_ahead_sigma = sqrt(decay * conditional_variance(n) + &
            (1.0_dp - decay) * x(n)**2)
      else
         result%status = quarks_invalid_input
         result%message = 'unknown volatility model'
         result%var = nan_value()
         result%es = nan_value()
         return
      end if
      standardized = x / max(sigma, sqrt(tiny(1.0_dp)))
      do i = 1, bootstrap_size
         if (present(rng)) then
            index = random_index(rng, n)
         else
            call random_number(u)
            index = min(n, 1 + int(u * real(n, dp)))
         end if
         bootstrap_loss(i) = -standardized(index) * one_ahead_sigma
      end do
      result%var = quantile_type7(bootstrap_loss, probability)
      tail_count = count(bootstrap_loss > result%var)
      if (tail_count == 0) then
         result%es = nan_value()
         result%status = quarks_empty_tail
         result%message = 'no bootstrap losses strictly above VaR'
      else
         result%es = sum(bootstrap_loss, mask=bootstrap_loss > result%var) / &
            real(tail_count, dp)
      end if
   end function fhs

end module quarks_risk
