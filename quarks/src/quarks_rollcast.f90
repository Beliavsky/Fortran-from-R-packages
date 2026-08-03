module quarks_rollcast
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use quarks_kinds, only : dp
   use quarks_rng, only : rng_state
   use quarks_risk, only : hs, vwhs, fhs
   use quarks_smoothing, only : smooth_scale
   use quarks_stats, only : nan_value
   use quarks_types, only : risk_result, rollcast_result, quarks_ok, &
      quarks_invalid_input, method_plain, method_age, method_vwhs, method_fhs, &
      volatility_ewma, smooth_none
   implicit none
   private

   public :: rollcast

contains

   function rollcast(x, p, method, volatility_model, lambda, nout, nwin, &
      nboot, smoothing, smoothing_bandwidth, rng, garch_max_iterations) &
      result(result)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: p, lambda, smoothing_bandwidth
      integer, intent(in), optional :: method, volatility_model, nout, nwin
      integer, intent(in), optional :: nboot, smoothing, garch_max_iterations
      type(rng_state), intent(inout), optional :: rng
      type(rollcast_result) :: result
      type(risk_result) :: current
      real(dp), allocatable :: window(:), standardized(:)
      real(dp) :: probability, decay, scale
      integer :: selected_method, selected_model, out_count, window_size
      integer :: bootstrap_size, selected_smoothing, max_iterations
      integer :: n, nin, forecasts, i, first, last, smooth_status

      probability = 0.975_dp
      selected_method = method_plain
      selected_model = volatility_ewma
      out_count = 0
      window_size = 250
      bootstrap_size = 10000
      selected_smoothing = smooth_none
      max_iterations = 1200
      if (present(p)) probability = p
      if (present(method)) selected_method = method
      if (present(volatility_model)) selected_model = volatility_model
      if (present(nout)) out_count = nout
      if (present(nwin)) window_size = nwin
      if (present(nboot)) bootstrap_size = nboot
      if (present(smoothing)) selected_smoothing = smoothing
      if (present(garch_max_iterations)) max_iterations = garch_max_iterations
      if (present(lambda)) then
         decay = lambda
      else if (selected_method == method_age) then
         decay = 0.98_dp
      else
         decay = 0.94_dp
      end if
      result%p = probability
      result%method = selected_method
      result%volatility_model = selected_model
      result%smoothing = selected_smoothing
      result%nout = out_count
      result%nwin = window_size
      result%nboot = merge(bootstrap_size, 0, selected_method == method_fhs)
      n = size(x)
      if (n <= 1 .or. probability <= 0.0_dp .or. probability >= 1.0_dp .or. &
          out_count < 0 .or. window_size <= 1 .or. &
          window_size + out_count > n .or. any(.not. ieee_is_finite(x)) .or. &
          selected_method < method_plain .or. selected_method > method_fhs) then
         allocate(result%var(0), result%es(0), result%xout(0))
         result%status = quarks_invalid_input
         result%message = 'invalid rolling-forecast inputs'
         return
      end if
      if (selected_method == method_fhs .and. bootstrap_size <= 0) then
         allocate(result%var(0), result%es(0), result%xout(0))
         result%status = quarks_invalid_input
         result%message = 'nboot must be positive for filtered historical simulation'
         return
      end if
      nin = n - out_count
      forecasts = max(out_count, 1)
      allocate(result%var(forecasts), result%es(forecasts), result%xout(out_count))
      if (out_count > 0) result%xout = x(nin + 1:n)
      result%var = nan_value()
      result%es = nan_value()
      do i = 1, forecasts
         if (out_count == 0) then
            last = n
         else
            last = nin + i - 1
         end if
         first = last - window_size + 1
         allocate(window(window_size))
         window = x(first:last)
         if (selected_smoothing == smooth_none) then
            allocate(standardized(window_size))
            standardized = window
            scale = 1.0_dp
         else if (present(smoothing_bandwidth)) then
            call smooth_scale(window, selected_smoothing, standardized, scale, &
               bandwidth=smoothing_bandwidth, status=smooth_status)
         else
            call smooth_scale(window, selected_smoothing, standardized, scale, &
               status=smooth_status)
         end if
         select case (selected_method)
         case (method_plain)
            current = hs(standardized, probability, method_plain, decay)
         case (method_age)
            current = hs(standardized, probability, method_age, decay)
         case (method_vwhs)
            current = vwhs(standardized, probability, selected_model, decay, &
               max_iterations)
         case (method_fhs)
            if (present(rng)) then
               current = fhs(standardized, probability, selected_model, decay, &
                  bootstrap_size, rng, max_iterations)
            else
               current = fhs(standardized, probability, selected_model, decay, &
                  bootstrap_size, garch_max_iterations=max_iterations)
            end if
         end select
         result%var(i) = current%var * scale
         result%es(i) = current%es * scale
         if (current%status /= quarks_ok .and. result%status == quarks_ok) then
            result%status = current%status
            result%message = current%message
         end if
         deallocate(window, standardized)
      end do
   end function rollcast

end module quarks_rollcast
