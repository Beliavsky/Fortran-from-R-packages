! SPDX-License-Identifier: Artistic-2.0
module ldhmm_modeling
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ldhmm_distribution, only : ecld_create, ecld_pdf, ecld_mean, ecld_sd, &
      ecld_kurtosis
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : absolute_acf, finite_count, log_sum_exp, normal_quantile, &
      population_kurtosis, population_skewness, quiet_nan, remove_outliers, &
      sample_sd, simple_moving_average, vector_mean
   use ldhmm_parameters, only : ldhmm_validate
   use ldhmm_status, only : LDHMM_SUCCESS, LDHMM_INVALID_ARGUMENT, LDHMM_NUMERICAL_ERROR
   use ldhmm_types, only : ecld_type, ldhmm_model
   implicit none
   private

   public :: ldhmm_state_distribution, ldhmm_state_distributions, ldhmm_state_pdf
   public :: ldhmm_mllk, ldhmm_log_forward, ldhmm_log_backward
   public :: ldhmm_decode, ldhmm_viterbi, ldhmm_conditional_prob
   public :: ldhmm_forecast_prob, ldhmm_forecast_state, ldhmm_forecast_volatility
   public :: ldhmm_pseudo_residuals, ldhmm_ld_stats, ldhmm_calc_stats_from_obs
   public :: ldhmm_decode_stats_history, ldhmm_sma, ldhmm_abs_acf
   public :: ldhmm_drop_outliers

contains

   function ldhmm_state_distribution(model, state) result(distribution)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in) :: state
      type(ecld_type) :: distribution
      real(dp) :: mu, sigma, lambda

      if (state < 1 .or. state > model%m) then
         distribution = ecld_create(lambda=1.0_dp, sigma=1.0_dp, mu=0.0_dp)
         return
      end if
      mu = min(1.0e5_dp, max(-1.0e5_dp, model%param(state,1)))
      sigma = min(1.0e5_dp, max(1.0e-8_dp, model%param(state,2)))
      if (model%param_nbr == 3) then
         lambda = min(10.0_dp, max(1.0e-2_dp, model%param(state,3)))
      else
         lambda = 1.0_dp
      end if
      distribution = ecld_create(lambda=lambda, sigma=sigma, mu=mu)
   end function ldhmm_state_distribution

   function ldhmm_state_distributions(model, states) result(distributions)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in), optional :: states(:)
      type(ecld_type), allocatable :: distributions(:)
      integer :: i

      if (present(states)) then
         allocate(distributions(size(states)))
         do i = 1, size(states)
            distributions(i) = ldhmm_state_distribution(model, states(i))
         end do
      else
         allocate(distributions(model%m))
         do i = 1, model%m
            distributions(i) = ldhmm_state_distribution(model, i)
         end do
      end if
   end function ldhmm_state_distributions

   function ldhmm_state_pdf(model, states, x) result(pdf)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in) :: states(:)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: pdf(:, :)
      type(ecld_type) :: distribution
      integer :: i

      allocate(pdf(size(states),size(x)))
      do i = 1, size(states)
         distribution = ldhmm_state_distribution(model, states(i))
         pdf(i,:) = ecld_pdf(distribution, x)
      end do
   end function ldhmm_state_pdf

   real(dp) function ldhmm_mllk(model, x, status) result(value)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: log_alpha(:, :)
      integer :: local_status, nobs

      nobs = size(x)
      call ldhmm_log_forward(model, x, log_alpha, local_status)
      if (local_status /= LDHMM_SUCCESS .or. nobs == 0) then
         value = huge(1.0_dp)
         if (present(status)) status = local_status
         return
      end if
      value = -log_sum_exp(log_alpha(:,nobs))
      if (.not. ieee_is_finite(value)) then
         value = huge(1.0_dp)
         local_status = LDHMM_NUMERICAL_ERROR
      end if
      if (present(status)) status = local_status
   end function ldhmm_mllk

   subroutine ldhmm_log_forward(model, x, log_alpha, status)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: log_alpha(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: log_pdf(:, :), work(:)
      integer :: i, j, t, local_status

      local_status = ldhmm_validate(model)
      if (local_status /= LDHMM_SUCCESS .or. size(x) == 0) then
         allocate(log_alpha(0,0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      call emission_log_matrix(model, x, log_pdf)
      allocate(log_alpha(model%m,size(x)), work(model%m))
      do j = 1, model%m
         log_alpha(j,1) = safe_log(model%delta(j)) + log_pdf(j,1)
      end do
      do t = 2, size(x)
         do j = 1, model%m
            do i = 1, model%m
               work(i) = log_alpha(i,t-1) + safe_log(model%gamma(i,j))
            end do
            log_alpha(j,t) = log_pdf(j,t) + log_sum_exp(work)
         end do
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end subroutine ldhmm_log_forward

   subroutine ldhmm_log_backward(model, x, log_beta, status)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: log_beta(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: log_pdf(:, :), work(:)
      integer :: i, j, t, local_status

      local_status = ldhmm_validate(model)
      if (local_status /= LDHMM_SUCCESS .or. size(x) == 0) then
         allocate(log_beta(0,0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      call emission_log_matrix(model, x, log_pdf)
      allocate(log_beta(model%m,size(x)), work(model%m))
      log_beta(:,size(x)) = 0.0_dp
      do t = size(x)-1, 1, -1
         do i = 1, model%m
            do j = 1, model%m
               work(j) = safe_log(model%gamma(i,j)) + log_pdf(j,t+1) + &
                  log_beta(j,t+1)
            end do
            log_beta(i,t) = log_sum_exp(work)
         end do
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end subroutine ldhmm_log_backward

   function ldhmm_decode(input_model, x, do_global, do_stats, status) result(model)
      type(ldhmm_model), intent(in) :: input_model
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: do_global, do_stats
      integer, intent(out), optional :: status
      type(ldhmm_model) :: model
      real(dp), allocatable :: log_alpha(:, :), log_beta(:, :)
      real(dp) :: log_likelihood, column_sum
      logical :: global_decode, calculate_stats
      integer :: i, local_status, nobs

      nobs = size(x)
      global_decode = .true.
      calculate_stats = .true.
      if (present(do_global)) global_decode = do_global
      if (present(do_stats)) calculate_stats = do_stats
      model = input_model
      call clear_decoding(model)
      call ldhmm_log_forward(model, x, log_alpha, local_status)
      if (local_status /= LDHMM_SUCCESS) then
         if (present(status)) status = local_status
         return
      end if
      call ldhmm_log_backward(model, x, log_beta, local_status)
      log_likelihood = log_sum_exp(log_alpha(:,nobs))
      allocate(model%observations(size(x)), model%states_prob(model%m,size(x)))
      allocate(model%states_local(size(x)))
      model%observations = x
      do i = 1, size(x)
         model%states_prob(:,i) = exp(log_alpha(:,i)+log_beta(:,i)-log_likelihood)
         column_sum = sum(model%states_prob(:,i))
         if (column_sum > 0.0_dp) model%states_prob(:,i) = model%states_prob(:,i)/column_sum
         model%states_local(i) = maxloc(model%states_prob(:,i), dim=1)
      end do
      if (global_decode) model%states_global = ldhmm_viterbi(model, x, local_status)
      if (calculate_stats) then
         model%states_local_stats = ldhmm_calc_stats_from_obs(model, use_local=.true.)
         if (global_decode) then
            model%states_global_stats = ldhmm_calc_stats_from_obs(model, use_local=.false.)
         end if
      end if
      if (present(status)) status = local_status
   end function ldhmm_decode

   function ldhmm_viterbi(model, x, status) result(states)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: status
      integer, allocatable :: states(:)
      real(dp), allocatable :: log_pdf(:, :), score(:, :), candidates(:)
      integer, allocatable :: previous(:, :)
      integer :: i, j, t

      if (ldhmm_validate(model) /= LDHMM_SUCCESS .or. size(x) == 0) then
         allocate(states(0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      call emission_log_matrix(model, x, log_pdf)
      allocate(score(model%m,size(x)), previous(model%m,size(x)), candidates(model%m))
      allocate(states(size(x)))
      previous = 0
      do j = 1, model%m
         score(j,1) = safe_log(model%delta(j)) + log_pdf(j,1)
      end do
      do t = 2, size(x)
         do j = 1, model%m
            do i = 1, model%m
               candidates(i) = score(i,t-1) + safe_log(model%gamma(i,j))
            end do
            previous(j,t) = maxloc(candidates, dim=1)
            score(j,t) = candidates(previous(j,t)) + log_pdf(j,t)
         end do
      end do
      states(size(x)) = maxloc(score(:,size(x)), dim=1)
      do t = size(x)-1, 1, -1
         states(t) = previous(states(t+1),t+1)
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_viterbi

   function ldhmm_conditional_prob(model, x, xc, status) result(density)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:), xc(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: density(:, :), pdf_grid(:, :)
      real(dp), allocatable :: log_alpha(:, :), log_beta(:, :), log_previous(:)
      real(dp), allocatable :: log_weights(:), transition_work(:), weights(:)
      integer :: i, j, t, local_status

      call ldhmm_log_forward(model, x, log_alpha, local_status)
      if (local_status /= LDHMM_SUCCESS) then
         allocate(density(0,0))
         if (present(status)) status = local_status
         return
      end if
      call ldhmm_log_backward(model, x, log_beta, local_status)
      pdf_grid = ldhmm_state_pdf(model, [(i,i=1,model%m)], xc)
      allocate(density(size(xc),size(x)), log_previous(model%m))
      allocate(log_weights(model%m), transition_work(model%m), weights(model%m))
      do t = 1, size(x)
         if (t == 1) then
            do i = 1, model%m
               log_previous(i) = safe_log(model%delta(i))
            end do
         else
            log_previous = log_alpha(:,t-1)
         end if
         do j = 1, model%m
            do i = 1, model%m
               transition_work(i) = log_previous(i) + safe_log(model%gamma(i,j))
            end do
            log_weights(j) = log_sum_exp(transition_work) + log_beta(j,t)
         end do
         weights = exp(log_weights-log_sum_exp(log_weights))
         density(:,t) = matmul(transpose(pdf_grid), weights)
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_conditional_prob

   function ldhmm_forecast_state(model, x, horizon, status) result(probabilities)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: horizon
      integer, intent(out), optional :: status
      real(dp), allocatable :: probabilities(:, :), log_alpha(:, :), phi(:)
      integer :: h, local_status, nobs

      nobs = size(x)
      call ldhmm_log_forward(model, x, log_alpha, local_status)
      if (local_status /= LDHMM_SUCCESS .or. horizon < 1) then
         allocate(probabilities(0,0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      allocate(probabilities(model%m,horizon), phi(model%m))
      phi = exp(log_alpha(:,nobs)-log_sum_exp(log_alpha(:,nobs)))
      do h = 1, horizon
         phi = matmul(phi, model%gamma)
         probabilities(:,h) = phi / sum(phi)
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_forecast_state

   function ldhmm_forecast_prob(model, x, xf, horizon, status) result(density)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:), xf(:)
      integer, intent(in) :: horizon
      integer, intent(out), optional :: status
      real(dp), allocatable :: density(:, :), state_forecast(:, :), pdf_future(:, :)
      integer :: h, i, local_status

      state_forecast = ldhmm_forecast_state(model, x, horizon, local_status)
      if (local_status /= LDHMM_SUCCESS) then
         allocate(density(0,0))
         if (present(status)) status = local_status
         return
      end if
      pdf_future = ldhmm_state_pdf(model, [(i,i=1,model%m)], xf)
      allocate(density(horizon,size(xf)))
      do h = 1, horizon
         density(h,:) = matmul(state_forecast(:,h), pdf_future)
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_forecast_prob

   function ldhmm_pseudo_residuals(model, x, grid_length, status) result(residuals)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: grid_length
      integer, intent(out), optional :: status
      real(dp), allocatable :: residuals(:), grid(:), conditional(:, :)
      real(dp) :: lower, upper, dx, cumulative, probability, data_range
      integer :: i, j, ngrid, local_status

      ngrid = 1000
      if (present(grid_length)) ngrid = max(20, grid_length)
      if (size(x) == 0) then
         allocate(residuals(0))
         if (present(status)) status = LDHMM_INVALID_ARGUMENT
         return
      end if
      data_range = max(maxval(x)-minval(x), 1.0e-6_dp)
      lower = minval(x) - 0.25_dp*data_range
      upper = maxval(x) + 0.25_dp*data_range
      allocate(grid(ngrid), residuals(size(x)))
      dx = (upper-lower)/real(ngrid-1,dp)
      do i = 1, ngrid
         grid(i) = lower + real(i-1,dp)*dx
      end do
      conditional = ldhmm_conditional_prob(model, x, grid, local_status)
      if (local_status /= LDHMM_SUCCESS) then
         residuals = quiet_nan()
         if (present(status)) status = local_status
         return
      end if
      do i = 1, size(x)
         j = count(grid <= x(i))
         j = max(1, min(ngrid-1, j))
         cumulative = 0.0_dp
         if (j > 1) cumulative = sum(conditional(1:j-1,i))*dx
         probability = cumulative + 0.5_dp*conditional(j,i)*dx
         probability = min(1.0_dp-1.0e-12_dp, max(1.0e-12_dp, probability))
         residuals(i) = normal_quantile(probability)
      end do
      if (present(status)) status = LDHMM_SUCCESS
   end function ldhmm_pseudo_residuals

   function ldhmm_ld_stats(model, annualize, days_per_year) result(statistics)
      type(ldhmm_model), intent(in) :: model
      logical, intent(in), optional :: annualize
      integer, intent(in), optional :: days_per_year
      real(dp), allocatable :: statistics(:, :)
      type(ecld_type) :: distribution
      logical :: use_annual
      integer :: i, days

      use_annual = .false.
      days = 252
      if (present(annualize)) use_annual = annualize
      if (present(days_per_year)) days = days_per_year
      allocate(statistics(model%m,4))
      do i = 1, model%m
         distribution = ldhmm_state_distribution(model, i)
         statistics(i,1) = ecld_mean(distribution)
         statistics(i,2) = ecld_sd(distribution)
         statistics(i,3) = ecld_kurtosis(distribution)
         statistics(i,4) = statistics(i,1)/statistics(i,2)**2
      end do
      if (use_annual) then
         statistics(:,1) = statistics(:,1)*real(days,dp)
         statistics(:,2) = statistics(:,2)*sqrt(real(days,dp))*100.0_dp
         statistics(:,4) = statistics(:,1)/(statistics(:,2)/100.0_dp)**2
      end if
   end function ldhmm_ld_stats

   function ldhmm_calc_stats_from_obs(model, drop, use_local) result(statistics)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in), optional :: drop
      logical, intent(in), optional :: use_local
      real(dp), allocatable :: statistics(:, :), selected(:), cleaned(:)
      logical :: local_states
      integer :: i, j, n, remove_n

      local_states = .true.
      remove_n = 0
      if (present(use_local)) local_states = use_local
      if (present(drop)) remove_n = max(0, drop)
      allocate(statistics(model%m,6))
      statistics = quiet_nan()
      if (.not. allocated(model%observations)) return
      do i = 1, model%m
         if (local_states) then
            if (.not. allocated(model%states_local)) cycle
            n = count(model%states_local == i)
            allocate(selected(n))
            j = 0
            call gather_state_values(model%observations, model%states_local, i, selected, j)
         else
            if (.not. allocated(model%states_global)) cycle
            n = count(model%states_global == i)
            allocate(selected(n))
            j = 0
            call gather_state_values(model%observations, model%states_global, i, selected, j)
         end if
         cleaned = remove_outliers(selected, remove_n)
         statistics(i,1) = vector_mean(cleaned)
         statistics(i,2) = sample_sd(cleaned)
         statistics(i,3) = population_kurtosis(cleaned)
         statistics(i,4) = population_skewness(cleaned)
         statistics(i,5) = real(size(cleaned),dp)
         if (ieee_is_finite(statistics(i,2)) .and. statistics(i,2) > 0.0_dp) then
            statistics(i,6) = statistics(i,1)/statistics(i,2)**2
         end if
         deallocate(selected)
         if (allocated(cleaned)) deallocate(cleaned)
      end do
   end function ldhmm_calc_stats_from_obs

   function ldhmm_decode_stats_history(model, moving_average_order, annualize, &
      days_per_year) result(statistics)
      type(ldhmm_model), intent(in) :: model
      integer, intent(in), optional :: moving_average_order
      logical, intent(in), optional :: annualize
      integer, intent(in), optional :: days_per_year
      real(dp), allocatable :: statistics(:, :), state_stats(:, :), raw(:, :)
      real(dp), allocatable :: smooth(:)
      real(dp) :: mixture_mean, mixture_variance
      integer :: i, order, days
      logical :: use_annual

      order = 0
      days = 252
      use_annual = .false.
      if (present(moving_average_order)) order = max(0, moving_average_order)
      if (present(days_per_year)) days = days_per_year
      if (present(annualize)) use_annual = annualize
      if (.not. allocated(model%states_prob)) then
         allocate(statistics(0,0))
         return
      end if
      state_stats = ldhmm_ld_stats(model)
      allocate(raw(size(model%states_prob,2),3))
      do i = 1, size(raw,1)
         mixture_mean = dot_product(model%states_prob(:,i), state_stats(:,1))
         mixture_variance = dot_product(model%states_prob(:,i), &
            (state_stats(:,1)-mixture_mean)**2+state_stats(:,2)**2)
         raw(i,1) = mixture_mean
         raw(i,2) = sqrt(max(0.0_dp, mixture_variance))
         raw(i,3) = dot_product(model%states_prob(:,i), state_stats(:,3))
      end do
      allocate(statistics(size(raw,1),3))
      do i = 1, 3
         smooth = simple_moving_average(raw(:,i), order)
         statistics(:,i) = smooth
      end do
      if (use_annual) then
         statistics(:,1) = statistics(:,1)*real(days,dp)
         statistics(:,2) = statistics(:,2)*sqrt(real(days,dp))*100.0_dp
      end if
   end function ldhmm_decode_stats_history

   function ldhmm_forecast_volatility(model, x, xf, moving_average_order, &
      days_per_year, status) result(forecast)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:), xf(:)
      integer, intent(in), optional :: moving_average_order, days_per_year
      integer, intent(out), optional :: status
      real(dp), allocatable :: forecast(:, :), extended(:), history(:, :)
      type(ldhmm_model) :: decoded
      integer :: i, order, days, local_status

      order = 0
      days = 252
      if (present(moving_average_order)) order = max(0, moving_average_order)
      if (present(days_per_year)) days = days_per_year
      local_status = LDHMM_SUCCESS
      allocate(forecast(2,size(xf)), extended(size(x)+1))
      forecast = quiet_nan()
      if (size(x) > 0) extended(1:size(x)) = x
      do i = 1, size(xf)
         extended(size(extended)) = xf(i)
         decoded = ldhmm_decode(model, extended, do_global=.false., &
            do_stats=.false., status=local_status)
         forecast(1,i) = xf(i)
         if (local_status /= LDHMM_SUCCESS) exit
         history = ldhmm_decode_stats_history(decoded, order)
         if (size(history,1) == 0) then
            local_status = LDHMM_NUMERICAL_ERROR
            exit
         end if
         forecast(2,i) = history(size(history,1),2)*sqrt(real(days,dp))*100.0_dp
      end do
      if (present(status)) status = local_status
   end function ldhmm_forecast_volatility

   function ldhmm_sma(x, order, na_backfill) result(moving_average)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: order
      logical, intent(in), optional :: na_backfill
      real(dp), allocatable :: moving_average(:)
      if (present(na_backfill)) then
         moving_average = simple_moving_average(x, order, na_backfill)
      else
         moving_average = simple_moving_average(x, order)
      end if
   end function ldhmm_sma

   function ldhmm_abs_acf(x, lag_max, drop) result(acf)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lag_max
      integer, intent(in), optional :: drop
      real(dp), allocatable :: acf(:)
      if (present(drop)) then
         acf = absolute_acf(x, lag_max, drop)
      else
         acf = absolute_acf(x, lag_max)
      end if
   end function ldhmm_abs_acf

   function ldhmm_drop_outliers(x, drop) result(values)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: drop
      real(dp), allocatable :: values(:)
      integer :: remove_n
      remove_n = 1
      if (present(drop)) remove_n = max(0, drop)
      values = remove_outliers(x, remove_n)
   end function ldhmm_drop_outliers

   subroutine emission_log_matrix(model, x, log_pdf)
      type(ldhmm_model), intent(in) :: model
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: log_pdf(:, :)
      type(ecld_type) :: distribution
      real(dp) :: density
      integer :: i, t

      allocate(log_pdf(model%m,size(x)))
      do i = 1, model%m
         distribution = ldhmm_state_distribution(model, i)
         do t = 1, size(x)
            if (ieee_is_finite(x(t))) then
               density = ecld_pdf(distribution, x(t))
               log_pdf(i,t) = safe_log(density)
            else
               log_pdf(i,t) = 0.0_dp
            end if
         end do
      end do
   end subroutine emission_log_matrix

   elemental real(dp) function safe_log(x) result(value)
      real(dp), intent(in) :: x
      value = log(max(x, 1.0e-300_dp))
   end function safe_log

   subroutine clear_decoding(model)
      type(ldhmm_model), intent(inout) :: model
      if (allocated(model%observations)) deallocate(model%observations)
      if (allocated(model%states_prob)) deallocate(model%states_prob)
      if (allocated(model%states_local)) deallocate(model%states_local)
      if (allocated(model%states_global)) deallocate(model%states_global)
      if (allocated(model%states_local_stats)) deallocate(model%states_local_stats)
      if (allocated(model%states_global_stats)) deallocate(model%states_global_stats)
   end subroutine clear_decoding

   subroutine gather_state_values(observations, states, target, selected, n_selected)
      real(dp), intent(in) :: observations(:)
      integer, intent(in) :: states(:), target
      real(dp), intent(out) :: selected(:)
      integer, intent(out) :: n_selected
      integer :: i
      n_selected = 0
      do i = 1, min(size(observations),size(states))
         if (states(i) /= target) cycle
         n_selected = n_selected + 1
         selected(n_selected) = observations(i)
      end do
   end subroutine gather_state_values

end module ldhmm_modeling
