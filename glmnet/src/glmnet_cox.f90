! SPDX-License-Identifier: GPL-2.0-only
module glmnet_cox
   use glmnet_kinds, only : dp, glmnet_eps, glmnet_huge
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument, &
      glmnet_nonfinite_input, glmnet_invalid_response, glmnet_max_iterations, &
      glmnet_no_events, glmnet_numerical_failure
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, &
      glmnet_family_cox, family_name
   use glmnet_utils, only : normalize_weights, prepare_design, make_lambda_sequence, &
      all_finite_vector, all_finite_matrix, soft_threshold, sort_indices_real
   implicit none
   private
   public :: fit_cox_path, cox_loss_gradient, coxnet_deviance, cox_gradient
   public :: concordance_index
contains
   subroutine fit_cox_path(x, start_time, stop_time, event, result, control, &
      weights_in, offset_in, strata_in, lambda_in, penalty_factor_in, lower_in, &
      upper_in, excluded_in, efron)
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:)
      integer, intent(in) :: event(:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      integer, intent(in), optional :: strata_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:), efron
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), offset(:), xw(:,:), xm(:), xs(:)
      real(dp), allocatable :: penalty(:), lower(:), upper(:), lambda(:), beta(:)
      real(dp), allocatable :: gradient(:), candidate(:), candidate_gradient(:)
      logical, allocatable :: usable(:), excluded(:)
      integer, allocatable :: strata(:)
      real(dp) :: null_loss, loss, candidate_loss, objective, candidate_objective
      real(dp) :: lambda_max, alpha_sequence, step, denominator, max_change
      integer :: n, p, j, l, iter, status, backtrack
      logical :: use_efron, converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      use_efron = .false.
      if (present(efron)) use_efron = efron
      n = size(x, 1)
      p = size(x, 2)
      call initialize_result(result, n, p)
      result%efron = use_efron
      if (n < 2 .or. p < 1 .or. size(start_time) /= n .or. &
          size(stop_time) /= n .or. size(event) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (.not. all_finite_matrix(x) .or. .not. all_finite_vector(start_time) .or. &
          .not. all_finite_vector(stop_time)) then
         result%status = glmnet_nonfinite_input
         return
      end if
      if (any(start_time < 0.0_dp) .or. any(stop_time <= start_time) .or. &
          any(event < 0) .or. any(event > 1)) then
         result%status = glmnet_invalid_response
         return
      end if
      if (sum(event) < 1) then
         result%status = glmnet_no_events
         return
      end if
      call normalize_weights(n, weights_in, weights, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(offset(n), strata(n))
      if (present(offset_in)) then
         if (size(offset_in) /= n .or. .not. all_finite_vector(offset_in)) then
            result%status = glmnet_invalid_argument
            return
         end if
         offset = offset_in
      else
         offset = 0.0_dp
      end if
      if (present(strata_in)) then
         if (size(strata_in) /= n .or. any(strata_in < 1)) then
            result%status = glmnet_invalid_argument
            return
         end if
         strata = strata_in
      else
         strata = 1
      end if
      call prepare_design(x, weights, ctl, xw, xm, xs, usable, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call initialize_constraints(p, penalty_factor_in, lower_in, upper_in, &
         excluded_in, xs, usable, penalty, lower, upper, excluded, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(beta(p), gradient(p), candidate(p), candidate_gradient(p))
      beta = 0.0_dp
      call cox_loss_gradient(xw, start_time, stop_time, event, strata, weights, &
         offset, beta, use_efron, null_loss, gradient, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      alpha_sequence = max(ctl%alpha, 1.0e-3_dp)
      lambda_max = 0.0_dp
      do j = 1, p
         if (.not. excluded(j)) lambda_max = max(lambda_max, abs(gradient(j)) / &
            max(alpha_sequence * penalty(j), glmnet_eps))
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_result(result, size(lambda), p)
      result%lambda = lambda
      result%nulldev = 2.0_dp * null_loss
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = .false.
      do l = 1, size(lambda)
         call cox_loss_gradient(xw, start_time, stop_time, event, strata, weights, &
            offset, beta, use_efron, loss, gradient, status)
         if (status /= glmnet_success) then
            result%status = status
            call truncate_result(result, max(l - 1, 1))
            return
         end if
         objective = loss + penalty_value(beta, penalty, lambda(l), ctl%alpha)
         step = initial_step_size(xw, weights, lambda(l), ctl%alpha, penalty)
         converged = .false.
         do iter = 1, ctl%max_iterations
            do j = 1, p
               if (excluded(j)) then
                  candidate(j) = 0.0_dp
               else
                  denominator = 1.0_dp + step * lambda(l) * &
                     (1.0_dp - ctl%alpha) * penalty(j)
                  candidate(j) = soft_threshold(beta(j) - step * gradient(j), &
                     step * lambda(l) * ctl%alpha * penalty(j)) / denominator
                  candidate(j) = min(max(candidate(j), lower(j)), upper(j))
               end if
            end do
            do backtrack = 1, 50
               call cox_loss_gradient(xw, start_time, stop_time, event, strata, &
                  weights, offset, candidate, use_efron, candidate_loss, &
                  candidate_gradient, status)
               if (status /= glmnet_success) exit
               candidate_objective = candidate_loss + penalty_value(candidate, &
                  penalty, lambda(l), ctl%alpha)
               if (candidate_objective <= objective + 1.0e-12_dp) exit
               step = 0.5_dp * step
               do j = 1, p
                  if (excluded(j)) then
                     candidate(j) = 0.0_dp
                  else
                     denominator = 1.0_dp + step * lambda(l) * &
                        (1.0_dp - ctl%alpha) * penalty(j)
                     candidate(j) = soft_threshold(beta(j) - step * gradient(j), &
                        step * lambda(l) * ctl%alpha * penalty(j)) / denominator
                     candidate(j) = min(max(candidate(j), lower(j)), upper(j))
                  end if
               end do
            end do
            if (status /= glmnet_success) exit
            max_change = maxval(abs(candidate - beta))
            beta = candidate
            loss = candidate_loss
            gradient = candidate_gradient
            objective = candidate_objective
            result%npasses = result%npasses + 1
            if (max_change <= ctl%threshold * (1.0_dp + maxval(abs(beta)))) then
               converged = .true.
               exit
            end if
            step = min(step * 1.1_dp, 1.0_dp)
         end do
         if (status /= glmnet_success) then
            result%status = glmnet_numerical_failure
            call truncate_result(result, l)
            return
         end if
         result%objective(l) = objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, &
            (null_loss - loss) / max(abs(null_loss), glmnet_eps)))
         result%iterations(l) = min(iter, ctl%max_iterations)
         result%converged(l) = converged
         result%beta(:, 1, l) = beta / xs
         result%intercept(1, l) = 0.0_dp
         result%df(l) = count(abs(result%beta(:, 1, l)) > 1.0e-12_dp)
         if (.not. converged .and. result%status == glmnet_success) &
            result%status = glmnet_max_iterations
      end do
   end subroutine fit_cox_path

   subroutine cox_loss_gradient(x, start_time, stop_time, event, strata, weights, &
      offset, beta, efron, loss, gradient, status)
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:), weights(:), offset(:), beta(:)
      integer, intent(in) :: event(:), strata(:)
      logical, intent(in) :: efron
      real(dp), intent(out) :: loss, gradient(size(beta))
      integer, intent(out) :: status
      real(dp), allocatable :: eta(:), exp_eta(:), event_times(:)
      real(dp), allocatable :: risk_x(:), event_exp_x(:)
      real(dp) :: risk_sum, event_exp_sum, event_weight, event_eta
      real(dp) :: denominator, fraction, tie_weight, total_event_weight
      integer :: s, ns, t_index, i, m, d
      integer, allocatable :: strata_values(:)
      logical :: in_risk, is_event

      status = glmnet_success
      loss = 0.0_dp
      gradient = 0.0_dp
      allocate(eta(size(start_time)), exp_eta(size(start_time)), risk_x(size(beta)), &
         event_exp_x(size(beta)))
      eta = matmul(x, beta) + offset
      exp_eta = exp(max(min(eta, 40.0_dp), -40.0_dp))
      call unique_integers(strata, strata_values)
      total_event_weight = sum(weights * real(event, dp))
      if (total_event_weight <= glmnet_eps) then
         status = glmnet_no_events
         return
      end if
      ns = size(strata_values)
      do s = 1, ns
         call event_times_for_stratum(stop_time, event, strata, strata_values(s), event_times)
         do t_index = 1, size(event_times)
            risk_sum = 0.0_dp
            risk_x = 0.0_dp
            event_exp_sum = 0.0_dp
            event_exp_x = 0.0_dp
            event_weight = 0.0_dp
            event_eta = 0.0_dp
            m = 0
            do i = 1, size(start_time)
               if (strata(i) /= strata_values(s)) cycle
               in_risk = start_time(i) < event_times(t_index) .and. &
                  stop_time(i) >= event_times(t_index)
               if (in_risk) then
                  risk_sum = risk_sum + weights(i) * exp_eta(i)
                  risk_x = risk_x + weights(i) * exp_eta(i) * x(i, :)
               end if
               is_event = event(i) == 1 .and. &
                  abs(stop_time(i) - event_times(t_index)) <= &
                  100.0_dp * glmnet_eps * max(1.0_dp, abs(event_times(t_index)))
               if (is_event) then
                  m = m + 1
                  event_weight = event_weight + weights(i)
                  event_eta = event_eta + weights(i) * eta(i)
                  event_exp_sum = event_exp_sum + weights(i) * exp_eta(i)
                  event_exp_x = event_exp_x + weights(i) * exp_eta(i) * x(i, :)
                  gradient = gradient - weights(i) * x(i, :)
               end if
            end do
            if (risk_sum <= glmnet_eps .or. m < 1) then
               status = glmnet_numerical_failure
               return
            end if
            loss = loss - event_eta
            if (efron .and. m > 1) then
               tie_weight = event_weight / real(m, dp)
               do d = 0, m - 1
                  fraction = real(d, dp) / real(m, dp)
                  denominator = risk_sum - fraction * event_exp_sum
                  if (denominator <= glmnet_eps) then
                     status = glmnet_numerical_failure
                     return
                  end if
                  loss = loss + tie_weight * log(denominator)
                  gradient = gradient + tie_weight * &
                     (risk_x - fraction * event_exp_x) / denominator
               end do
            else
               loss = loss + event_weight * log(risk_sum)
               gradient = gradient + event_weight * risk_x / risk_sum
            end if
         end do
      end do
      loss = loss / total_event_weight
      gradient = gradient / total_event_weight
   end subroutine cox_loss_gradient

   subroutine cox_gradient(x, start_time, stop_time, event, beta, gradient, status, &
      weights, offset, strata, efron)
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:), beta(:)
      integer, intent(in) :: event(:)
      real(dp), intent(out) :: gradient(size(beta))
      integer, intent(out) :: status
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: strata(:)
      logical, intent(in), optional :: efron
      real(dp), allocatable :: w(:), o(:)
      integer, allocatable :: s(:)
      real(dp) :: loss
      logical :: ef
      call local_survival_defaults(size(x, 1), weights, offset, strata, efron, w, o, s, ef, status)
      if (status /= glmnet_success) return
      call cox_loss_gradient(x, start_time, stop_time, event, s, w, o, beta, ef, loss, gradient, status)
   end subroutine cox_gradient

   subroutine coxnet_deviance(x, start_time, stop_time, event, beta, deviance, status, &
      weights, offset, strata, efron)
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:), beta(:)
      integer, intent(in) :: event(:)
      real(dp), intent(out) :: deviance
      integer, intent(out) :: status
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: strata(:)
      logical, intent(in), optional :: efron
      real(dp), allocatable :: w(:), o(:), gradient(:)
      integer, allocatable :: s(:)
      real(dp) :: loss
      logical :: ef
      allocate(gradient(size(beta)))
      call local_survival_defaults(size(x, 1), weights, offset, strata, efron, w, o, s, ef, status)
      if (status /= glmnet_success) return
      call cox_loss_gradient(x, start_time, stop_time, event, s, w, o, beta, ef, loss, gradient, status)
      deviance = 2.0_dp * loss
   end subroutine coxnet_deviance

   function concordance_index(time, event, score, weights) result(value)
      real(dp), intent(in) :: time(:), score(:)
      integer, intent(in) :: event(:)
      real(dp), intent(in), optional :: weights(:)
      real(dp) :: value
      real(dp), allocatable :: w(:)
      real(dp) :: comparable, concordant, pair_weight
      integer :: i, j, n
      n = size(time)
      allocate(w(n))
      if (present(weights)) then
         if (size(weights) == n) then
            w = weights
         else
            w = 1.0_dp
         end if
      else
         w = 1.0_dp
      end if
      comparable = 0.0_dp
      concordant = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            pair_weight = w(i) * w(j)
            if (event(i) == 1 .and. time(i) < time(j)) then
               comparable = comparable + pair_weight
               if (abs(score(i) - score(j)) <= 10.0_dp * glmnet_eps) then
                  concordant = concordant + 0.5_dp * pair_weight
               else if (score(i) > score(j)) then
                  concordant = concordant + pair_weight
               end if
            else if (event(j) == 1 .and. time(j) < time(i)) then
               comparable = comparable + pair_weight
               if (abs(score(i) - score(j)) <= 10.0_dp * glmnet_eps) then
                  concordant = concordant + 0.5_dp * pair_weight
               else if (score(j) > score(i)) then
                  concordant = concordant + pair_weight
               end if
            end if
         end do
      end do
      if (comparable <= glmnet_eps) then
         value = 0.5_dp
      else
         value = concordant / comparable
      end if
   end function concordance_index

   subroutine local_survival_defaults(n, weights_in, offset_in, strata_in, efron_in, &
      weights, offset, strata, efron, status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: weights_in(:), offset_in(:)
      integer, intent(in), optional :: strata_in(:)
      logical, intent(in), optional :: efron_in
      real(dp), allocatable, intent(out) :: weights(:), offset(:)
      integer, allocatable, intent(out) :: strata(:)
      logical, intent(out) :: efron
      integer, intent(out) :: status
      call normalize_weights(n, weights_in, weights, status)
      if (status /= glmnet_success) return
      allocate(offset(n), strata(n))
      offset = 0.0_dp
      strata = 1
      efron = .false.
      if (present(offset_in)) then
         if (size(offset_in) /= n) then
            status = glmnet_invalid_argument
            return
         end if
         offset = offset_in
      end if
      if (present(strata_in)) then
         if (size(strata_in) /= n) then
            status = glmnet_invalid_argument
            return
         end if
         strata = strata_in
      end if
      if (present(efron_in)) efron = efron_in
   end subroutine local_survival_defaults

   subroutine unique_integers(values, unique_values)
      integer, intent(in) :: values(:)
      integer, allocatable, intent(out) :: unique_values(:)
      integer, allocatable :: work(:)
      integer :: i, j, count_unique, key
      allocate(work(size(values)))
      work = values
      do i = 2, size(work)
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      count_unique = 1
      do i = 2, size(work)
         if (work(i) /= work(count_unique)) then
            count_unique = count_unique + 1
            work(count_unique) = work(i)
         end if
      end do
      allocate(unique_values(count_unique))
      unique_values = work(:count_unique)
   end subroutine unique_integers

   subroutine event_times_for_stratum(stop_time, event, strata, stratum, times)
      real(dp), intent(in) :: stop_time(:)
      integer, intent(in) :: event(:), strata(:), stratum
      real(dp), allocatable, intent(out) :: times(:)
      real(dp), allocatable :: candidates(:)
      integer, allocatable :: order(:)
      integer :: i, count_event, count_unique
      count_event = count(event == 1 .and. strata == stratum)
      allocate(candidates(count_event))
      count_event = 0
      do i = 1, size(stop_time)
         if (event(i) == 1 .and. strata(i) == stratum) then
            count_event = count_event + 1
            candidates(count_event) = stop_time(i)
         end if
      end do
      if (count_event == 0) then
         allocate(times(0))
         return
      end if
      call sort_indices_real(candidates, order)
      count_unique = 1
      do i = 2, count_event
         if (abs(candidates(order(i)) - candidates(order(count_unique))) > &
             100.0_dp * glmnet_eps * max(1.0_dp, abs(candidates(order(i))))) then
            count_unique = count_unique + 1
            order(count_unique) = order(i)
         end if
      end do
      allocate(times(count_unique))
      do i = 1, count_unique
         times(i) = candidates(order(i))
      end do
   end subroutine event_times_for_stratum

   pure function penalty_value(beta, penalty, lambda, alpha) result(value)
      real(dp), intent(in) :: beta(:), penalty(:), lambda, alpha
      real(dp) :: value
      value = lambda * sum(penalty * (alpha * abs(beta) + &
         0.5_dp * (1.0_dp - alpha) * beta ** 2))
   end function penalty_value

   pure function initial_step_size(x, weights, lambda, alpha, penalty) result(step)
      real(dp), intent(in) :: x(:,:), weights(:), lambda, alpha, penalty(:)
      real(dp) :: step, scale
      integer :: j
      scale = 0.0_dp
      do j = 1, size(x, 2)
         scale = max(scale, sum(weights * x(:, j) ** 2))
      end do
      step = 1.0_dp / max(scale + lambda * (1.0_dp - alpha) * maxval(penalty), 1.0_dp)
   end function initial_step_size

   subroutine initialize_constraints(p, penalty_in, lower_in, upper_in, excluded_in, &
      scale, usable, penalty, lower, upper, excluded, status)
      integer, intent(in) :: p
      real(dp), intent(in), optional :: penalty_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      real(dp), intent(in) :: scale(:)
      logical, intent(in) :: usable(:)
      real(dp), allocatable, intent(out) :: penalty(:), lower(:), upper(:)
      logical, allocatable, intent(out) :: excluded(:)
      integer, intent(out) :: status
      status = glmnet_success
      allocate(penalty(p), lower(p), upper(p), excluded(p))
      penalty = 1.0_dp
      lower = -glmnet_huge
      upper = glmnet_huge
      excluded = .not. usable
      if (present(penalty_in)) then
         if (size(penalty_in) /= p .or. .not. all_finite_vector(penalty_in) .or. &
             any(penalty_in < 0.0_dp)) then
            status = glmnet_invalid_argument
            return
         end if
         penalty = penalty_in
      end if
      if (present(lower_in)) then
         if (size(lower_in) /= p .or. .not. all_finite_vector(lower_in)) then
            status = glmnet_invalid_argument
            return
         end if
         lower = lower_in * scale
      end if
      if (present(upper_in)) then
         if (size(upper_in) /= p .or. .not. all_finite_vector(upper_in)) then
            status = glmnet_invalid_argument
            return
         end if
         upper = upper_in * scale
      end if
      if (any(lower > upper)) then
         status = glmnet_invalid_argument
         return
      end if
      if (present(excluded_in)) then
         if (size(excluded_in) /= p) then
            status = glmnet_invalid_argument
            return
         end if
         excluded = excluded .or. excluded_in
      end if
   end subroutine initialize_constraints

   subroutine initialize_result(result, n, p)
      type(glmnet_path_result), intent(out) :: result
      integer, intent(in) :: n, p
      result%family_code = glmnet_family_cox
      result%family = family_name(glmnet_family_cox)
      result%status = glmnet_success
      result%nobs = n
      result%nvars = p
      result%nout = 1
      result%nlambda = 0
   end subroutine initialize_result

   subroutine allocate_result(result, nlambda, p)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: nlambda, p
      result%nlambda = nlambda
      allocate(result%lambda(nlambda), result%intercept(1, nlambda), &
         result%beta(p, 1, nlambda), result%dev_ratio(nlambda), &
         result%objective(nlambda), result%df(nlambda), result%iterations(nlambda), &
         result%converged(nlambda))
      result%intercept = 0.0_dp
      result%beta = 0.0_dp
      result%dev_ratio = 0.0_dp
      result%objective = 0.0_dp
      result%df = 0
      result%iterations = 0
      result%converged = .false.
   end subroutine allocate_result

   subroutine truncate_result(result, keep)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: keep
      real(dp), allocatable :: lambda(:), intercept(:,:), beta(:,:,:), dev(:), objective(:)
      integer, allocatable :: df(:), iterations(:)
      logical, allocatable :: converged(:)
      if (keep >= result%nlambda) return
      allocate(lambda(keep), intercept(1, keep), beta(result%nvars, 1, keep), &
         dev(keep), objective(keep), df(keep), iterations(keep), converged(keep))
      lambda = result%lambda(:keep)
      intercept = result%intercept(:, :keep)
      beta = result%beta(:, :, :keep)
      dev = result%dev_ratio(:keep)
      objective = result%objective(:keep)
      df = result%df(:keep)
      iterations = result%iterations(:keep)
      converged = result%converged(:keep)
      call move_alloc(lambda, result%lambda)
      call move_alloc(intercept, result%intercept)
      call move_alloc(beta, result%beta)
      call move_alloc(dev, result%dev_ratio)
      call move_alloc(objective, result%objective)
      call move_alloc(df, result%df)
      call move_alloc(iterations, result%iterations)
      call move_alloc(converged, result%converged)
      result%nlambda = keep
   end subroutine truncate_result
end module glmnet_cox
