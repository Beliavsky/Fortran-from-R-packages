! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_regression
  use robstattm_kinds, only : dp
  use robstattm_types, only : robstattm_control, regression_result, linear_test_result, &
    model_selection_result, &
    robstattm_success, robstattm_invalid_argument, robstattm_singular, robstattm_no_convergence
  use robstattm_psi, only : rho_value, rho_prime, rho_second, rho_weight, scale_m, &
    inverse_robust_r_squared, tuning_for_efficiency
  use robstattm_utils, only : matrix_inverse, mean_value, outer_product, lower_string
  use robustbase_lmrob, only : lmrob_result, lmrob_s_fit, lmrob_lar_fit
  use robustbase_linalg, only : least_squares
  use rrcov_stats, only : chi_square_cdf, f_cdf
  implicit none
  private
  public :: lmrobdet_control, lmrobm_control, mm_py_fit, sm_py_fit, lmrob_m_fit
  public :: refine_sm, dcml_fit, dcml_covariance, robust_rfpe, robust_linear_test
  public :: robust_lar_fit, least_squares_fit
  public :: lmrobdet_mm, lmrobdet_dcml, lmrobdet_lin_test, lmrobdet_mm_rfpe
  public :: stepwise_rfpe
contains
  function lmrobdet_control(bb, efficiency, family, max_iter, tolerance, n_resample) result(control)
    real(dp), intent(in), optional :: bb, efficiency, tolerance
    character(len=*), intent(in), optional :: family
    integer, intent(in), optional :: max_iter, n_resample
    type(robstattm_control) :: control
    if (present(bb)) control%bb = min(0.5_dp, max(1.0e-4_dp, bb))
    if (present(efficiency)) control%efficiency = min(0.9999_dp, max(0.05_dp, efficiency))
    if (present(family)) control%family = lower_string(trim(adjustl(family)))
    if (present(max_iter)) control%max_iter = max(1, max_iter)
    if (present(tolerance)) control%tolerance = max(tolerance, epsilon(1.0_dp))
    if (present(n_resample)) control%n_resample = max(1, n_resample)
    control%tuning_psi = tuning_for_efficiency(control%efficiency, control%family)
    select case (trim(control%family))
    case ('bisquare', 'tukey')
      control%tuning_chi = 1.54764_dp
    case ('huber')
      control%tuning_chi = max(0.5_dp, tuning_for_efficiency(0.7_dp, 'huber'))
    case default
      control%tuning_chi = max(0.5_dp, tuning_for_efficiency(0.7_dp, control%family))
    end select
  end function lmrobdet_control

  function lmrobm_control(bb, efficiency, family, max_iter, tolerance) result(control)
    real(dp), intent(in), optional :: bb, efficiency, tolerance
    character(len=*), intent(in), optional :: family
    integer, intent(in), optional :: max_iter
    type(robstattm_control) :: control
    control = lmrobdet_control(bb=bb, efficiency=efficiency, family=family, &
      max_iter=max_iter, tolerance=tolerance)
  end function lmrobm_control

  subroutine least_squares_fit(x, y, result)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    real(dp), allocatable :: beta(:), xtx(:, :), inverse_xtx(:, :)
    real(dp) :: sse, centered_sse, mean_y
    integer :: n, p, status, i
    n = size(x, 1)
    p = size(x, 2)
    if (size(y) /= n .or. n < p .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    allocate(beta(p))
    call least_squares(x, y, beta, status)
    if (status /= 0) then
      result%status = robstattm_singular
      return
    end if
    result%coefficients = beta
    result%fitted = matmul(x, beta)
    result%residuals = y - result%fitted
    allocate(result%weights(n))
    result%weights = 1.0_dp
    sse = sum(result%residuals ** 2)
    result%scale = sqrt(sse / real(max(1, n - p), dp))
    xtx = matmul(transpose(x), x)
    inverse_xtx = matrix_inverse(xtx, status)
    allocate(result%covariance(p, p), result%standard_errors(p))
    if (status == 0) then
      result%covariance = result%scale ** 2 * inverse_xtx
      do i = 1, p
        result%standard_errors(i) = sqrt(max(result%covariance(i, i), 0.0_dp))
      end do
    else
      result%covariance = 0.0_dp
      result%standard_errors = 0.0_dp
    end if
    mean_y = sum(y) / real(n, dp)
    centered_sse = sum((y - mean_y) ** 2)
    if (centered_sse > tiny(1.0_dp)) result%r_squared = 1.0_dp - sse / centered_sse
    result%adjusted_r_squared = adjusted_r_squared(result%r_squared, n, p)
    result%objective = sse
    result%rank = p
    result%iterations = 1
    result%converged = status == 0
    result%status = merge(robstattm_success, robstattm_singular, status == 0)
    result%method = 'least squares'
  end subroutine least_squares_fit

  subroutine lmrobdet_mm(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call mm_py_fit(x, y, result, control)
    result%method = 'lmrobdetMM'
  end subroutine lmrobdet_mm

  subroutine lmrobdet_dcml(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    type(regression_result) :: robust_fit, ls_fit
    call mm_py_fit(x, y, robust_fit, control)
    if (.not. allocated(robust_fit%coefficients)) then
      result%status = robust_fit%status
      return
    end if
    call least_squares_fit(x, y, ls_fit)
    if (.not. allocated(ls_fit%coefficients)) then
      result%status = ls_fit%status
      return
    end if
    call dcml_fit(x, y, robust_fit, ls_fit, result, control)
    result%method = 'lmrobdetDCML'
  end subroutine lmrobdet_dcml

  subroutine lmrobdet_lin_test(full_fit, restricted_fit, result, control)
    type(regression_result), intent(in) :: full_fit, restricted_fit
    type(linear_test_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call robust_linear_test(full_fit, restricted_fit, control, result)
  end subroutine lmrobdet_lin_test

  function lmrobdet_mm_rfpe(fit, control, scale, minimum_rho, penalty) result(value)
    type(regression_result), intent(in) :: fit
    type(robstattm_control), intent(in), optional :: control
    real(dp), intent(in), optional :: scale
    real(dp), intent(out), optional :: minimum_rho, penalty
    real(dp) :: value
    value = robust_rfpe(fit, control, scale, minimum_rho, penalty)
  end function lmrobdet_mm_rfpe

  subroutine stepwise_rfpe(x, y, result, control, direction, force_first, max_steps)
    real(dp), intent(in) :: x(:, :), y(:)
    type(model_selection_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    character(len=*), intent(in), optional :: direction
    logical, intent(in), optional :: force_first
    integer, intent(in), optional :: max_steps
    type(robstattm_control) :: ctl
    type(regression_result) :: current_fit, candidate_fit
    logical, allocatable :: selected(:), candidate_mask(:)
    real(dp), allocatable :: xsub(:, :), history(:)
    real(dp) :: current_value, candidate_value, best_value, common_scale
    integer :: p, j, best_column, step, limit
    logical :: protect_first, changed
    character(len=16) :: mode

    p = size(x, 2)
    if (size(y) /= size(x, 1) .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    mode = 'backward'
    if (present(direction)) mode = lower_string(trim(adjustl(direction)))
    protect_first = .true.
    if (present(force_first)) protect_first = force_first
    limit = p
    if (present(max_steps)) limit = max(1, max_steps)
    allocate(selected(p), candidate_mask(p), history(limit + 1))
    history = huge(1.0_dp)

    if (trim(mode) == 'forward') then
      selected = .false.
      if (protect_first) then
        selected(1) = .true.
      else
        selected(1) = .true.
      end if
    else
      selected = .true.
    end if
    call select_columns(x, selected, xsub)
    call mm_py_fit(xsub, y, current_fit, ctl)
    if (.not. allocated(current_fit%coefficients)) then
      result%status = current_fit%status
      return
    end if
    common_scale = current_fit%scale
    current_value = robust_rfpe(current_fit, ctl, scale=common_scale)
    history(1) = current_value

    do step = 1, limit
      best_value = current_value
      best_column = 0
      changed = .false.

      if (trim(mode) /= 'forward') then
        do j = 1, p
          if (.not. selected(j)) cycle
          if (protect_first .and. j == 1) cycle
          if (count(selected) <= 1) cycle
          candidate_mask = selected
          candidate_mask(j) = .false.
          call select_columns(x, candidate_mask, xsub)
          call mm_py_fit(xsub, y, candidate_fit, ctl)
          if (.not. allocated(candidate_fit%coefficients)) cycle
          candidate_value = robust_rfpe(candidate_fit, ctl, scale=common_scale)
          if (candidate_value < best_value - 1.0e-10_dp) then
            best_value = candidate_value
            best_column = -j
          end if
        end do
      end if

      if (trim(mode) /= 'backward') then
        do j = 1, p
          if (selected(j)) cycle
          candidate_mask = selected
          candidate_mask(j) = .true.
          call select_columns(x, candidate_mask, xsub)
          call mm_py_fit(xsub, y, candidate_fit, ctl)
          if (.not. allocated(candidate_fit%coefficients)) cycle
          candidate_value = robust_rfpe(candidate_fit, ctl, scale=common_scale)
          if (candidate_value < best_value - 1.0e-10_dp) then
            best_value = candidate_value
            best_column = j
          end if
        end do
      end if

      if (best_column < 0) then
        selected(-best_column) = .false.
        changed = .true.
      else if (best_column > 0) then
        selected(best_column) = .true.
        changed = .true.
      end if
      if (.not. changed) exit
      current_value = best_value
      history(step + 1) = current_value
      result%iterations = step
    end do

    allocate(result%selected_columns(count(selected)))
    result%selected_columns = pack([(j, j=1,p)], selected)
    result%criterion = current_value
    result%criterion_history = history(1:result%iterations + 1)
    result%converged = .true.
    result%status = robstattm_success
    result%direction = trim(mode)
  end subroutine stepwise_rfpe

  subroutine mm_py_fit(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    type(robstattm_control) :: ctl
    type(lmrob_result) :: initial
    real(dp), allocatable :: beta(:), weights(:), residuals(:)
    integer :: iteration, info, n, p
    real(dp) :: delta
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    n = size(x, 1)
    p = size(x, 2)
    if (size(y) /= n .or. n < p .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    call lmrob_s_fit(x, y, initial, n_resample=ctl%n_resample, sampling='nonsingular', &
      tuning_chi=ctl%tuning_chi, bb=corrected_breakdown(ctl, n, p), &
      max_refine=ctl%refine_iter, tol=ctl%tolerance)
    allocate(beta(p), weights(n), residuals(n))
    beta = initial%coefficients
    call fixed_scale_m_fit(x, y, beta, max(initial%scale, 1.0e-14_dp), ctl, &
      weights, iteration, info)
    residuals = y - matmul(x, beta)
    call fill_regression_result(x, beta, residuals, weights, initial%scale, &
      iteration + initial%iterations, info == 0 .and. initial%converged, ctl, result)
    result%method = 'MM-PY'
    delta = robust_r_squared_value(y, residuals, result%scale, ctl)
    result%r_squared = delta
    result%adjusted_r_squared = adjusted_r_squared(delta, n, p)
  end subroutine mm_py_fit

  subroutine sm_py_fit(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call mm_py_fit(x, y, result, control)
    result%method = 'SM-PY'
  end subroutine sm_py_fit

  subroutine lmrob_m_fit(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call mm_py_fit(x, y, result, control)
    result%method = 'M'
  end subroutine lmrob_m_fit

  subroutine robust_lar_fit(x, y, result, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    type(lmrob_result) :: fit
    integer :: mi
    real(dp) :: tt
    mi = 1000
    if (present(max_iter)) mi = max_iter
    tt = 1.0e-8_dp
    if (present(tol)) tt = tol
    call lmrob_lar_fit(x, y, fit, max_iter=mi, tol=tt)
    call copy_lmrob_result(fit, result)
    result%method = 'LAR'
  end subroutine robust_lar_fit

  subroutine refine_sm(x, y, initial_beta, initial_scale, result, control, step, max_iter)
    real(dp), intent(in) :: x(:, :), y(:), initial_beta(:), initial_scale
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    character(len=*), intent(in), optional :: step
    integer, intent(in), optional :: max_iter
    type(robstattm_control) :: ctl
    real(dp), allocatable :: beta(:), weights(:), residuals(:)
    real(dp) :: scale
    integer :: mi, iteration, info
    character(len=8) :: mode
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    mi = ctl%refine_iter
    if (present(max_iter)) mi = max_iter
    mode = 'M'
    if (present(step)) mode = step
    beta = initial_beta
    scale = initial_scale
    allocate(weights(size(y)), residuals(size(y)))
    if (trim(mode) == 'S') then
      call refine_s_scale(x, y, beta, scale, ctl, weights, iteration, info, mi)
    else
      call fixed_scale_m_fit(x, y, beta, max(scale, 1.0e-14_dp), ctl, weights, iteration, info, mi)
    end if
    residuals = y - matmul(x, beta)
    call fill_regression_result(x, beta, residuals, weights, scale, iteration, info == 0, ctl, result)
    result%method = 'refine-' // trim(mode)
  end subroutine refine_sm

  subroutine dcml_fit(x, y, robust_fit, ls_fit, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(in) :: robust_fit, ls_fit
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    type(robstattm_control) :: ctl
    real(dp), allocatable :: cc_matrix(:, :), difference(:)
    real(dp) :: robust_scale, d, delta_s, mixing, sw
    integer :: n, p, info, i
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    n = size(x, 1)
    p = size(x, 2)
    if (.not. allocated(robust_fit%coefficients) .or. .not. allocated(ls_fit%coefficients)) then
      result%status = robstattm_invalid_argument
      return
    end if
    robust_scale = scale_m(robust_fit%residuals, corrected_breakdown(ctl, n, p), &
      ctl%family, ctl%tuning_chi, ctl%max_iter, ctl%tolerance)
    sw = max(sum(robust_fit%weights), tiny(1.0_dp))
    cc_matrix = matmul(transpose(x), x * spread(robust_fit%weights, 2, p)) / sw
    difference = robust_fit%coefficients - ls_fit%coefficients
    d = dot_product(difference, matmul(cc_matrix, difference)) / max(robust_scale ** 2, tiny(1.0_dp))
    delta_s = 0.3_dp * real(p, dp) / real(n, dp)
    if (d <= tiny(1.0_dp)) then
      mixing = 1.0_dp
    else
      mixing = min(1.0_dp, sqrt(delta_s / d))
    end if
    allocate(result%coefficients(p), result%fitted(n), result%residuals(n), result%weights(n))
    result%coefficients = mixing * ls_fit%coefficients + (1.0_dp - mixing) * robust_fit%coefficients
    result%fitted = matmul(x, result%coefficients)
    result%residuals = y - result%fitted
    result%weights = robust_fit%weights
    result%scale = scale_m(result%residuals, corrected_breakdown(ctl, n, p), ctl%family, &
      ctl%tuning_chi, ctl%max_iter, ctl%tolerance)
    call dcml_covariance(ls_fit%residuals, robust_fit%residuals, cc_matrix, robust_scale, &
      mixing, p, n, ctl, result%covariance, info)
    result%covariance = result%covariance / real(n, dp)
    allocate(result%standard_errors(p))
    do i = 1, p
      result%standard_errors(i) = sqrt(max(result%covariance(i, i), 0.0_dp))
    end do
    result%mixing = mixing
    result%rank = p
    result%iterations = robust_fit%iterations
    result%converged = robust_fit%converged
    result%status = merge(robstattm_success, robstattm_singular, info == 0)
    result%method = 'DCML'
  end subroutine dcml_fit

  subroutine dcml_covariance(residuals_ls, residuals_robust, cc_matrix, robust_scale, mixing, &
      p, n, control, covariance, status)
    real(dp), intent(in) :: residuals_ls(:), residuals_robust(:), cc_matrix(:, :)
    real(dp), intent(in) :: robust_scale, mixing
    integer, intent(in) :: p, n
    type(robstattm_control), intent(in) :: control
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: rr(:), psi(:), psip(:), inverse_cc(:, :)
    real(dp) :: t, c0, a1, b0, scale_ls, value
    integer :: i
    allocate(rr(size(residuals_robust)), psi(size(residuals_robust)), psip(size(residuals_robust)))
    rr = residuals_robust / max(robust_scale, tiny(1.0_dp))
    do i = 1, size(rr)
      psi(i) = rho_prime(rr(i), control%family, control%tuning_psi)
      psip(i) = rho_second(rr(i), control%family, control%tuning_psi)
    end do
    t = 1.0_dp - mixing
    c0 = mean_value(psi * residuals_ls)
    a1 = mean_value(psi * psi)
    b0 = mean_value(psip)
    scale_ls = scale_m(residuals_ls, corrected_breakdown(control, n, p), control%family, &
      control%tuning_chi, control%max_iter, control%tolerance)
    value = t * t * robust_scale * robust_scale * a1 / max(b0 * b0, tiny(1.0_dp)) + &
      scale_ls * scale_ls * (1.0_dp - t) ** 2 + &
      2.0_dp * t * (1.0_dp - t) * robust_scale * c0 / max(b0, sqrt(tiny(1.0_dp)))
    inverse_cc = matrix_inverse(cc_matrix, status)
    if (status == 0) then
      covariance = value * inverse_cc
    else
      allocate(covariance(p, p))
      covariance = 0.0_dp
    end if
  end subroutine dcml_covariance

  function robust_rfpe(fit, control, scale, minimum_rho, penalty) result(value)
    type(regression_result), intent(in) :: fit
    type(robstattm_control), intent(in), optional :: control
    real(dp), intent(in), optional :: scale
    real(dp), intent(out), optional :: minimum_rho, penalty
    real(dp) :: value, s, a2, b2, d2
    type(robstattm_control) :: ctl
    real(dp), allocatable :: u(:), psi(:), psip(:), rho_values(:)
    integer :: i, p, n
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    s = fit%scale
    if (present(scale)) s = scale
    n = size(fit%residuals)
    p = size(fit%coefficients)
    allocate(u(n), psi(n), psip(n), rho_values(n))
    u = fit%residuals / max(s, tiny(1.0_dp))
    do i = 1, n
      rho_values(i) = rho_value(u(i), ctl%family, ctl%tuning_psi, .true.)
      psi(i) = rho_prime(u(i), ctl%family, ctl%tuning_psi, .true.)
      psip(i) = rho_second(u(i), ctl%family, ctl%tuning_psi, .true.)
    end do
    a2 = mean_value(rho_values)
    b2 = real(p, dp) * mean_value(psi * psi)
    d2 = mean_value(psip)
    if (d2 <= 0.0_dp) then
      value = huge(1.0_dp)
    else
      value = a2 + b2 / (d2 * real(n, dp))
    end if
    if (present(minimum_rho)) minimum_rho = a2
    if (present(penalty)) penalty = b2 / max(d2 * real(n, dp), tiny(1.0_dp))
  end function robust_rfpe

  subroutine robust_linear_test(full_fit, restricted_fit, control, result)
    type(regression_result), intent(in) :: full_fit, restricted_fit
    type(robstattm_control), intent(in), optional :: control
    type(linear_test_result), intent(out) :: result
    type(robstattm_control) :: ctl
    real(dp), allocatable :: u_full(:), u_restricted(:), a(:), b(:), c(:), d(:)
    real(dp) :: scale
    integer :: i, n, p, q
    ctl = lmrobdet_control()
    if (present(control)) ctl = control
    n = size(full_fit%residuals)
    p = size(full_fit%coefficients)
    q = p - size(restricted_fit%coefficients)
    if (q <= 0 .or. n <= p) then
      result%status = robstattm_invalid_argument
      return
    end if
    scale = max(full_fit%scale, tiny(1.0_dp))
    u_full = full_fit%residuals / scale
    u_restricted = restricted_fit%residuals / scale
    allocate(a(n), b(n), c(n), d(n))
    do i = 1, n
      a(i) = rho_value(u_restricted(i), ctl%family, ctl%tuning_psi, .true.)
      b(i) = rho_value(u_full(i), ctl%family, ctl%tuning_psi, .true.)
      c(i) = rho_second(u_full(i), ctl%family, ctl%tuning_psi, .true.)
      d(i) = rho_prime(u_full(i), ctl%family, ctl%tuning_psi, .true.) ** 2
    end do
    result%statistic = 2.0_dp * (sum(a) - sum(b)) * sum(c) / max(sum(d), tiny(1.0_dp))
    result%df1 = q
    result%df2 = n - p
    result%chi_square_p_value = 1.0_dp - chi_square_cdf(max(result%statistic, 0.0_dp), real(q, dp))
    result%f_p_value = 1.0_dp - f_cdf(max(result%statistic / real(q, dp), 0.0_dp), &
      real(q, dp), real(n - p, dp))
    result%status = robstattm_success
  end subroutine robust_linear_test

  subroutine fixed_scale_m_fit(x, y, beta, scale, control, weights, iterations, status, max_iter)
    real(dp), intent(in) :: x(:, :), y(:), scale
    real(dp), intent(inout) :: beta(:)
    type(robstattm_control), intent(in) :: control
    real(dp), intent(out) :: weights(:)
    integer, intent(out) :: iterations, status
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: residuals(:), new_beta(:), xw(:, :), yw(:)
    real(dp) :: delta
    integer :: i, j, mi, info
    mi = control%max_iter
    if (present(max_iter)) mi = max_iter
    allocate(residuals(size(y)), new_beta(size(beta)), xw(size(x,1),size(x,2)), yw(size(y)))
    status = 0
    do iterations = 1, mi
      residuals = y - matmul(x, beta)
      do i = 1, size(y)
        weights(i) = rho_weight(residuals(i) / scale, control%family, control%tuning_psi)
      end do
      if (count(weights > 1.0e-12_dp) < size(beta)) then
        status = 1
        return
      end if
      do j = 1, size(x, 2)
        xw(:, j) = x(:, j) * sqrt(weights)
      end do
      yw = y * sqrt(weights)
      call least_squares(xw, yw, new_beta, info)
      if (info /= 0) then
        status = info
        return
      end if
      delta = maxval(abs(new_beta - beta))
      beta = new_beta
      if (delta <= control%tolerance * (1.0_dp + maxval(abs(beta)))) return
    end do
    status = 1
  end subroutine fixed_scale_m_fit

  subroutine refine_s_scale(x, y, beta, scale, control, weights, iterations, status, max_iter)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), intent(inout) :: beta(:), scale
    type(robstattm_control), intent(in) :: control
    real(dp), intent(out) :: weights(:)
    integer, intent(out) :: iterations, status
    integer, intent(in) :: max_iter
    real(dp), allocatable :: residuals(:)
    integer :: i
    status = 0
    do iterations = 1, max_iter
      residuals = y - matmul(x, beta)
      scale = scale_m(residuals, control%bb, control%family, control%tuning_chi, &
        control%max_iter, control%tolerance)
      call fixed_scale_m_fit(x, y, beta, max(scale, 1.0e-14_dp), control, weights, i, status, 1)
      if (status /= 0) return
    end do
  end subroutine refine_s_scale

  subroutine fill_regression_result(x, beta, residuals, weights, scale, iterations, converged, control, result)
    real(dp), intent(in) :: x(:, :), beta(:), residuals(:), weights(:), scale
    integer, intent(in) :: iterations
    logical, intent(in) :: converged
    type(robstattm_control), intent(in) :: control
    type(regression_result), intent(out) :: result
    real(dp), allocatable :: a(:, :), b(:, :), ainv(:, :), psi(:), psip(:)
    integer :: i, n, p, status
    n = size(x, 1)
    p = size(x, 2)
    allocate(result%coefficients(p), result%fitted(n), result%residuals(n), result%weights(n))
    result%coefficients = beta
    result%fitted = matmul(x, beta)
    result%residuals = residuals
    result%weights = weights
    result%scale = scale
    result%iterations = iterations
    result%rank = p
    result%converged = converged
    allocate(a(p,p), b(p,p), psi(n), psip(n))
    a = 0.0_dp
    b = 0.0_dp
    do i = 1, n
      psi(i) = rho_prime(residuals(i) / max(scale, tiny(1.0_dp)), control%family, control%tuning_psi)
      psip(i) = rho_second(residuals(i) / max(scale, tiny(1.0_dp)), control%family, control%tuning_psi)
      a = a + psip(i) * outer_product(x(i,:), x(i,:))
      b = b + psi(i) * psi(i) * outer_product(x(i,:), x(i,:))
    end do
    a = a / real(n, dp)
    b = b / real(n, dp)
    ainv = matrix_inverse(a, status)
    allocate(result%covariance(p,p), result%standard_errors(p))
    if (status == 0) then
      result%covariance = scale * scale * matmul(matmul(ainv, b), ainv) / real(n, dp)
      do i = 1, p
        result%standard_errors(i) = sqrt(max(result%covariance(i,i), 0.0_dp))
      end do
      result%status = robstattm_success
    else
      result%covariance = 0.0_dp
      result%standard_errors = 0.0_dp
      result%status = robstattm_singular
    end if
    result%objective = sum([(rho_value(residuals(i) / max(scale, tiny(1.0_dp)), &
      control%family, control%tuning_psi, .true.), i=1,n)])
  end subroutine fill_regression_result

  subroutine copy_lmrob_result(source, target)
    type(lmrob_result), intent(in) :: source
    type(regression_result), intent(out) :: target
    target%coefficients = source%coefficients
    target%fitted = source%fitted
    target%residuals = source%residuals
    target%weights = source%weights
    target%covariance = source%covariance
    target%standard_errors = source%standard_errors
    target%scale = source%scale
    target%objective = source%objective
    target%iterations = source%iterations
    target%rank = size(source%coefficients)
    target%converged = source%converged
    target%status = merge(robstattm_success, robstattm_no_convergence, source%converged)
  end subroutine copy_lmrob_result

  subroutine select_columns(x, mask, selected_x)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in) :: mask(:)
    real(dp), allocatable, intent(out) :: selected_x(:, :)
    integer :: j, k
    allocate(selected_x(size(x, 1), count(mask)))
    k = 0
    do j = 1, size(x, 2)
      if (.not. mask(j)) cycle
      k = k + 1
      selected_x(:, k) = x(:, j)
    end do
  end subroutine select_columns

  function corrected_breakdown(control, n, p) result(value)
    type(robstattm_control), intent(in) :: control
    integer, intent(in) :: n, p
    real(dp) :: value
    value = control%bb
    if (control%finite_sample_correction) value = value * max(0.05_dp, 1.0_dp - real(p,dp)/real(n,dp))
  end function corrected_breakdown

  function robust_r_squared_value(y, residuals, scale, control) result(value)
    real(dp), intent(in) :: y(:), residuals(:), scale
    type(robstattm_control), intent(in) :: control
    real(dp) :: value, center, s2, s02, raw
    integer :: i
    center = sum(y) / real(size(y), dp)
    s2 = 0.0_dp
    s02 = 0.0_dp
    do i = 1, size(y)
      s2 = s2 + rho_value(residuals(i)/max(scale,tiny(1.0_dp)), control%family, control%tuning_psi, .true.)
      s02 = s02 + rho_value((y(i)-center)/max(scale,tiny(1.0_dp)), control%family, control%tuning_psi, .true.)
    end do
    s2 = s2 / real(size(y), dp)
    s02 = s02 / real(size(y), dp)
    raw = (s02 - s2) / max(s02 * (1.0_dp - s2), tiny(1.0_dp))
    value = inverse_robust_r_squared(raw, control%family, control%tuning_psi)
  end function robust_r_squared_value

  pure function adjusted_r_squared(r2, n, p) result(value)
    real(dp), intent(in) :: r2
    integer, intent(in) :: n, p
    real(dp) :: value
    if (n <= p) then
      value = r2
    else
      value = real(n-1,dp)/real(n-p,dp)*r2-real(p-1,dp)/real(n-p,dp)
    end if
  end function adjusted_r_squared
end module robstattm_regression
