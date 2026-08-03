! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp
  use segmented_status
  use segmented_types
  use segmented_utils
  implicit none
  private
  public :: fit_segmented_lm, fit_stepmented_lm
  public :: fit_segmented_glm, fit_stepmented_glm
  public :: segmented_model_matrix
contains
  subroutine segmented_model_matrix(x, z, breakpoints, kind, design, status)
    real(dp), intent(in) :: x(:,:), z(:,:), breakpoints(:)
    integer, intent(in) :: kind
    real(dp), allocatable, intent(out) :: design(:,:)
    integer, intent(out) :: status
    if (size(x, 1) /= size(z, 1) .or. size(z, 2) /= size(breakpoints) .or. &
        size(x, 2) < 1 .or. size(z, 2) < 1) then
      allocate(design(0, 0))
      status = SEG_DIMENSION_ERROR
      return
    end if
    if (kind /= SEGMENTED_CONTINUOUS .and. kind /= SEGMENTED_STEP) then
      allocate(design(0, 0))
      status = SEG_INVALID_ARGUMENT
      return
    end if
    call make_design(x, z, breakpoints, kind, design)
    status = SEG_SUCCESS
  end subroutine segmented_model_matrix

  subroutine fit_segmented_lm(y, x, z, psi0, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    real(dp), allocatable :: psi(:), candidate(:), lower(:), upper(:)
    real(dp), allocatable :: design(:,:), aug(:,:), u(:,:), v(:,:), beta(:), cov(:,:)
    real(dp), allocatable :: fitted(:), residuals(:), w(:), histpsi(:,:), histobj(:)
    real(dp) :: rss, sigma, loglik, current, candidate_obj, step, change, denom
    integer :: n, p, m, status, iter, ls, j
    logical :: accepted

    ctl = segmented_control()
    if (present(control)) ctl = control
    call validate_inputs(y, x, z, psi0, ctl, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    n = size(y)
    p = size(x, 2)
    m = size(z, 2)
    allocate(w(n), psi(m), candidate(m), histpsi(ctl%max_iter + 1, m))
    allocate(histobj(ctl%max_iter + 1))
    w = 1.0_dp
    if (present(weights)) w = weights
    call breakpoint_limits(z, ctl%lower_quantile, ctl%upper_quantile, lower, upper, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    psi = psi0
    call clamp_breakpoints(z, lower, upper, psi)
    call make_design(x, z, psi, SEGMENTED_CONTINUOUS, design)
    call weighted_lm_fit(y, design, beta, cov, fitted, residuals, current, sigma, &
        loglik, status, w, offset)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    histpsi(1, :) = psi
    histobj(1) = current
    do iter = 1, ctl%max_iter
      call hinge_matrix(z, psi, u, v)
      allocate(aug(n, p + 2 * m))
      aug(:, :p) = x
      aug(:, p + 1:p + m) = u
      aug(:, p + m + 1:) = v
      call weighted_lm_fit(y, aug, beta, cov, fitted, residuals, rss, sigma, &
          loglik, status, w, offset)
      deallocate(aug)
      if (status /= SEG_SUCCESS) exit
      candidate = psi
      do j = 1, m
        denom = beta(p + j)
        if (abs(denom) > sqrt(epsilon(1.0_dp))) then
          candidate(j) = psi(j) + beta(p + m + j) / denom
        end if
      end do
      call clamp_breakpoints(z, lower, upper, candidate)
      accepted = .false.
      step = 1.0_dp
      do ls = 1, ctl%max_line_search
        candidate = psi + step * (candidate - psi)
        call clamp_breakpoints(z, lower, upper, candidate)
        call make_design(x, z, candidate, SEGMENTED_CONTINUOUS, design)
        call weighted_lm_fit(y, design, beta, cov, fitted, residuals, candidate_obj, &
            sigma, loglik, status, w, offset)
        if (status == SEG_SUCCESS .and. candidate_obj <= current + &
            100.0_dp * epsilon(1.0_dp) * max(1.0_dp, current)) then
          accepted = .true.
          exit
        end if
        step = 0.5_dp * step
      end do
      if (.not. accepted) then
        status = SEG_SUCCESS
        exit
      end if
      change = maxval(abs(candidate - psi)) / max(1.0_dp, maxval(abs(psi)))
      psi = candidate
      current = candidate_obj
      histpsi(iter + 1, :) = psi
      histobj(iter + 1) = current
      if (ctl%verbose) write(*, '(a,i0,a,es12.4)') 'segmented LM iteration ', iter, &
          ', RSS=', current
      if (change <= ctl%breakpoint_tolerance) exit
    end do
    call finalize_fit(y, x, z, psi, SEGMENTED_CONTINUOUS, FAMILY_GAUSSIAN, &
        result, w, offset, iter, histpsi, histobj)
    if (result%status == SEG_SUCCESS) then
      result%converged = iter <= ctl%max_iter
      if (.not. result%converged) result%status = SEG_MAX_ITER
    end if
  end subroutine fit_segmented_lm

  subroutine fit_stepmented_lm(y, x, z, psi0, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_step_model(y, x, z, psi0, FAMILY_GAUSSIAN, result, weights, offset, control)
  end subroutine fit_stepmented_lm

  subroutine fit_segmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    real(dp), allocatable :: psi(:), candidate(:), lower(:), upper(:), design(:,:), aug(:,:)
    real(dp), allocatable :: u(:,:), v(:,:), beta(:), cov(:,:), fitted(:), residuals(:)
    real(dp), allocatable :: w(:), histpsi(:,:), histobj(:)
    real(dp) :: current, candidate_obj, loglik, step, change, denom
    integer :: n, p, m, status, iter, ls, j
    logical :: accepted

    ctl = segmented_control()
    if (present(control)) ctl = control
    call validate_inputs(y, x, z, psi0, ctl, status)
    if (status /= SEG_SUCCESS .or. family == FAMILY_GAUSSIAN) then
      if (family == FAMILY_GAUSSIAN .and. status == SEG_SUCCESS) then
        call fit_segmented_lm(y, x, z, psi0, result, weights, offset, ctl)
      else
        result%status = status
      end if
      return
    end if
    if (family /= FAMILY_BINOMIAL .and. family /= FAMILY_POISSON) then
      result%status = SEG_INVALID_ARGUMENT
      return
    end if
    n = size(y)
    p = size(x, 2)
    m = size(z, 2)
    allocate(w(n), psi(m), candidate(m), histpsi(ctl%max_iter + 1, m))
    allocate(histobj(ctl%max_iter + 1))
    w = 1.0_dp
    if (present(weights)) w = weights
    call breakpoint_limits(z, ctl%lower_quantile, ctl%upper_quantile, lower, upper, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    psi = psi0
    call clamp_breakpoints(z, lower, upper, psi)
    call make_design(x, z, psi, SEGMENTED_CONTINUOUS, design)
    call glm_fit_fixed(y, design, family, beta, cov, fitted, residuals, current, &
        loglik, status, w, offset, ctl%glm_max_iter, ctl%tolerance)
    if (status /= SEG_SUCCESS .and. status /= SEG_MAX_ITER) then
      result%status = status
      return
    end if
    histpsi(1, :) = psi
    histobj(1) = current
    do iter = 1, ctl%max_iter
      call hinge_matrix(z, psi, u, v)
      allocate(aug(n, p + 2 * m))
      aug(:, :p) = x
      aug(:, p + 1:p + m) = u
      aug(:, p + m + 1:) = v
      call glm_fit_fixed(y, aug, family, beta, cov, fitted, residuals, candidate_obj, &
          loglik, status, w, offset, ctl%glm_max_iter, ctl%tolerance)
      deallocate(aug)
      if (status /= SEG_SUCCESS .and. status /= SEG_MAX_ITER) exit
      candidate = psi
      do j = 1, m
        denom = beta(p + j)
        if (abs(denom) > sqrt(epsilon(1.0_dp))) then
          candidate(j) = psi(j) + beta(p + m + j) / denom
        end if
      end do
      call clamp_breakpoints(z, lower, upper, candidate)
      accepted = .false.
      step = 1.0_dp
      do ls = 1, ctl%max_line_search
        candidate = psi + step * (candidate - psi)
        call clamp_breakpoints(z, lower, upper, candidate)
        call make_design(x, z, candidate, SEGMENTED_CONTINUOUS, design)
        call glm_fit_fixed(y, design, family, beta, cov, fitted, residuals, &
            candidate_obj, loglik, status, w, offset, ctl%glm_max_iter, ctl%tolerance)
        if ((status == SEG_SUCCESS .or. status == SEG_MAX_ITER) .and. &
            candidate_obj <= current + 1.0e-10_dp * max(1.0_dp, current)) then
          accepted = .true.
          exit
        end if
        step = 0.5_dp * step
      end do
      if (.not. accepted) then
        status = SEG_SUCCESS
        exit
      end if
      change = maxval(abs(candidate - psi)) / max(1.0_dp, maxval(abs(psi)))
      psi = candidate
      current = candidate_obj
      histpsi(iter + 1, :) = psi
      histobj(iter + 1) = current
      if (ctl%verbose) write(*, '(a,i0,a,es12.4)') 'segmented GLM iteration ', &
          iter, ', deviance=', current
      if (change <= ctl%breakpoint_tolerance) exit
    end do
    call finalize_fit(y, x, z, psi, SEGMENTED_CONTINUOUS, family, result, w, &
        offset, iter, histpsi, histobj, ctl)
    if (result%status == SEG_SUCCESS) then
      result%converged = iter <= ctl%max_iter
      if (.not. result%converged) result%status = SEG_MAX_ITER
    end if
  end subroutine fit_segmented_glm

  subroutine fit_stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_step_model(y, x, z, psi0, family, result, weights, offset, control)
  end subroutine fit_stepmented_glm

  subroutine fit_step_model(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    real(dp), allocatable :: psi(:), trial(:), lower(:), upper(:), design(:,:), sorted_z(:)
    real(dp), allocatable :: beta(:), cov(:,:), fitted(:), residuals(:), w(:)
    real(dp), allocatable :: histpsi(:,:), histobj(:)
    real(dp) :: current, objective, loglik, best, q, change, oldpsi
    integer :: n, m, status, iter, j, g, ng, idx
    logical :: improved

    ctl = segmented_control()
    if (present(control)) ctl = control
    call validate_inputs(y, x, z, psi0, ctl, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    if (family < FAMILY_GAUSSIAN .or. family > FAMILY_POISSON) then
      result%status = SEG_INVALID_ARGUMENT
      return
    end if
    n = size(y)
    m = size(z, 2)
    allocate(w(n), psi(m), trial(m), histpsi(ctl%max_iter + 1, m))
    allocate(histobj(ctl%max_iter + 1))
    w = 1.0_dp
    if (present(weights)) w = weights
    call breakpoint_limits(z, ctl%lower_quantile, ctl%upper_quantile, lower, upper, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    psi = psi0
    call clamp_breakpoints(z, lower, upper, psi)
    call make_design(x, z, psi, SEGMENTED_STEP, design)
    call evaluate_fixed(y, design, family, w, offset, ctl, beta, cov, fitted, &
        residuals, current, loglik, status)
    if (status /= SEG_SUCCESS .and. status /= SEG_MAX_ITER) then
      result%status = status
      return
    end if
    histpsi(1, :) = psi
    histobj(1) = current
    ng = min(max(12, ctl%grid_points), n)
    do iter = 1, ctl%max_iter
      improved = .false.
      change = 0.0_dp
      do j = 1, m
        oldpsi = psi(j)
        best = current
        trial = psi
        allocate(sorted_z(n))
        sorted_z = z(:, j)
        call sort_real(sorted_z)
        do g = 1, ng
          idx = 1 + int(real(g - 1, dp) / real(max(1, ng - 1), dp) * real(n - 2, dp))
          q = 0.5_dp * (sorted_z(idx) + sorted_z(idx + 1))
          q = min(upper(j), max(lower(j), q))
          trial(j) = q
          call clamp_breakpoints(z, lower, upper, trial)
          call make_design(x, z, trial, SEGMENTED_STEP, design)
          call evaluate_fixed(y, design, family, w, offset, ctl, beta, cov, &
              fitted, residuals, objective, loglik, status)
          if ((status == SEG_SUCCESS .or. status == SEG_MAX_ITER) .and. &
              objective < best - 1.0e-10_dp * max(1.0_dp, abs(best))) then
            best = objective
            psi(j) = trial(j)
            improved = .true.
          end if
          trial = psi
        end do
        deallocate(sorted_z)
        change = max(change, abs(psi(j) - oldpsi) / max(1.0_dp, abs(oldpsi)))
        current = best
      end do
      histpsi(iter + 1, :) = psi
      histobj(iter + 1) = current
      if (ctl%verbose) write(*, '(a,i0,a,es12.4)') 'stepmented iteration ', iter, &
          ', objective=', current
      if (.not. improved .or. change <= ctl%breakpoint_tolerance) exit
    end do
    call finalize_fit(y, x, z, psi, SEGMENTED_STEP, family, result, w, offset, &
        iter, histpsi, histobj, ctl)
    if (result%status == SEG_SUCCESS) then
      result%converged = iter <= ctl%max_iter
      if (.not. result%converged) result%status = SEG_MAX_ITER
    end if
  end subroutine fit_step_model

  subroutine evaluate_fixed(y, design, family, weights, offset, ctl, beta, covariance, &
      fitted, residuals, objective, loglik, status)
    real(dp), intent(in) :: y(:), design(:,:), weights(:)
    integer, intent(in) :: family
    real(dp), intent(in), optional :: offset(:)
    type(segmented_control), intent(in) :: ctl
    real(dp), allocatable, intent(out) :: beta(:), covariance(:,:), fitted(:), residuals(:)
    real(dp), intent(out) :: objective, loglik
    integer, intent(out) :: status
    real(dp) :: sigma
    if (family == FAMILY_GAUSSIAN) then
      call weighted_lm_fit(y, design, beta, covariance, fitted, residuals, objective, &
          sigma, loglik, status, weights, offset)
    else
      call glm_fit_fixed(y, design, family, beta, covariance, fitted, residuals, &
          objective, loglik, status, weights, offset, ctl%glm_max_iter, ctl%tolerance)
    end if
  end subroutine evaluate_fixed

  subroutine finalize_fit(y, x, z, psi, kind, family, result, weights, offset, &
      iterations, histpsi, histobj, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi(:), weights(:)
    integer, intent(in) :: kind, family, iterations
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: offset(:)
    real(dp), intent(in) :: histpsi(:,:), histobj(:)
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    real(dp), allocatable :: design(:,:), aug(:,:), u(:,:), v(:,:)
    real(dp), allocatable :: beta(:), cov(:,:), fitted(:), residuals(:)
    real(dp), allocatable :: abeta(:), acov(:,:), afit(:), ares(:)
    real(dp) :: objective, loglik, aobj, allog, denom
    integer :: p, m, n, status, used, j, k
    ctl = segmented_control()
    if (present(control)) ctl = control
    n = size(y)
    p = size(x, 2)
    m = size(z, 2)
    call make_design(x, z, psi, kind, design)
    call evaluate_fixed(y, design, family, weights, offset, ctl, beta, cov, fitted, &
        residuals, objective, loglik, status)
    result%status = status
    if (status /= SEG_SUCCESS .and. status /= SEG_MAX_ITER) return
    result%coefficients = beta
    result%covariance = cov
    result%fitted = fitted
    result%residuals = residuals
    result%weights = weights
    result%breakpoints = psi
    allocate(result%breakpoint_se(m), source = 0.0_dp)
    result%objective = objective
    result%log_likelihood = loglik
    result%family = family
    result%kind = kind
    result%n_base = p
    result%n_break = m
    result%iterations = min(iterations, size(histobj) - 1)
    result%sigma = sqrt(max(0.0_dp, sum(weights * residuals**2) / &
        real(max(1, n - size(beta)), dp)))
    k = size(beta) + m
    result%aic = -2.0_dp * loglik + 2.0_dp * real(k, dp)
    result%bic = -2.0_dp * loglik + log(real(n, dp)) * real(k, dp)
    if (kind == SEGMENTED_CONTINUOUS) then
      call hinge_matrix(z, psi, u, v)
      allocate(aug(n, p + 2 * m))
      aug(:, :p) = x
      aug(:, p + 1:p + m) = u
      aug(:, p + m + 1:) = v
      call evaluate_fixed(y, aug, family, weights, offset, ctl, abeta, acov, afit, &
          ares, aobj, allog, status)
      if (status == SEG_SUCCESS .or. status == SEG_MAX_ITER) then
        do j = 1, m
          denom = abs(abeta(p + j))
          if (denom > sqrt(epsilon(1.0_dp))) then
            result%breakpoint_se(j) = sqrt(max(0.0_dp, acov(p + m + j, p + m + j))) / denom
          else
            result%breakpoint_se(j) = huge(1.0_dp)
          end if
        end do
      end if
    else
      call step_breakpoint_se(y, x, z, psi, family, weights, offset, ctl, &
          result%breakpoint_se)
    end if
    used = min(size(histobj), max(1, result%iterations + 1))
    allocate(result%history_breakpoints(used, m), result%history_objective(used))
    result%history_breakpoints = histpsi(:used, :)
    result%history_objective = histobj(:used)
    result%status = SEG_SUCCESS
  end subroutine finalize_fit

  subroutine step_breakpoint_se(y, x, z, psi, family, weights, offset, ctl, se)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi(:), weights(:)
    integer, intent(in) :: family
    real(dp), intent(in), optional :: offset(:)
    type(segmented_control), intent(in) :: ctl
    real(dp), intent(out) :: se(:)
    real(dp), allocatable :: trial(:), design(:,:), beta(:), cov(:,:), fitted(:), residuals(:)
    real(dp) :: base, loglik, objl, objr, h, curvature
    integer :: j, status
    call make_design(x, z, psi, SEGMENTED_STEP, design)
    call evaluate_fixed(y, design, family, weights, offset, ctl, beta, cov, fitted, &
        residuals, base, loglik, status)
    allocate(trial(size(psi)))
    trial = psi
    do j = 1, size(psi)
      h = max(1.0e-4_dp, 0.01_dp * (maxval(z(:, j)) - minval(z(:, j))))
      trial = psi
      trial(j) = psi(j) - h
      call make_design(x, z, trial, SEGMENTED_STEP, design)
      call evaluate_fixed(y, design, family, weights, offset, ctl, beta, cov, fitted, &
          residuals, objl, loglik, status)
      trial(j) = psi(j) + h
      call make_design(x, z, trial, SEGMENTED_STEP, design)
      call evaluate_fixed(y, design, family, weights, offset, ctl, beta, cov, fitted, &
          residuals, objr, loglik, status)
      curvature = max(0.0_dp, (objl - 2.0_dp * base + objr) / h**2)
      if (curvature > epsilon(1.0_dp)) then
        se(j) = sqrt(2.0_dp / curvature)
      else
        se(j) = huge(1.0_dp)
      end if
    end do
  end subroutine step_breakpoint_se

  subroutine validate_inputs(y, x, z, psi0, ctl, status)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    type(segmented_control), intent(in) :: ctl
    integer, intent(out) :: status
    if (size(y) < 5 .or. size(x, 1) /= size(y) .or. size(z, 1) /= size(y) .or. &
        size(z, 2) /= size(psi0) .or. size(x, 2) < 1 .or. size(z, 2) < 1 .or. &
        size(y) <= size(x, 2) + 2 * size(z, 2)) then
      status = SEG_DIMENSION_ERROR
    else if (ctl%max_iter < 1 .or. ctl%lower_quantile < 0.0_dp .or. &
        ctl%upper_quantile > 1.0_dp .or. ctl%lower_quantile >= ctl%upper_quantile) then
      status = SEG_INVALID_ARGUMENT
    else if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(x)) .or. &
        .not. all(ieee_is_finite(z)) .or. .not. all(ieee_is_finite(psi0))) then
      status = SEG_NONFINITE
    else
      status = SEG_SUCCESS
    end if
  end subroutine validate_inputs
end module segmented_fit
