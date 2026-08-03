! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_mixed
  use nlme_kinds, only : dp
  use nlme_lme, only : fit_lme
  use nlme_status, only : NLME_SUCCESS, NLME_MAX_ITER
  use nlme_types, only : lme_result
  use segmented_status
  use segmented_types
  use segmented_utils, only : breakpoint_limits, clamp_breakpoints, hinge_matrix, make_design
  implicit none
  private
  public :: fit_segmented_lme
contains
  subroutine fit_segmented_lme(y, x, z, psi0, random_design, group, result, options, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:), random_design(:,:)
    integer, intent(in) :: group(:)
    type(segmented_lme_result), intent(out) :: result
    type(segmented_lme_options), intent(in), optional :: options
    type(segmented_control), intent(in), optional :: control
    type(segmented_control) :: ctl
    type(segmented_lme_options) :: opt
    type(lme_result) :: augmented_fit, candidate_fit, final_fit
    real(dp), allocatable :: psi(:), candidate(:), lower(:), upper(:), u(:,:), v(:,:)
    real(dp), allocatable :: design(:,:), augmented(:,:), histpsi(:,:), histobj(:)
    real(dp) :: current, trial_objective, step, change, denom
    integer :: n, p, m, status, iter, ls, j, used
    logical :: accepted

    ctl = segmented_control()
    if (present(control)) ctl = control
    opt = segmented_lme_options()
    if (present(options)) opt = options
    n = size(y)
    p = size(x, 2)
    m = size(z, 2)
    if (size(x, 1) /= n .or. size(z, 1) /= n .or. size(random_design, 1) /= n .or. &
        size(group) /= n .or. size(psi0) /= m .or. p < 1 .or. m < 1 .or. &
        size(random_design, 2) < 1 .or. n <= p + 2 * m) then
      result%status = SEG_DIMENSION_ERROR
      return
    end if
    call breakpoint_limits(z, ctl%lower_quantile, ctl%upper_quantile, lower, upper, status)
    if (status /= SEG_SUCCESS) then
      result%status = status
      return
    end if
    allocate(psi(m), candidate(m), histpsi(ctl%max_iter + 1, m))
    allocate(histobj(ctl%max_iter + 1))
    psi = psi0
    call clamp_breakpoints(z, lower, upper, psi)
    call make_design(x, z, psi, SEGMENTED_CONTINUOUS, design)
    call run_lme(y, design, random_design, group, opt, final_fit)
    if (.not. valid_lme(final_fit)) then
      result%status = SEG_NLME_ERROR
      result%fit = final_fit
      return
    end if
    current = -final_fit%log_likelihood
    histpsi(1, :) = psi
    histobj(1) = current
    do iter = 1, ctl%max_iter
      call hinge_matrix(z, psi, u, v)
      allocate(augmented(n, p + 2 * m))
      augmented(:, :p) = x
      augmented(:, p + 1:p + m) = u
      augmented(:, p + m + 1:) = v
      call run_lme(y, augmented, random_design, group, opt, augmented_fit)
      deallocate(augmented)
      if (.not. valid_lme(augmented_fit)) exit
      candidate = psi
      do j = 1, m
        denom = augmented_fit%beta(p + j)
        if (abs(denom) > sqrt(epsilon(1.0_dp))) then
          candidate(j) = psi(j) + augmented_fit%beta(p + m + j) / denom
        end if
      end do
      call clamp_breakpoints(z, lower, upper, candidate)
      step = 1.0_dp
      accepted = .false.
      do ls = 1, ctl%max_line_search
        candidate = psi + step * (candidate - psi)
        call clamp_breakpoints(z, lower, upper, candidate)
        call make_design(x, z, candidate, SEGMENTED_CONTINUOUS, design)
        call run_lme(y, design, random_design, group, opt, candidate_fit)
        if (valid_lme(candidate_fit)) then
          trial_objective = -candidate_fit%log_likelihood
          if (trial_objective <= current + 1.0e-8_dp * max(1.0_dp, abs(current))) then
            accepted = .true.
            exit
          end if
        end if
        step = 0.5_dp * step
      end do
      if (.not. accepted) exit
      change = maxval(abs(candidate - psi)) / max(1.0_dp, maxval(abs(psi)))
      psi = candidate
      current = trial_objective
      final_fit = candidate_fit
      histpsi(iter + 1, :) = psi
      histobj(iter + 1) = current
      if (ctl%verbose) write(*, '(a,i0,a,es12.4)') 'segmented LME iteration ', iter, &
          ', negative logLik=', current
      if (change <= ctl%breakpoint_tolerance) exit
    end do
    call make_design(x, z, psi, SEGMENTED_CONTINUOUS, design)
    call run_lme(y, design, random_design, group, opt, final_fit)
    result%fit = final_fit
    result%breakpoints = psi
    allocate(result%breakpoint_se(m), source = huge(1.0_dp))
    call hinge_matrix(z, psi, u, v)
    allocate(augmented(n, p + 2 * m))
    augmented(:, :p) = x
    augmented(:, p + 1:p + m) = u
    augmented(:, p + m + 1:) = v
    call run_lme(y, augmented, random_design, group, opt, augmented_fit)
    if (valid_lme(augmented_fit)) then
      do j = 1, m
        denom = abs(augmented_fit%beta(p + j))
        if (denom > sqrt(epsilon(1.0_dp))) then
          result%breakpoint_se(j) = sqrt(max(0.0_dp, &
              augmented_fit%beta_cov(p + m + j, p + m + j))) / denom
        end if
      end do
    end if
    result%iterations = min(iter, ctl%max_iter)
    used = min(size(histobj), max(1, result%iterations + 1))
    allocate(result%history_breakpoints(used, m), result%history_objective(used))
    result%history_breakpoints = histpsi(:used, :)
    result%history_objective = histobj(:used)
    if (valid_lme(final_fit)) then
      result%status = SEG_SUCCESS
      result%converged = iter <= ctl%max_iter
    else
      result%status = SEG_NLME_ERROR
    end if
  end subroutine fit_segmented_lme

  subroutine run_lme(y, design, random_design, group, options, fit)
    real(dp), intent(in) :: y(:), design(:,:), random_design(:,:)
    integer, intent(in) :: group(:)
    type(segmented_lme_options), intent(in) :: options
    type(lme_result), intent(out) :: fit
    call fit_lme(y, design, random_design, group, fit, random=options%random, &
        correlation=options%correlation, variance=options%variance, control=options%control)
  end subroutine run_lme

  pure logical function valid_lme(fit)
    type(lme_result), intent(in) :: fit
    valid_lme = (fit%status == NLME_SUCCESS .or. fit%status == NLME_MAX_ITER) .and. &
        allocated(fit%beta) .and. fit%log_likelihood > -huge(1.0_dp) / 10.0_dp
  end function valid_lme
end module segmented_mixed
