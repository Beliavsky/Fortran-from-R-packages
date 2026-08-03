! SPDX-License-Identifier: Apache-2.0
module intraday_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use intraday_kinds, only : dp
  use intraday_types, only : volume_parameters, parameter_mask, volume_model_spec, &
                             volume_fit_control, volume_model, intraday_ok, &
                             intraday_invalid_input, intraday_not_converged, &
                             intraday_numerical_failure
  use intraday_utils, only : clean_volume_data, default_parameters, parameters_valid, &
                             parameter_change, all_parameters_fixed, copy_selected_parameters
  use intraday_kalman, only : uniss_em_update
  implicit none
  private

  public :: fit_volume

contains

  subroutine fit_volume(data, model, spec, control)
    real(dp), intent(in) :: data(:, :)
    type(volume_model), intent(out) :: model
    type(volume_model_spec), intent(in), optional :: spec
    type(volume_fit_control), intent(in), optional :: control

    real(dp), allocatable :: clean(:, :), log_data(:, :)
    type(volume_model_spec) :: use_spec
    type(volume_fit_control) :: ctl
    type(volume_parameters) :: current, next, first, second, candidate
    type(volume_parameters), allocatable :: work_history(:)
    integer :: status, iter, hist_count
    character(len=160) :: message
    real(dp) :: change

    model%status = intraday_ok
    model%message = ''
    ctl = volume_fit_control()
    if (present(control)) ctl = control
    if (ctl%maxit < 1 .or. ctl%abstol <= 0.0_dp .or. .not. ieee_is_finite(ctl%abstol)) then
      model%status = intraday_invalid_input
      model%message = 'maxit must be positive and abstol must be finite and positive'
      return
    end if

    call clean_volume_data(data, clean, status=status, message=message)
    if (status /= intraday_ok) then
      model%status = status
      model%message = trim(message)
      return
    end if
    allocate(log_data(size(clean, 1), size(clean, 2)))
    log_data = log(clean)

    use_spec = volume_model_spec()
    if (present(spec)) use_spec = spec
    current = default_parameters(log_data)
    call apply_specification(current, use_spec, size(clean, 1), status, message)
    if (status /= intraday_ok) then
      model%status = status
      model%message = trim(message)
      return
    end if
    model%fixed = use_spec%is_fixed

    if (all_parameters_fixed(model%fixed)) then
      model%par = current
      model%converged = .true.
      model%iterations = 0
      model%final_change = 0.0_dp
      model%status = intraday_ok
      model%message = 'all parameters supplied as fixed values'
      if (ctl%save_history) then
        allocate(model%history(1))
        model%history(1) = current
      end if
      return
    end if

    if (ctl%save_history) then
      allocate(work_history(ctl%maxit + 1))
      work_history(1) = current
      hist_count = 1
    else
      hist_count = 0
    end if

    model%converged = .false.
    model%status = intraday_not_converged
    model%message = 'maximum iterations reached before convergence'

    do iter = 1, ctl%maxit
      if (ctl%acceleration) then
        call uniss_em_update(log_data, current, model%fixed, first, status, message)
        if (status /= intraday_ok) exit
        call uniss_em_update(log_data, first, model%fixed, second, status, message)
        if (status /= intraday_ok) exit
        call accelerate_parameters(current, first, second, model%fixed, candidate)
        if (.not. parameters_valid(candidate, size(clean, 1))) candidate = second
        change = parameter_change(first, second)
        next = candidate
      else
        call uniss_em_update(log_data, current, model%fixed, next, status, message)
        if (status /= intraday_ok) exit
        change = parameter_change(current, next)
      end if

      current = next
      if (ctl%save_history) then
        hist_count = hist_count + 1
        work_history(hist_count) = current
      end if
      if (ctl%verbose > 0) then
        if ((.not. ctl%acceleration .and. mod(iter, 25) == 0) .or. &
            (ctl%acceleration .and. mod(iter, 5) == 0)) then
          write(*, '(a,i0,a,es12.4)') 'iter: ', iter, ' change: ', change
        end if
      end if
      if (change < ctl%abstol) then
        model%converged = .true.
        model%status = intraday_ok
        model%message = 'converged'
        exit
      end if
    end do

    if (status /= intraday_ok) then
      model%status = intraday_numerical_failure
      model%message = trim(message)
    end if
    model%par = current
    model%iterations = min(iter, ctl%maxit)
    model%final_change = change

    if (ctl%save_history) then
      allocate(model%history(hist_count))
      model%history = work_history(1:hist_count)
    end if
  end subroutine fit_volume

  subroutine apply_specification(par, spec, n_bin, status, message)
    type(volume_parameters), intent(inout) :: par
    type(volume_model_spec), intent(in) :: spec
    integer, intent(in) :: n_bin
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    status = intraday_ok
    message = ''

    call copy_if_requested(par, spec%initial, spec%has_initial, n_bin, status, message)
    if (status /= intraday_ok) return
    call copy_if_requested(par, spec%fixed, spec%is_fixed, n_bin, status, message)
    if (status /= intraday_ok) return

    if (.not. parameters_valid(par, n_bin)) then
      status = intraday_invalid_input
      message = 'initial or fixed parameters are invalid'
    end if
  end subroutine apply_specification

  subroutine copy_if_requested(target, source, mask, n_bin, status, message)
    type(volume_parameters), intent(inout) :: target
    type(volume_parameters), intent(in) :: source
    type(parameter_mask), intent(in) :: mask
    integer, intent(in) :: n_bin
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    status = intraday_ok
    message = ''
    if (mask%phi) then
      if (.not. allocated(source%phi)) then
        status = intraday_invalid_input
        message = 'phi was selected but was not allocated'
        return
      end if
      if (size(source%phi) /= n_bin) then
        status = intraday_invalid_input
        message = 'phi length must equal the number of intraday bins'
        return
      end if
    end if
    call copy_selected_parameters(target, source, mask)
  end subroutine copy_if_requested

  subroutine accelerate_parameters(current, first, second, fixed, accelerated)
    type(volume_parameters), intent(in) :: current, first, second
    type(parameter_mask), intent(in) :: fixed
    type(volume_parameters), intent(out) :: accelerated
    real(dp) :: r, v, step, rnorm, vnorm
    integer :: i, j

    accelerated = current

    if (.not. fixed%phi) then
      rnorm = sqrt(sum((first%phi - current%phi)**2))
      vnorm = sqrt(sum((second%phi - 2.0_dp * first%phi + current%phi)**2))
      if (vnorm > 1.0e-8_dp) then
        step = -rnorm / vnorm
        accelerated%phi = current%phi - 2.0_dp * step * (first%phi - current%phi) + &
          step**2 * (second%phi - 2.0_dp * first%phi + current%phi)
      else
        accelerated%phi = second%phi
      end if
    end if

    if (.not. fixed%a_eta) accelerated%a_eta = scalar_accel(current%a_eta, first%a_eta, second%a_eta)
    if (.not. fixed%a_mu) accelerated%a_mu = scalar_accel(current%a_mu, first%a_mu, second%a_mu)
    if (.not. fixed%var_eta) accelerated%var_eta = scalar_accel(current%var_eta, first%var_eta, second%var_eta)
    if (.not. fixed%var_mu) accelerated%var_mu = scalar_accel(current%var_mu, first%var_mu, second%var_mu)
    if (.not. fixed%r) accelerated%r = scalar_accel(current%r, first%r, second%r)

    if (.not. fixed%x0) then
      do i = 1, 2
        accelerated%x0(i) = scalar_accel(current%x0(i), first%x0(i), second%x0(i))
      end do
    end if
    if (.not. fixed%v0) then
      do j = 1, 2
        do i = 1, 2
          accelerated%v0(i, j) = scalar_accel(current%v0(i, j), first%v0(i, j), second%v0(i, j))
        end do
      end do
      accelerated%v0 = 0.5_dp * (accelerated%v0 + transpose(accelerated%v0))
    end if

    if (accelerated%r <= 0.0_dp .or. .not. ieee_is_finite(accelerated%r)) accelerated%r = second%r
    if (accelerated%var_eta < 0.0_dp .or. .not. ieee_is_finite(accelerated%var_eta)) &
      accelerated%var_eta = second%var_eta
    if (accelerated%var_mu < 0.0_dp .or. .not. ieee_is_finite(accelerated%var_mu)) &
      accelerated%var_mu = second%var_mu

  contains

    real(dp) function scalar_accel(a, b, c) result(value)
      real(dp), intent(in) :: a, b, c
      r = b - a
      v = c - 2.0_dp * b + a
      if (abs(v) > 1.0e-8_dp) then
        step = -abs(r) / abs(v)
        value = a - step * r
      else
        value = c
      end if
      if (.not. ieee_is_finite(value)) value = c
    end function scalar_accel

  end subroutine accelerate_parameters

end module intraday_fit
