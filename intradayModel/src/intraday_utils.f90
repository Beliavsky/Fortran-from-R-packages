! SPDX-License-Identifier: Apache-2.0
module intraday_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use intraday_kinds, only : dp
  use intraday_types, only : volume_parameters, parameter_mask, volume_model_spec, &
                             error_metrics, intraday_ok, intraday_invalid_input
  implicit none
  private

  public :: initialize_volume_spec, clean_volume_data, default_parameters
  public :: parameters_valid, parameter_change, compute_error_metrics
  public :: all_parameters_fixed, copy_selected_parameters

contains

  subroutine initialize_volume_spec(spec, n_bin)
    type(volume_model_spec), intent(out) :: spec
    integer, intent(in) :: n_bin

    if (allocated(spec%fixed%phi)) deallocate(spec%fixed%phi)
    if (allocated(spec%initial%phi)) deallocate(spec%initial%phi)
    allocate(spec%fixed%phi(max(0, n_bin)), spec%initial%phi(max(0, n_bin)))
    spec%fixed%phi = 0.0_dp
    spec%initial%phi = 0.0_dp
    spec%is_fixed = parameter_mask()
    spec%has_initial = parameter_mask()
  end subroutine initialize_volume_spec

  subroutine clean_volume_data(data, cleaned, kept_columns, status, message)
    real(dp), intent(in) :: data(:, :)
    real(dp), allocatable, intent(out) :: cleaned(:, :)
    integer, allocatable, intent(out), optional :: kept_columns(:)
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    logical, allocatable :: keep(:)
    integer :: j, k, n_keep

    if (present(status)) status = intraday_ok
    if (present(message)) message = ''

    if (size(data, 1) < 1 .or. size(data, 2) < 1) then
      allocate(cleaned(0, 0))
      if (present(kept_columns)) allocate(kept_columns(0))
      if (present(status)) status = intraday_invalid_input
      if (present(message)) message = 'data must have at least one bin and one day'
      return
    end if

    allocate(keep(size(data, 2)))
    do j = 1, size(data, 2)
      keep(j) = all(ieee_is_finite(data(:, j)))
      if (keep(j)) keep(j) = all(data(:, j) > 0.0_dp)
    end do
    n_keep = count(keep)
    if (n_keep < 1) then
      allocate(cleaned(size(data, 1), 0))
      if (present(kept_columns)) allocate(kept_columns(0))
      if (present(status)) status = intraday_invalid_input
      if (present(message)) message = 'no complete positive finite trading day remains'
      return
    end if

    allocate(cleaned(size(data, 1), n_keep))
    if (present(kept_columns)) allocate(kept_columns(n_keep))
    k = 0
    do j = 1, size(data, 2)
      if (keep(j)) then
        k = k + 1
        cleaned(:, k) = data(:, j)
        if (present(kept_columns)) kept_columns(k) = j
      end if
    end do
  end subroutine clean_volume_data

  function default_parameters(log_data) result(par)
    real(dp), intent(in) :: log_data(:, :)
    type(volume_parameters) :: par
    integer :: i
    real(dp) :: overall

    overall = sum(log_data) / real(size(log_data), dp)
    par%x0 = [overall, 0.0_dp]
    allocate(par%phi(size(log_data, 1)))
    do i = 1, size(log_data, 1)
      par%phi(i) = sum(log_data(i, :)) / real(size(log_data, 2), dp) - overall
    end do
  end function default_parameters

  pure logical function parameters_valid(par, n_bin) result(ok)
    type(volume_parameters), intent(in) :: par
    integer, intent(in) :: n_bin
    real(dp) :: detv

    ok = allocated(par%phi)
    if (.not. ok) return
    ok = size(par%phi) == n_bin
    if (.not. ok) return
    ok = ieee_is_finite(par%a_eta) .and. ieee_is_finite(par%a_mu) .and. &
         ieee_is_finite(par%var_eta) .and. ieee_is_finite(par%var_mu) .and. &
         ieee_is_finite(par%r) .and. all(ieee_is_finite(par%x0)) .and. &
         all(ieee_is_finite(par%v0)) .and. all(ieee_is_finite(par%phi))
    if (.not. ok) return
    ok = par%var_eta >= 0.0_dp .and. par%var_mu >= 0.0_dp .and. par%r > 0.0_dp
    if (.not. ok) return
    ok = abs(par%v0(1, 2) - par%v0(2, 1)) <= 1.0e-10_dp * &
         max(1.0_dp, maxval(abs(par%v0)))
    if (.not. ok) return
    detv = par%v0(1, 1) * par%v0(2, 2) - par%v0(1, 2) * par%v0(2, 1)
    ok = par%v0(1, 1) >= 0.0_dp .and. par%v0(2, 2) >= 0.0_dp .and. detv >= -1.0e-14_dp
  end function parameters_valid

  real(dp) function parameter_change(a, b) result(value)
    type(volume_parameters), intent(in) :: a, b
    real(dp) :: ss

    ss = (a%a_eta - b%a_eta)**2 + (a%a_mu - b%a_mu)**2 + &
         (a%var_eta - b%var_eta)**2 + (a%var_mu - b%var_mu)**2 + &
         (a%r - b%r)**2 + sum((a%x0 - b%x0)**2) + sum((a%v0 - b%v0)**2)
    if (allocated(a%phi) .and. allocated(b%phi)) then
      if (size(a%phi) == size(b%phi)) ss = ss + sum((a%phi - b%phi)**2)
    end if
    value = sqrt(max(0.0_dp, ss))
  end function parameter_change

  function compute_error_metrics(reference, predicted) result(metrics)
    real(dp), intent(in) :: reference(:), predicted(:)
    type(error_metrics) :: metrics
    logical, allocatable :: nonzero(:)
    integer :: n

    n = min(size(reference), size(predicted))
    if (n < 1) return
    metrics%mae = sum(abs(predicted(1:n) - reference(1:n))) / real(n, dp)
    metrics%rmse = sqrt(sum((predicted(1:n) - reference(1:n))**2) / real(n, dp))
    allocate(nonzero(n))
    nonzero = abs(reference(1:n)) > tiny(1.0_dp)
    if (any(nonzero)) then
      metrics%mape = sum(abs((predicted(1:n) - reference(1:n)) / &
                         merge(reference(1:n), 1.0_dp, nonzero)), mask=nonzero) / &
                         real(count(nonzero), dp)
    else
      metrics%mape = huge(1.0_dp)
    end if
  end function compute_error_metrics

  logical function all_parameters_fixed(mask) result(answer)
    type(parameter_mask), intent(in) :: mask
    answer = mask%a_eta .and. mask%a_mu .and. mask%var_eta .and. mask%var_mu .and. &
             mask%r .and. mask%phi .and. mask%x0 .and. mask%v0
  end function all_parameters_fixed

  subroutine copy_selected_parameters(target, source, mask)
    type(volume_parameters), intent(inout) :: target
    type(volume_parameters), intent(in) :: source
    type(parameter_mask), intent(in) :: mask

    if (mask%a_eta) target%a_eta = source%a_eta
    if (mask%a_mu) target%a_mu = source%a_mu
    if (mask%var_eta) target%var_eta = source%var_eta
    if (mask%var_mu) target%var_mu = source%var_mu
    if (mask%r) target%r = source%r
    if (mask%x0) target%x0 = source%x0
    if (mask%v0) target%v0 = source%v0
    if (mask%phi) then
      if (allocated(target%phi)) deallocate(target%phi)
      allocate(target%phi(size(source%phi)))
      target%phi = source%phi
    end if
  end subroutine copy_selected_parameters

end module intraday_utils
