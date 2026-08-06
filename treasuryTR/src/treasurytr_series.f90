! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

module treasurytr_series
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
  use treasurytr_kinds, only : dp
  use treasurytr_math, only : period_total_return
  implicit none
  private

  integer, parameter, public :: tt_success = 0
  integer, parameter, public :: tt_err_size = 1
  integer, parameter, public :: tt_err_maturity = 2
  integer, parameter, public :: tt_err_scale = 3

  public :: carry_forward
  public :: percent_to_decimal
  public :: prepare_yields
  public :: total_return

  interface carry_forward
    module procedure carry_forward_vector
    module procedure carry_forward_matrix
  end interface carry_forward

  interface prepare_yields
    module procedure prepare_yields_vector
    module procedure prepare_yields_matrix
  end interface prepare_yields

  interface total_return
    module procedure total_return_vector
    module procedure total_return_matrix
  end interface total_return

contains

  pure elemental function percent_to_decimal(value) result(out)
    real(dp), intent(in) :: value
    real(dp) :: out

    out = value / 100.0_dp
  end function percent_to_decimal

  pure function carry_forward_vector(x) result(out)
    real(dp), intent(in) :: x(:)
    real(dp) :: out(size(x))
    integer :: i

    out = x
    do i = 2, size(out)
      if (ieee_is_nan(out(i)) .and. .not. ieee_is_nan(out(i - 1))) out(i) = out(i - 1)
    end do
  end function carry_forward_vector

  pure function carry_forward_matrix(x) result(out)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: out(size(x, 1), size(x, 2))
    integer :: j

    do j = 1, size(x, 2)
      out(:, j) = carry_forward_vector(x(:, j))
    end do
  end function carry_forward_matrix

  pure function prepare_yields_vector(x, na_locf, percent_adjust) result(out)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: na_locf
    logical, intent(in), optional :: percent_adjust
    real(dp) :: out(size(x))
    logical :: do_locf, do_percent

    do_locf = .true.
    do_percent = .true.
    if (present(na_locf)) do_locf = na_locf
    if (present(percent_adjust)) do_percent = percent_adjust

    out = x
    if (do_locf) out = carry_forward(out)
    if (do_percent) out = percent_to_decimal(out)
  end function prepare_yields_vector

  pure function prepare_yields_matrix(x, na_locf, percent_adjust) result(out)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in), optional :: na_locf
    logical, intent(in), optional :: percent_adjust
    real(dp) :: out(size(x, 1), size(x, 2))
    logical :: do_locf, do_percent

    do_locf = .true.
    do_percent = .true.
    if (present(na_locf)) do_locf = na_locf
    if (present(percent_adjust)) do_percent = percent_adjust

    out = x
    if (do_locf) out = carry_forward(out)
    if (do_percent) out = percent_to_decimal(out)
  end function prepare_yields_matrix

  function total_return_vector(yields, maturity, scale, mdur, convex, &
      source_compatible, status, message) result(out)
    real(dp), intent(in) :: yields(:)
    real(dp), intent(in) :: maturity
    real(dp), intent(in), optional :: scale
    real(dp), intent(in), optional :: mdur(:)
    real(dp), intent(in), optional :: convex(:)
    logical, intent(in), optional :: source_compatible
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    real(dp) :: out(size(yields))
    real(dp) :: periods_per_year
    logical :: use_source
    integer :: i

    call initialize_status(status, message)
    out = quiet_nan()

    if (maturity <= 0.0_dp) then
      call fail_status(tt_err_maturity, 'maturity must be positive', status, message)
      return
    end if

    periods_per_year = 261.0_dp
    if (present(scale)) periods_per_year = scale
    if (periods_per_year <= 0.0_dp) then
      call fail_status(tt_err_scale, 'scale must be positive', status, message)
      return
    end if

    if (present(mdur)) then
      if (size(mdur) /= size(yields)) then
        call fail_status(tt_err_size, 'mdur must have the same size as yields', status, message)
        return
      end if
    end if
    if (present(convex)) then
      if (size(convex) /= size(yields)) then
        call fail_status(tt_err_size, 'convex must have the same size as yields', status, message)
        return
      end if
    end if

    use_source = .true.
    if (present(source_compatible)) use_source = source_compatible
    if (size(yields) == 0) return

    do i = 2, size(yields)
      if (present(mdur) .and. present(convex)) then
        out(i) = period_total_return(yields(i), yields(i - 1), maturity, periods_per_year, &
          mdur(i), convex(i), use_source)
      else if (present(mdur)) then
        out(i) = period_total_return(yields(i), yields(i - 1), maturity, periods_per_year, &
          mdur_current=mdur(i), source_compatible=use_source)
      else if (present(convex)) then
        out(i) = period_total_return(yields(i), yields(i - 1), maturity, periods_per_year, &
          convex_current=convex(i), source_compatible=use_source)
      else
        out(i) = period_total_return(yields(i), yields(i - 1), maturity, periods_per_year, &
          source_compatible=use_source)
      end if
    end do
  end function total_return_vector

  function total_return_matrix(yields, maturity, scale, mdur, convex, &
      source_compatible, status, message) result(out)
    real(dp), intent(in) :: yields(:, :)
    real(dp), intent(in) :: maturity
    real(dp), intent(in), optional :: scale
    real(dp), intent(in), optional :: mdur(:, :)
    real(dp), intent(in), optional :: convex(:, :)
    logical, intent(in), optional :: source_compatible
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    real(dp) :: out(size(yields, 1), size(yields, 2))
    integer :: j, local_status
    character(len=160) :: local_message

    call initialize_status(status, message)
    out = quiet_nan()

    if (present(mdur)) then
      if (any(shape(mdur) /= shape(yields))) then
        call fail_status(tt_err_size, 'mdur must have the same shape as yields', status, message)
        return
      end if
    end if
    if (present(convex)) then
      if (any(shape(convex) /= shape(yields))) then
        call fail_status(tt_err_size, 'convex must have the same shape as yields', status, message)
        return
      end if
    end if

    do j = 1, size(yields, 2)
      if (present(mdur) .and. present(convex)) then
        out(:, j) = total_return_vector(yields(:, j), maturity, scale, mdur(:, j), &
          convex(:, j), source_compatible, local_status, local_message)
      else if (present(mdur)) then
        out(:, j) = total_return_vector(yields(:, j), maturity, scale, mdur=mdur(:, j), &
          source_compatible=source_compatible, status=local_status, message=local_message)
      else if (present(convex)) then
        out(:, j) = total_return_vector(yields(:, j), maturity, scale, convex=convex(:, j), &
          source_compatible=source_compatible, status=local_status, message=local_message)
      else
        out(:, j) = total_return_vector(yields(:, j), maturity, scale, &
          source_compatible=source_compatible, status=local_status, message=local_message)
      end if
      if (local_status /= tt_success) then
        call fail_status(local_status, trim(local_message), status, message)
        return
      end if
    end do
  end function total_return_matrix

  subroutine initialize_status(status, message)
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    if (present(status)) status = tt_success
    if (present(message)) message = ''
  end subroutine initialize_status

  subroutine fail_status(code, text, status, message)
    integer, intent(in) :: code
    character(len=*), intent(in) :: text
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    if (present(status)) status = code
    if (present(message)) message = text
  end subroutine fail_status

  pure function quiet_nan() result(value)
    real(dp) :: value

    value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan

end module treasurytr_series
