! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_data
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use strand_kinds, only : dp
  use strand_linalg, only : correlation
  implicit none
  private

  type, public :: cross_section_state
    integer, allocatable :: id(:)
    real(dp), allocatable :: values(:, :)
    logical, allocatable :: carry_forward(:)
    integer, allocatable :: data_date(:)
    integer :: current_date = 0
  end type cross_section_state

  type, public :: cross_section_stats
    integer :: input_rows = 0
    integer :: carry_forward_rows = 0
    integer, allocatable :: missing(:)
    real(dp), allocatable :: prior_correlation(:)
  end type cross_section_stats

  type, public :: portfolio_state
    integer, allocatable :: internal_shares(:, :)
    integer, allocatable :: external_shares(:, :)
  end type portfolio_state

  public :: update_cross_section, compute_cross_section_stats
  public :: initialize_portfolio, consolidated_shares, apply_adjustment_ratio

contains

  function update_cross_section(previous, new_id, new_values, date, carry_values, &
    carry_columns, replace_values, replace_columns) result(current)
    type(cross_section_state), intent(in), optional :: previous
    integer, intent(in) :: new_id(:), date
    real(dp), intent(in) :: new_values(:, :)
    real(dp), intent(in), optional :: carry_values(:), replace_values(:)
    integer, intent(in), optional :: carry_columns(:), replace_columns(:)
    type(cross_section_state) :: current
    integer, allocatable :: ids(:), source(:)
    integer :: nnew, nold, ncol, nout, i, j, idx

    nnew = size(new_id)
    ncol = size(new_values, 2)
    if (size(new_values, 1) /= nnew) error stop 'update_cross_section: row mismatch'
    if (has_duplicates(new_id)) error stop 'update_cross_section: duplicate ids'
    nold = 0
    if (present(previous)) then
      if (allocated(previous%id)) then
        nold = size(previous%id)
        if (size(previous%values, 2) /= ncol) error stop 'update_cross_section: column mismatch'
      end if
    end if
    allocate(ids(nnew + nold), source(nnew + nold))
    nout = nnew
    ids(1:nnew) = new_id
    source(1:nnew) = 0
    if (present(previous) .and. nold > 0) then
      do i = 1, nold
        if (.not. any(new_id == previous%id(i))) then
          nout = nout + 1
          ids(nout) = previous%id(i)
          source(nout) = i
        end if
      end do
    end if
    allocate(current%id(nout), current%values(nout, ncol), current%carry_forward(nout), current%data_date(nout))
    current%id = ids(1:nout)
    current%values = 0.0_dp
    current%carry_forward = .false.
    current%data_date = date
    current%current_date = date
    current%values(1:nnew, :) = new_values
    if (present(previous)) then
      do i = nnew + 1, nout
        idx = source(i)
        current%values(i, :) = previous%values(idx, :)
        current%carry_forward(i) = .true.
        current%data_date(i) = previous%data_date(idx)
      end do
    end if
    if (present(carry_columns) .and. present(carry_values)) then
      if (size(carry_columns) /= size(carry_values)) error stop 'update_cross_section: carry defaults mismatch'
      do j = 1, size(carry_columns)
        if (carry_columns(j) < 1 .or. carry_columns(j) > ncol) error stop 'update_cross_section: bad carry column'
        do i = 1, nout
          if (current%carry_forward(i)) current%values(i, carry_columns(j)) = carry_values(j)
        end do
      end do
    end if
    if (present(replace_columns) .and. present(replace_values)) then
      if (size(replace_columns) /= size(replace_values)) error stop 'update_cross_section: replacements mismatch'
      do j = 1, size(replace_columns)
        if (replace_columns(j) < 1 .or. replace_columns(j) > ncol) error stop 'update_cross_section: bad replace column'
        do i = 1, nout
          if (.not. ieee_is_finite(current%values(i, replace_columns(j)))) then
            current%values(i, replace_columns(j)) = replace_values(j)
          end if
        end do
      end do
    end if
  end function update_cross_section

  function compute_cross_section_stats(current, previous) result(stats)
    type(cross_section_state), intent(in) :: current
    type(cross_section_state), intent(in), optional :: previous
    type(cross_section_stats) :: stats
    real(dp), allocatable :: x(:), y(:)
    integer :: ncol, i, j, k, matched

    if (.not. allocated(current%values)) return
    ncol = size(current%values, 2)
    stats%input_rows = size(current%values, 1)
    stats%carry_forward_rows = count(current%carry_forward)
    allocate(stats%missing(ncol), stats%prior_correlation(ncol))
    do j = 1, ncol
      stats%missing(j) = count(.not. ieee_is_finite(current%values(:, j)))
      stats%prior_correlation(j) = 0.0_dp
    end do
    if (.not. present(previous)) return
    if (.not. allocated(previous%values)) return
    allocate(x(min(size(current%id), size(previous%id))), y(min(size(current%id), size(previous%id))))
    do j = 1, ncol
      matched = 0
      do i = 1, size(current%id)
        do k = 1, size(previous%id)
          if (current%id(i) == previous%id(k)) then
            if (ieee_is_finite(current%values(i, j)) .and. ieee_is_finite(previous%values(k, j))) then
              matched = matched + 1
              x(matched) = current%values(i, j)
              y(matched) = previous%values(k, j)
            end if
            exit
          end if
        end do
      end do
      if (matched >= 2) stats%prior_correlation(j) = correlation(x(1:matched), y(1:matched))
    end do
  end function compute_cross_section_stats

  subroutine initialize_portfolio(portfolio, nsecurity, nstrategy, initial_shares)
    type(portfolio_state), intent(out) :: portfolio
    integer, intent(in) :: nsecurity, nstrategy
    integer, intent(in), optional :: initial_shares(:, :)

    allocate(portfolio%internal_shares(nsecurity, nstrategy), portfolio%external_shares(nsecurity, nstrategy))
    portfolio%internal_shares = 0
    portfolio%external_shares = 0
    if (present(initial_shares)) then
      if (any(shape(initial_shares) /= [nsecurity, nstrategy])) error stop 'initialize_portfolio: dimension mismatch'
      portfolio%external_shares = initial_shares
    end if
  end subroutine initialize_portfolio

  function consolidated_shares(portfolio) result(shares)
    type(portfolio_state), intent(in) :: portfolio
    integer, allocatable :: shares(:, :)
    allocate(shares(size(portfolio%internal_shares, 1), size(portfolio%internal_shares, 2)))
    shares = portfolio%internal_shares + portfolio%external_shares
  end function consolidated_shares

  subroutine apply_adjustment_ratio(portfolio, ratio)
    type(portfolio_state), intent(inout) :: portfolio
    real(dp), intent(in) :: ratio(:)
    integer :: i, j

    if (size(ratio) /= size(portfolio%internal_shares, 1) .or. any(ratio <= 0.0_dp)) then
      error stop 'apply_adjustment_ratio: invalid ratio'
    end if
    do j = 1, size(portfolio%internal_shares, 2)
      do i = 1, size(ratio)
        portfolio%internal_shares(i, j) = nint(real(portfolio%internal_shares(i, j), dp) / ratio(i))
        portfolio%external_shares(i, j) = nint(real(portfolio%external_shares(i, j), dp) / ratio(i))
      end do
    end do
  end subroutine apply_adjustment_ratio

  logical function has_duplicates(x) result(duplicates)
    integer, intent(in) :: x(:)
    integer :: i
    duplicates = .false.
    do i = 2, size(x)
      if (any(x(1:i - 1) == x(i))) then
        duplicates = .true.
        return
      end if
    end do
  end function has_duplicates

end module strand_data
