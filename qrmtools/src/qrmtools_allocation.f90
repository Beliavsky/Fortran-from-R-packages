! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_allocation
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : allocation_result
  use qrmtools_stats, only : variance_value
  use qrmtools_risk, only : var_np
  implicit none
  private

  public :: alloc_ellip, conditioning_var, alloc_np

contains

  function alloc_ellip(total, location, scale) result(allocation)
    real(dp), intent(in) :: total
    real(dp), intent(in) :: location(:)
    real(dp), intent(in) :: scale(:,:)
    real(dp), allocatable :: allocation(:)
    real(dp), allocatable :: rows(:)

    rows = sum(scale, dim=2)
    allocate(allocation(size(location)))
    allocation = location + rows / sum(rows) * total
  end function alloc_ellip

  function conditioning_var(x, lower_level, upper_level) result(conditional)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in) :: lower_level
    real(dp), intent(in), optional :: upper_level
    real(dp), allocatable :: conditional(:,:)
    real(dp), allocatable :: sums(:)
    real(dp) :: upper
    real(dp) :: low_value
    real(dp) :: high_value
    integer :: i
    integer :: n

    upper = 1.0_dp
    if (present(upper_level)) upper = upper_level
    sums = sum(x, dim=2)
    low_value = var_np(sums, lower_level)
    high_value = var_np(sums, upper)
    n = count(sums > low_value .and. sums <= high_value)
    allocate(conditional(n, size(x,2)))

    n = 0
    do i = 1, size(x,1)
      if (sums(i) > low_value .and. sums(i) <= high_value) then
        n = n + 1
        conditional(n,:) = x(i,:)
      end if
    end do
  end function conditioning_var

  function alloc_np(x, lower_level, upper_level, include_conditional) result(output)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in) :: lower_level
    real(dp), intent(in), optional :: upper_level
    logical, intent(in), optional :: include_conditional
    type(allocation_result) :: output
    real(dp), allocatable :: conditional(:,:)
    real(dp), allocatable :: sums(:)
    real(dp) :: upper
    real(dp) :: low_value
    real(dp) :: high_value
    logical :: include
    integer :: i
    integer :: j
    integer :: n

    upper = 1.0_dp
    if (present(upper_level)) upper = upper_level
    sums = sum(x, dim=2)
    low_value = var_np(sums, lower_level)
    high_value = var_np(sums, upper)
    n = count(sums > low_value .and. sums <= high_value)
    allocate(conditional(n, size(x,2)))

    n = 0
    do i = 1, size(x,1)
      if (sums(i) > low_value .and. sums(i) <= high_value) then
        n = n + 1
        conditional(n,:) = x(i,:)
      end if
    end do

    output%n = size(conditional,1)
    if (output%n == 0) then
      output%message = 'The conditioning set is empty.'
      return
    end if

    allocate(output%allocation(size(x,2)))
    allocate(output%standard_error(size(x,2)))
    output%allocation = sum(conditional, dim=1) / real(output%n, dp)
    do j = 1, size(x,2)
      output%standard_error(j) = sqrt(variance_value(conditional(:,j)) / real(output%n, dp))
    end do

    include = .false.
    if (present(include_conditional)) include = include_conditional
    if (include) output%conditional = conditional
    output%ok = .true.
  end function alloc_np

end module qrmtools_allocation
