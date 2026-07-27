! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
! This program is free software under GNU GPL version 2 only.
program fit_csv
  use pa
  implicit none
  character(len=512) :: filename, mode, aggregation, line
  integer :: unit, ios, n, i, status, nlevels(1)
  integer, allocatable :: period(:), category(:), categorical(:, :), gs(:), ge(:)
  real(dp), allocatable :: wb(:), wp(:), ret(:), value(:), growth(:), numeric(:, :), design(:, :)
  type(brinson_multi_result) :: bm
  type(attribution_summary) :: bs
  type(regression_multi_result) :: rm
  type(regression_summary) :: rs

  call get_command_argument(1, filename)
  call get_command_argument(2, mode)
  call get_command_argument(3, aggregation)
  if (len_trim(filename) == 0 .or. len_trim(mode) == 0) then
    write(*,'(a)') 'usage: fit_csv file.csv brinson|regression [arithmetic|geometric|linking]'
    error stop 1
  end if
  if (len_trim(aggregation) == 0) aggregation = 'geometric'

  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'cannot open CSV file'
  read(unit,'(a)',iostat=ios) line
  n = 0
  do
    read(unit,'(a)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) > 0) n = n + 1
  end do
  rewind(unit)
  read(unit,'(a)') line
  allocate(period(n),category(n),wb(n),wp(n),ret(n),value(n),growth(n))
  do i = 1, n
    read(unit,'(a)',iostat=ios) line
    if (ios /= 0) error stop 'unexpected end of CSV file'
    read(line,*,iostat=ios) period(i),category(i),wb(i),wp(i),ret(i),value(i),growth(i)
    if (ios /= 0) error stop 'invalid CSV row'
  end do
  close(unit)

  select case (trim(mode))
  case ('brinson')
    call fit_brinson_multi(period, category, wb, wp, ret, bm)
    if (bm%status /= 0) error stop 'Brinson fit failed'
    call summarize_brinson_multi(bm, trim(aggregation), bs)
    if (bs%status == 1) error stop 'invalid aggregation method'
    write(*,'(a)') 'allocation selection interaction active_return'
    write(*,'(4es18.9)') bs%aggregate

  case ('regression')
    allocate(categorical(n,1),numeric(n,2))
    categorical(:,1) = category
    numeric(:,1) = value
    numeric(:,2) = growth
    nlevels(1) = maxval(category)
    call build_design_matrix(categorical, nlevels, numeric, design, gs, ge, status)
    if (status /= 0) error stop 'design-matrix construction failed'
    call fit_regression_multi(period, ret, design, wb, wp, rm)
    if (rm%status /= 0) error stop 'regression fit failed'
    call summarize_regression_multi(rm, gs, ge, trim(aggregation), rs)
    if (rs%status == 1) error stop 'invalid aggregation method'
    write(*,'(a)') 'sector value growth residual [portfolio benchmark] active_return'
    write(*,'(*(es18.9,1x))') rs%aggregate

  case default
    error stop 'mode must be brinson or regression'
  end select
end program fit_csv
