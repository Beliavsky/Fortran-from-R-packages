! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

program fit_csv
  use betategarch, only : dp, tegarch_fit_result, tegarch_fit, tegarch_bic_per_observation
  implicit none

  character(len=1024) :: filename, arg
  real(dp), allocatable :: y(:)
  type(tegarch_fit_result) :: fit
  integer :: column, components, asym_int, skew_int, i
  logical :: asym, skewed

  if (command_argument_count() < 1) then
    print '(a)', "usage: fit_csv FILE [COLUMN=1] [COMPONENTS=1] [ASYM=1] [SKEW=1]"
    stop 1
  end if

  call get_command_argument(1, filename)
  column = 1
  components = 1
  asym_int = 1
  skew_int = 1
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) column
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) components
  end if
  if (command_argument_count() >= 4) then
    call get_command_argument(4, arg)
    read(arg, *) asym_int
  end if
  if (command_argument_count() >= 5) then
    call get_command_argument(5, arg)
    read(arg, *) skew_int
  end if
  asym = asym_int /= 0
  skewed = skew_int /= 0

  call read_csv_column(trim(filename), column, y)
  if (size(y) < 2) error stop "fit_csv: fewer than two numeric observations"
  call tegarch_fit(y, components, asym, skewed, fit, compute_hessian=.true., &
    max_iterations=2500, tolerance=1.0e-7_dp)

  print '(a,1x,i0)', "observations:", size(y)
  print '(a,1x,f18.8)', "log-likelihood:", fit%log_likelihood
  print '(a,1x,f18.8)', "BIC per observation:", &
    tegarch_bic_per_observation(fit%log_likelihood, size(y), size(fit%free_parameters))
  print '(a,1x,i0)', "convergence code:", fit%convergence
  print '(a,1x,i0)', "iterations:", fit%iterations
  print '(a,1x,i0)', "evaluations:", fit%evaluations
  print '(a)', "free parameter estimates:"
  do i = 1, size(fit%free_parameters)
    print '(i3,1x,es20.12)', i, fit%free_parameters(i)
  end do
  print '(a,1x,l1)', "Hessian available:", fit%hessian_available
  print '(a,1x,l1)', "covariance available:", fit%covariance_available

contains

  subroutine read_csv_column(path, selected_column, values)
    character(len=*), intent(in) :: path
    integer, intent(in) :: selected_column
    real(dp), allocatable, intent(out) :: values(:)

    character(len=4096) :: line
    real(dp) :: value
    integer :: unit, ios, n, j
    logical :: ok

    if (selected_column < 1) error stop "read_csv_column: column must be positive"
    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "read_csv_column: cannot open input file"
    n = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      call parse_csv_value(line, selected_column, value, ok)
      if (ok) n = n + 1
    end do
    rewind(unit)
    allocate(values(n))
    j = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      call parse_csv_value(line, selected_column, value, ok)
      if (ok) then
        j = j + 1
        values(j) = value
      end if
    end do
    close(unit)
  end subroutine read_csv_column

  subroutine parse_csv_value(line, selected_column, value, ok)
    character(len=*), intent(in) :: line
    integer, intent(in) :: selected_column
    real(dp), intent(out) :: value
    logical, intent(out) :: ok

    character(len=:), allocatable :: token
    integer :: start, comma, current_column, ios, line_length

    ok = .false.
    value = 0.0_dp
    line_length = len_trim(line)
    if (line_length == 0) return
    if (line(1:1) == "#") return

    start = 1
    current_column = 1
    do
      comma = index(line(start:line_length), ",")
      if (comma == 0) then
        token = adjustl(line(start:line_length))
      else
        token = adjustl(line(start:start + comma - 2))
      end if
      if (current_column == selected_column) then
        token = trim(token)
        if (len(token) >= 2) then
          if (token(1:1) == '"' .and. token(len(token):len(token)) == '"') then
            token = token(2:len(token)-1)
          end if
        end if
        read(token, *, iostat=ios) value
        ok = ios == 0
        return
      end if
      if (comma == 0) return
      start = start + comma
      current_column = current_column + 1
      if (start > line_length) return
    end do
  end subroutine parse_csv_value

end program fit_csv
