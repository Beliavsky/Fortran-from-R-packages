! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
program fit_csv
   use bayesgarch_kinds, only : dp
   use bayesgarch_rng, only : seed_rng
   use bayesgarch_sample, only : form_posterior_sample, posterior_mean, posterior_sd
   use bayesgarch_sampler, only : bayesgarch_control, bayesgarch_result, run_bayesgarch
   implicit none

   character(len=1024) :: filename
   character(len=128) :: argument
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: sample(:, :)
   real(dp) :: means(4)
   real(dp) :: sds(4)
   type(bayesgarch_control) :: control
   type(bayesgarch_result) :: result
   integer :: burn
   integer :: thin
   integer :: j

   if (command_argument_count() < 1) then
      write(*, '(a)') "usage: fit_csv FILE [N_ITER=2000] [BURN=500] [THIN=5] [N_CHAINS=1]"
      error stop 2
   end if
   call get_command_argument(1, filename)
   control = bayesgarch_control()
   control%n_iter = 2000
   control%n_chains = 1
   control%start = [0.01_dp, 0.10_dp, 0.70_dp, 20.0_dp]
   control%enforce_stationarity = .true.
   burn = 500
   thin = 5

   if (command_argument_count() >= 2) then
      call get_command_argument(2, argument)
      read(argument, *) control%n_iter
   end if
   if (command_argument_count() >= 3) then
      call get_command_argument(3, argument)
      read(argument, *) burn
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, argument)
      read(argument, *) thin
   end if
   if (command_argument_count() >= 5) then
      call get_command_argument(5, argument)
      read(argument, *) control%n_chains
   end if
   if (burn < 0 .or. burn >= control%n_iter) error stop "fit_csv: invalid burn"
   if (thin < 1) error stop "fit_csv: thin must be positive"

   call read_numeric_column(trim(filename), y)
   call seed_rng(20260723)
   call run_bayesgarch(y, result, control=control)
   call form_posterior_sample(result, burn, thin, sample)
   means = posterior_mean(sample)
   sds = posterior_sd(sample)

   write(*, '(a,i0)') "observations: ", size(y)
   write(*, '(a,i0)') "posterior draws: ", size(sample, 1)
   write(*, '(a)') "parameter       post.mean         post.sd"
   do j = 1, 4
      write(*, '(a10,2(1x,es16.8))') parameter_name(j), means(j), sds(j)
   end do
   do j = 1, control%n_chains
      write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') "chain ", j, &
         ": alpha_accept=", result%alpha_accept(j), &
         " beta_accept=", result%beta_accept(j), &
         " nu_updates=", result%nu_updates(j), &
         " constraint_reject=", result%constraint_reject(j)
   end do

contains

   subroutine read_numeric_column(path, values)
      character(len=*), intent(in) :: path
      real(dp), allocatable, intent(out) :: values(:)
      character(len=4096) :: line
      real(dp) :: value
      integer :: count
      integer :: ios
      integer :: unit
      logical :: ok

      count = 0
      open(newunit=unit, file=path, status="old", action="read", iostat=ios)
      if (ios /= 0) error stop "fit_csv: cannot open input file"
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call parse_first_number(line, value, ok)
         if (ok) count = count + 1
      end do
      close(unit)
      if (count < 2) error stop "fit_csv: fewer than two numeric observations"

      allocate(values(count))
      count = 0
      open(newunit=unit, file=path, status="old", action="read", iostat=ios)
      if (ios /= 0) error stop "fit_csv: cannot reopen input file"
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call parse_first_number(line, value, ok)
         if (ok) then
            count = count + 1
            values(count) = value
         end if
      end do
      close(unit)
   end subroutine read_numeric_column

   subroutine parse_first_number(line, value, ok)
      character(len=*), intent(in) :: line
      real(dp), intent(out) :: value
      logical, intent(out) :: ok
      character(len=len(line)) :: work
      integer :: comma
      integer :: ios

      work = adjustl(line)
      comma = index(work, ',')
      if (comma > 0) work(comma:) = ' '
      read(work, *, iostat=ios) value
      ok = ios == 0
   end subroutine parse_first_number

   pure function parameter_name(index) result(name)
      integer, intent(in) :: index
      character(len=10) :: name

      select case (index)
      case (1)
         name = "alpha0"
      case (2)
         name = "alpha1"
      case (3)
         name = "beta"
      case (4)
         name = "nu"
      case default
         name = "unknown"
      end select
   end function parameter_name

end program fit_csv
