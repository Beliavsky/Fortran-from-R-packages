! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program fit_curve
  use rq_kinds, only: dp
  use rq_curves
  implicit none
  character(len=256) :: filename,mode,line
  real(dp),allocatable :: times(:),rates(:)
  type(curve_fit_result) :: result
  integer :: unit,ios,n,i,comma
  if(command_argument_count()<2) then
    write(*,'(a)') 'usage: fit_curve maturity_rate.csv ns|svensson'
    error stop 1
  end if
  call get_command_argument(1,filename)
  call get_command_argument(2,mode)
  open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
  if(ios/=0) error stop 'cannot open input CSV'
  n=0
  do
    read(unit,'(a)',iostat=ios) line
    if(ios/=0) exit
    comma=index(line,',')
    if(comma<=0) cycle
    block
      real(dp) :: a,b
      integer :: stat1,stat2
      read(line(:comma-1),*,iostat=stat1) a
      read(line(comma+1:),*,iostat=stat2) b
      if(stat1==0.and.stat2==0) n=n+1
    end block
  end do
  rewind(unit)
  allocate(times(n),rates(n)); i=0
  do
    read(unit,'(a)',iostat=ios) line
    if(ios/=0) exit
    comma=index(line,',')
    if(comma<=0) cycle
    block
      real(dp) :: a,b
      integer :: stat1,stat2
      read(line(:comma-1),*,iostat=stat1) a
      read(line(comma+1:),*,iostat=stat2) b
      if(stat1==0.and.stat2==0) then
        i=i+1; times(i)=a; rates(i)=b
      end if
    end block
  end do
  close(unit)
  if(i/=n.or.n<4) error stop 'insufficient numeric CSV rows'
  select case(trim(mode))
  case('ns')
    call fit_nelson_siegel(times,rates,result)
  case('svensson')
    call fit_svensson(times,rates,result)
  case default
    error stop 'mode must be ns or svensson'
  end select
  write(*,'(a,*(f14.8,1x))') 'parameters: ',result%parameters
  write(*,'(a,es14.6)') 'sse: ',result%sse
  write(*,'(a)') 'maturity,observed,fitted,residual'
  do i=1,n
    write(*,'(f10.4,3(",",f14.8))') times(i),rates(i),result%fitted(i),result%residuals(i)
  end do
end program fit_curve
