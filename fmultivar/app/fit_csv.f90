! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program fit_csv
  use fmultivar, only : dp, skew_fit_result, fit_multivariate_normal, &
    fit_skew_normal, fit_skew_t, fit_skew_cauchy
  implicit none
  character(len=1024) :: filename,method,arg
  real(dp),allocatable :: x(:,:)
  type(skew_fit_result) :: fit
  real(dp) :: fixed_nu
  integer :: max_iter,i,j
  logical :: use_fixed

  if(command_argument_count()<2)then
    write(*,'(a)')'Usage: fit_csv FILE METHOD [FIXED_NU] [MAX_ITER]'
    write(*,'(a)')'METHOD: normal, snorm, st, or cauchy. Use FIXED_NU <= 0 for free st nu.'
    error stop 2
  end if
  call get_command_argument(1,filename);call get_command_argument(2,method)
  fixed_nu=-1.0_dp;max_iter=1800
  if(command_argument_count()>=3)then
    call get_command_argument(3,arg);read(arg,*)fixed_nu
  end if
  if(command_argument_count()>=4)then
    call get_command_argument(4,arg);read(arg,*)max_iter
  end if
  call read_csv_numeric(trim(filename),x)
  select case(trim(adjustl(method)))
  case('normal','norm')
    fit=fit_multivariate_normal(x)
  case('snorm','skew-normal')
    fit=fit_skew_normal(x,max_iter,3.0e-6_dp)
  case('st','skew-t')
    use_fixed=fixed_nu>0.0_dp
    if(use_fixed)then
      fit=fit_skew_t(x,fixed_nu,max_iter,3.0e-6_dp)
    else
      fit=fit_skew_t(x,max_iter=max_iter,tol=3.0e-6_dp)
    end if
  case('cauchy','skew-cauchy')
    fit=fit_skew_cauchy(x,max_iter,5.0e-6_dp)
  case default
    error stop 'Unknown fitting method'
  end select
  write(*,'(a,l1)')'converged: ',fit%converged
  write(*,'(a,i0)')'iterations: ',fit%iterations
  write(*,'(a,f18.8)')'log_likelihood: ',fit%loglik
  write(*,'(a)',advance='no')'location:'
  do i=1,size(fit%location);write(*,'(1x,es14.6)',advance='no')fit%location(i);end do
  write(*,*)
  write(*,'(a)',advance='no')'alpha:'
  do i=1,size(fit%alpha);write(*,'(1x,es14.6)',advance='no')fit%alpha(i);end do
  write(*,*)
  write(*,'(a,es14.6)')'nu: ',fit%nu
  write(*,'(a)')'omega:'
  do i=1,size(fit%omega,1)
    write(*,'(*(1x,es14.6))')(fit%omega(i,j),j=1,size(fit%omega,2))
  end do
contains
  subroutine read_csv_numeric(path,data)
    character(len=*),intent(in)::path
    real(dp),allocatable,intent(out)::data(:,:)
    character(len=8192)::line,date_token
    integer::unit,ios,nrow,ncol,i,j
    real(dp),allocatable::row(:)
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'Could not open CSV file'
    read(unit,'(a)',iostat=ios)line
    if(ios/=0)error stop 'CSV header is missing'
    ncol=count([(line(i:i)==',',i=1,len_trim(line))])
    if(ncol<1)error stop 'CSV must contain Date plus at least one numeric column'
    nrow=0
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)>0)nrow=nrow+1
    end do
    rewind(unit);read(unit,'(a)')line
    allocate(data(nrow,ncol),row(ncol));i=0
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)==0)cycle
      call commas_to_spaces(line)
      read(line,*,iostat=ios)date_token,(row(j),j=1,ncol)
      if(ios/=0)error stop 'Malformed CSV data row'
      i=i+1;data(i,:)=row
    end do
    close(unit)
  end subroutine read_csv_numeric

  subroutine commas_to_spaces(line)
    character(len=*),intent(inout)::line
    integer::i
    do i=1,len_trim(line);if(line(i:i)==',')line(i:i)=' ';end do
  end subroutine commas_to_spaces
end program fit_csv
