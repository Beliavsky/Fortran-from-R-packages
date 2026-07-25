! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program fit_csv
  use fcopulae, only : dp, copula_fit_result, archm_fit, elliptical_fit, ev_fit
  implicit none
  character(len=1024)::filename,family,type_name,arg
  real(dp),allocatable::u(:),v(:)
  type(copula_fit_result)::fit
  integer::type_id,max_iter,i

  if(command_argument_count()<3)then
    write(*,'(a)')'Usage: fit_csv FILE FAMILY TYPE [MAX_ITER]'
    write(*,'(a)')'FAMILY: archm, elliptical, or ev'
    write(*,'(a)')'TYPE: 1..22; norm/cauchy/t/logistic/laplace/kotz/epower; or gumbel/galambos/husler.reiss/tawn/bb5'
    error stop 2
  end if
  call get_command_argument(1,filename)
  call get_command_argument(2,family)
  call get_command_argument(3,type_name)
  max_iter=1200
  if(command_argument_count()>=4)then;call get_command_argument(4,arg);read(arg,*)max_iter;end if
  call read_uv_csv(trim(filename),u,v)
  select case(trim(adjustl(family)))
  case('archm','archimedean')
    read(type_name,*)type_id
    fit=archm_fit(u,v,type_id,max_iter=max_iter)
  case('elliptical','ell')
    fit=elliptical_fit(u,v,trim(type_name),max_iter=max_iter)
  case('ev','extreme-value','extreme_value')
    fit=ev_fit(u,v,trim(type_name),max_iter=max_iter)
  case default
    error stop 'Unknown copula family'
  end select
  write(*,'(a,l1)')'converged: ',fit%converged
  write(*,'(a,i0)')'iterations: ',fit%iterations
  write(*,'(a,i0)')'evaluations: ',fit%evaluations
  write(*,'(a,es18.9)')'log_likelihood: ',fit%loglik
  write(*,'(a,es18.9)')'AIC: ',fit%aic
  write(*,'(a,es18.9)')'BIC: ',fit%bic
  write(*,'(a)',advance='no')'parameters:'
  do i=1,size(fit%param);write(*,'(1x,es18.9)',advance='no')fit%param(i);end do
  write(*,*)
contains
  subroutine read_uv_csv(path,x,y)
    character(len=*),intent(in)::path
    real(dp),allocatable,intent(out)::x(:),y(:)
    character(len=8192)::line,token
    integer::unit,ios,nrow,i,ncomma
    real(dp)::a,b
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'Could not open CSV file'
    read(unit,'(a)',iostat=ios)line
    if(ios/=0)error stop 'CSV header is missing'
    ncomma=count([(line(i:i)==',',i=1,len_trim(line))])
    if(ncomma<1.or.ncomma>2)error stop 'CSV must have U,V or Date,U,V columns'
    nrow=0
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)>0)nrow=nrow+1
    end do
    rewind(unit);read(unit,'(a)')line
    allocate(x(nrow),y(nrow));i=0
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(len_trim(line)==0)cycle
      call commas_to_spaces(line)
      if(ncomma==1)then
        read(line,*,iostat=ios)a,b
      else
        read(line,*,iostat=ios)token,a,b
      end if
      if(ios/=0)error stop 'Malformed CSV row'
      if(a<=0.0_dp.or.a>=1.0_dp.or.b<=0.0_dp.or.b>=1.0_dp)error stop 'U and V must be inside (0,1)'
      i=i+1;x(i)=a;y(i)=b
    end do
    close(unit)
  end subroutine read_uv_csv

  subroutine commas_to_spaces(line)
    character(len=*),intent(inout)::line
    integer::j
    do j=1,len_trim(line);if(line(j:j)==',')line(j:j)=' ';end do
  end subroutine commas_to_spaces
end program fit_csv
