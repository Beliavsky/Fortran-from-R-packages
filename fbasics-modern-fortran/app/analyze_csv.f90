! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program analyze_csv
  use fbasics
  implicit none
  character(len=512) :: filename, family, line
  real(dp), allocatable :: x(:), temp(:)
  real(dp) :: value
  integer :: unit, ios, n, capacity
  logical :: ok
  type(basic_stats_result) :: stats
  type(distribution_fit) :: fit
  type(stable_fit_result) :: stable_fit
  type(spline_density_fit) :: density_fit
  type(test_result) :: jb
  integer :: output_kind
  if(command_argument_count()<1)then
    write(*,'(a)')'usage: analyze_csv FILE [normal|student|nig|gld|fmkl|fm5|gh|hyp|sgh|snig|stable|stable-mle|ssd]'
    error stop 2
  end if
  call get_command_argument(1,filename)
  family='student'
  if(command_argument_count()>=2)call get_command_argument(2,family)
  capacity=1024;n=0;allocate(x(capacity))
  open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
  if(ios/=0)then;write(*,'(a)')'cannot open '//trim(filename);error stop 2;end if
  do
    read(unit,'(a)',iostat=ios)line
    if(ios/=0)exit
    call parse_last_numeric(line,value,ok)
    if(.not.ok)cycle
    n=n+1
    if(n>capacity)then
      capacity=2*capacity;allocate(temp(capacity));temp(1:n-1)=x(1:n-1);call move_alloc(temp,x)
    end if
    x(n)=value
  end do
  close(unit)
  if(n<5)then;write(*,'(a)')'fewer than five numeric observations';error stop 2;end if
  x=x(1:n)
  stats=basic_stats(x);jb=jarque_bera_test(x);output_kind=0
  select case(trim(adjustl(family)))
  case('normal');call fit_normal(x,fit)
  case('student','t');call fit_student(x,fit)
  case('nig');call fit_nig(x,fit)
  case('gld');call fit_gld_quantiles(x,fit)
  case('fmkl');call fit_gld_extended(x,'fmkl','rob',fit,max_iter=350)
  case('fm5');call fit_gld_extended(x,'fm5','rob',fit,max_iter=350)
  case('gh');call fit_gh(x,fit,max_iter=180)
  case('hyp');call fit_hyp(x,fit,max_iter=180)
  case('sgh');call fit_sgh(x,fit,max_iter=180)
  case('snig');call fit_snig(x,fit,max_iter=180)
  case('stable');call fit_stable(x,'ecf',stable_fit,max_iter=450);output_kind=1
  case('stable-mle');call fit_stable(x,'mle',stable_fit,max_iter=120);output_kind=1
  case('ssd');call fit_spline_density(x,density_fit,n_basis=12,grid_size=401,max_iter=100);output_kind=2
  case default;write(*,'(a)')'unknown family: '//trim(family);error stop 2
  end select
  write(*,'(a,i0)')'nobs: ',n
  write(*,'(a,es16.8)')'mean: ',stats%mean
  write(*,'(a,es16.8)')'stdev: ',stats%stdev
  write(*,'(a,es16.8)')'skewness: ',stats%skewness
  write(*,'(a,es16.8)')'excess kurtosis: ',stats%kurtosis
  write(*,'(a,es16.8)')'Jarque-Bera p-value: ',jb%p_value
  select case(output_kind)
  case(1)
    write(*,'(a)')'family: stable-S1'
    write(*,'(a,a)')'method: ',trim(stable_fit%method)
    write(*,'(a,4(1x,es16.8))')'parameters (alpha beta gamma delta):',stable_fit%alpha,stable_fit%beta,stable_fit%gamma,stable_fit%delta
    write(*,'(a,es16.8)')'objective: ',stable_fit%objective
    if(stable_fit%method=='mle')write(*,'(a,es16.8)')'log likelihood: ',stable_fit%loglik
    write(*,'(a,l1)')'converged: ',stable_fit%converged
  case(2)
    write(*,'(a)')'family: spline-density'
    write(*,'(a,es16.8)')'lambda: ',density_fit%lambda
    write(*,'(a,es16.8)')'log likelihood: ',density_fit%loglik
    write(*,'(a,i0)')'grid points: ',size(density_fit%grid)
    write(*,'(a,l1)')'converged: ',density_fit%converged
  case default
    write(*,'(a,a)')'family: ',trim(fit%family)
    write(*,'(a,*(1x,es16.8))')'parameters:',fit%parameters
    write(*,'(a,es16.8)')'log likelihood/objective score: ',fit%loglik
    write(*,'(a,l1)')'converged: ',fit%converged
  end select
contains
  subroutine parse_last_numeric(text,v,success)
    character(len=*),intent(in)::text
    real(dp),intent(out)::v
    logical,intent(out)::success
    integer::pos,stat
    character(len=:),allocatable::field
    pos=scan(trim(text),',',back=.true.)
    if(pos>0)then;field=adjustl(text(pos+1:));else;field=adjustl(text);end if
    read(field,*,iostat=stat)v
    success=stat==0
  end subroutine parse_last_numeric
end program analyze_csv
