! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program fit_csv
  use tsdyn
  implicit none
  character(len=1024)::filename,model_name,arg
  real(dp),allocatable::data(:,:),fc1(:),fc2(:,:)
  integer::p,info,ncol
  type(ar_model)::ar
  type(setar_model)::tar
  type(lstar_model)::star
  type(var_model)::vm
  type(vecm_model)::cm
  real(dp),allocatable::beta(:,:)

  if(command_argument_count()<2)then
    write(*,'(a)') 'usage: fit_csv FILE MODEL [ORDER]'
    write(*,'(a)') 'MODEL: ar, setar, lstar, var, vecm'
    error stop 2
  end if
  call get_command_argument(1,filename);call get_command_argument(2,model_name)
  p=1
  if(command_argument_count()>=3)then;call get_command_argument(3,arg);read(arg,*,iostat=info)p;if(info/=0)p=1;end if
  call read_dated_csv(trim(filename),data,ncol,info);if(info/=0)error stop 'could not read CSV'
  select case(trim(adjustl(model_name)))
  case('ar')
    call fit_ar(data(:,1),p,include_const,'level',ar,info);if(info/=0)error stop 'AR fit failed'
    call forecast_ar(ar,data(:,1),5,fc1,info);if(info/=0)error stop 'AR forecast failed'
    write(*,'(a,*(f12.6,1x))') 'coefficients ',ar%coefficients
    write(*,'(a,*(f12.6,1x))') 'forecast ',fc1
  case('setar')
    call fit_setar(data(:,1),[p,p],include_const,1,tar,info,ngrid=20,trim=0.10_dp);if(info/=0)error stop 'SETAR fit failed'
    call forecast_setar(tar,data(:,1),5,fc1,info=info);if(info/=0)error stop 'SETAR forecast failed'
    write(*,'(a,f12.6)') 'threshold ',tar%thresholds(1)
    write(*,'(a,*(f12.6,1x))') 'forecast ',fc1
  case('lstar')
    call fit_lstar(data(:,1),p,p,include_const,star,info,ngamma=8,nthreshold=20);if(info/=0)error stop 'LSTAR fit failed'
    call forecast_lstar(star,data(:,1),5,fc1,info);if(info/=0)error stop 'LSTAR forecast failed'
    write(*,'(a,2f12.6)') 'gamma threshold ',star%gamma,star%threshold
    write(*,'(a,*(f12.6,1x))') 'forecast ',fc1
  case('var')
    if(ncol<2)error stop 'VAR requires at least two value columns'
    call fit_var(data,p,include_const,vm,info);if(info/=0)error stop 'VAR fit failed'
    call forecast_var(vm,data,5,fc2,info);if(info/=0)error stop 'VAR forecast failed'
    write(*,'(a,f12.6)') 'loglik ',vm%loglik
    write(*,'(a,*(f12.6,1x))') 'first forecast ',fc2(1,:)
  case('vecm')
    if(ncol<2)error stop 'VECM requires at least two value columns'
    allocate(beta(ncol,1));beta=0.0_dp;beta(1,1)=1.0_dp;beta(2,1)=-1.0_dp
    call fit_vecm(data,p,1,include_const,'fixed',cm,info,beta_fixed=beta);if(info/=0)error stop 'VECM fit failed'
    call forecast_vecm(cm,data,5,fc2,info);if(info/=0)error stop 'VECM forecast failed'
    write(*,'(a,*(f12.6,1x))') 'alpha ',cm%alpha(:,1)
    write(*,'(a,*(f12.6,1x))') 'first forecast ',fc2(1,:)
  case default
    error stop 'unknown model'
  end select
contains
  subroutine read_dated_csv(path,x,k,istat)
    character(len=*),intent(in)::path
    real(dp),allocatable,intent(out)::x(:,:)
    integer,intent(out)::k,istat
    character(len=8192)::line,date_token
    integer::u,ios,nrow,i,j,ncomma
    real(dp),allocatable::row(:)
    open(newunit=u,file=path,status='old',action='read',iostat=ios);if(ios/=0)then;istat=ios;allocate(x(0,0));k=0;return;end if
    read(u,'(a)',iostat=ios)line;if(ios/=0)then;istat=ios;close(u);allocate(x(0,0));k=0;return;end if
    ncomma=count([(line(i:i)==',',i=1,len_trim(line))]);k=ncomma
    if(k<1)then;istat=-1;close(u);allocate(x(0,0));return;end if
    nrow=0
    do;read(u,'(a)',iostat=ios)line;if(ios/=0)exit;if(len_trim(line)>0)nrow=nrow+1;end do
    rewind(u);read(u,'(a)')line;allocate(x(nrow,k),row(k));i=0
    do
      read(u,'(a)',iostat=ios)line;if(ios/=0)exit;if(len_trim(line)==0)cycle
      do j=1,len_trim(line);if(line(j:j)==',')line(j:j)=' ';end do
      read(line,*,iostat=ios)date_token,row;if(ios/=0)then;istat=ios;close(u);return;end if
      i=i+1;x(i,:)=row
    end do
    close(u);istat=0
  end subroutine read_dated_csv
end program fit_csv
