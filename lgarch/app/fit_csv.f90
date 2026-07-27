! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program fit_csv
  use lgarch_kinds, only : dp
  use lgarch_univariate, only : fit_lgarch,lgarch_fit_result,LGARCH_LS,LGARCH_ML,LGARCH_CEX2
  use lgarch_multivariate, only : fit_mlgarch,mlgarch_fit_result
  implicit none
  character(len=512)::filename,mode,arg,line,date_token
  integer::argc,arch,garch,method,n,m,unit,ios,i
  real(dp),allocatable::data(:,:)
  type(lgarch_fit_result)::ufit
  type(mlgarch_fit_result)::mfit
  argc=command_argument_count()
  if(argc<2) then
    print '(a)','Usage: fit_csv FILE univariate [ls|ml|cex2] [arch] [garch]'
    print '(a)','   or: fit_csv FILE multivariate [arch] [garch]'
    stop 2
  end if
  call get_command_argument(1,filename); call get_command_argument(2,mode)
  arch=1; garch=1; method=LGARCH_LS
  if(trim(mode)=='univariate') then
    if(argc>=3) then
      call get_command_argument(3,arg)
      select case(trim(arg)); case('ls'); method=LGARCH_LS; case('ml'); method=LGARCH_ML; case('cex2'); method=LGARCH_CEX2; case default; error stop 'Unknown method'; end select
    end if
    if(argc>=4) then; call get_command_argument(4,arg); read(arg,*)arch; end if
    if(argc>=5) then; call get_command_argument(5,arg); read(arg,*)garch; end if
  else if(trim(mode)=='multivariate') then
    if(argc>=3) then; call get_command_argument(3,arg); read(arg,*)arch; end if
    if(argc>=4) then; call get_command_argument(4,arg); read(arg,*)garch; end if
  else
    error stop 'Mode must be univariate or multivariate'
  end if
  open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios); if(ios/=0) error stop 'Cannot open CSV file'
  read(unit,'(a)',iostat=ios)line; if(ios/=0) error stop 'CSV is empty'
  m=count_commas(trim(line))
  n=0
  do
    read(unit,'(a)',iostat=ios)line; if(ios/=0) exit
    if(len_trim(line)>0)n=n+1
  end do
  rewind(unit); read(unit,'(a)')line; allocate(data(n,m)); i=0
  do
    read(unit,'(a)',iostat=ios)line; if(ios/=0) exit
    if(len_trim(line)==0)cycle
    i=i+1; read(line,*,iostat=ios)date_token,data(i,:); if(ios/=0) error stop 'Invalid CSV row'
  end do
  close(unit)
  if(trim(mode)=='univariate') then
    if(m/=1) error stop 'Univariate mode requires exactly one numeric series after Date'
    call fit_lgarch(data(:,1),arch,garch,method,ufit,compute_vcov=.true.)
    print '(a,l1)','Converged: ',ufit%converged
    print '(a,*(1x,es14.6))','Log-GARCH parameters:',ufit%lgarch_par
    print '(a,es16.8)','Model log likelihood: ',ufit%loglik_model
    print '(a,es16.8)','ARMA residual sum of squares: ',ufit%rss
  else
    call fit_mlgarch(data,arch,garch,mfit,compute_vcov=.true.)
    print '(a,l1)','Converged: ',mfit%converged
    print '(a,*(1x,es14.6))','Multivariate log-GARCH parameters:',mfit%mlgarch_par
    print '(a,es16.8)','Model log likelihood: ',mfit%loglik_model
    print '(a,es16.8)','VARMA log likelihood: ',mfit%objective_varma
  end if
contains
  integer function count_commas(s) result(nc)
    character(len=*),intent(in)::s
    integer::j
    nc=0; do j=1,len_trim(s); if(s(j:j)==',')nc=nc+1; end do
  end function count_commas
end program fit_csv
