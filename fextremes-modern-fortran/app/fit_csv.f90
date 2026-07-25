! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program fit_csv
  use fextremes_kinds, only: dp
  use fextremes_csv, only: read_numeric_series
  use fextremes_stats, only: quantile_type1
  use fextremes_fit, only: gev_fit_result,gpd_fit_result,fit_gev,fit_gumbel,fit_gpd
  use fextremes_risk, only: risk_result,gpd_risk_measures
  implicit none
  character(len=512)::filename,mode,arg
  real(dp),allocatable::x(:)
  real(dp)::threshold
  logical::ok
  type(gev_fit_result)::gf
  type(gpd_fit_result)::pf
  type(risk_result)::risk
  integer::i
  if(command_argument_count()<2) then
    print '(a)','Usage: fit_csv FILE {gev|gumbel|gpd} [threshold]'
    stop 2
  end if
  call get_command_argument(1,filename); call get_command_argument(2,mode)
  call read_numeric_series(trim(filename),x,ok)
  if(.not.ok) error stop 'Could not read numeric data.'
  select case(trim(mode))
  case('gev')
    call fit_gev(x,gf,'mle')
    print '(a,l1)','converged: ',gf%converged
    print '(a,3(1x,es16.8))','xi mu beta:',gf%xi,gf%mu,gf%beta
    print '(a,1x,es16.8)','negative log likelihood:',gf%nll
  case('gumbel')
    call fit_gumbel(x,gf,'mle')
    print '(a,l1)','converged: ',gf%converged
    print '(a,2(1x,es16.8))','mu beta:',gf%mu,gf%beta
    print '(a,1x,es16.8)','negative log likelihood:',gf%nll
  case('gpd')
    threshold=quantile_type1(x,0.95_dp)
    if(command_argument_count()>=3) then
      call get_command_argument(3,arg); read(arg,*) threshold
    end if
    call fit_gpd(x,threshold,pf,'mle','observed')
    print '(a,l1)','converged: ',pf%converged
    print '(a,3(1x,es16.8))','threshold xi beta:',threshold,pf%xi,pf%beta
    call gpd_risk_measures(pf,[0.99_dp,0.995_dp,0.999_dp],risk)
    do i=1,size(risk%probability)
      print '(f8.4,2(1x,es16.8))',risk%probability(i),risk%value_at_risk(i),risk%expected_shortfall(i)
    end do
  case default
    error stop 'Unknown mode.'
  end select
end program fit_csv
