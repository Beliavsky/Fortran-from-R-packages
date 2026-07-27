! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program analyze_csv
  use kinds_mod, only: dp
  use csv_mod, only: read_dated_csv
  use statistics_mod, only: mean_value, sd_value, skewness_value, kurtosis_value
  use returns_mod, only: annualized_return, cumulative_return
  use drawdown_mod, only: max_drawdown, pain_index, ulcer_index
  use risk_mod, only: value_at_risk, expected_shortfall
  use performance_ratios_mod, only: annualized_sharpe_ratio, sortino_ratio, omega_ratio
  use performance_ratios_mod, only: calmar_ratio, information_ratio, tracking_error
  use capm_mod, only: sfm_result, sfm_fit
  implicit none
  character(len=1024)::filename,arg
  real(dp),allocatable::data(:,:),rf(:)
  real(dp)::scale
  integer::status,return_col,benchmark_col,n
  type(sfm_result)::fit

  if(command_argument_count()<1)then
    write(*,'(a)')'usage: analyze_csv FILE [return_column] [benchmark_column] [scale]'
    stop 2
  end if
  call get_command_argument(1,filename)
  return_col=1;benchmark_col=0;scale=252.0_dp
  if(command_argument_count()>=2)then;call get_command_argument(2,arg);read(arg,*)return_col;end if
  if(command_argument_count()>=3)then;call get_command_argument(3,arg);read(arg,*)benchmark_col;end if
  if(command_argument_count()>=4)then;call get_command_argument(4,arg);read(arg,*)scale;end if
  call read_dated_csv(trim(filename),data,status)
  if(status/=0)then;write(*,'(a,i0)')'CSV read error: ',status;stop 1;end if
  if(return_col<1 .or. return_col>size(data,2))then;write(*,'(a)')'invalid return column';stop 1;end if
  n=size(data,1);allocate(rf(n));rf=0.0_dp
  write(*,'(a,i0)')'observations: ',n
  write(*,'(a,es16.8)')'mean: ',mean_value(data(:,return_col))
  write(*,'(a,es16.8)')'standard deviation: ',sd_value(data(:,return_col))
  write(*,'(a,es16.8)')'skewness: ',skewness_value(data(:,return_col),1)
  write(*,'(a,es16.8)')'excess kurtosis: ',kurtosis_value(data(:,return_col),1,.true.)
  write(*,'(a,es16.8)')'cumulative return: ',cumulative_return(data(:,return_col))
  write(*,'(a,es16.8)')'annualized return: ',annualized_return(data(:,return_col),scale)
  write(*,'(a,es16.8)')'annualized Sharpe: ',annualized_sharpe_ratio(data(:,return_col),rf,scale,.false.)
  write(*,'(a,es16.8)')'Sortino: ',sortino_ratio(data(:,return_col),0.0_dp)
  write(*,'(a,es16.8)')'Omega: ',omega_ratio(data(:,return_col),0.0_dp)
  write(*,'(a,es16.8)')'95% modified VaR: ',value_at_risk(data(:,return_col),0.95_dp,'modified')
  write(*,'(a,es16.8)')'95% modified ES: ',expected_shortfall(data(:,return_col),0.95_dp,'modified')
  write(*,'(a,es16.8)')'maximum drawdown: ',max_drawdown(data(:,return_col))
  write(*,'(a,es16.8)')'pain index: ',pain_index(data(:,return_col))
  write(*,'(a,es16.8)')'ulcer index: ',ulcer_index(data(:,return_col))
  write(*,'(a,es16.8)')'Calmar ratio: ',calmar_ratio(data(:,return_col),scale)
  if(benchmark_col>=1 .and. benchmark_col<=size(data,2))then
    call sfm_fit(data(:,return_col),data(:,benchmark_col),rf,fit)
    write(*,'(a,es16.8)')'CAPM alpha: ',fit%alpha
    write(*,'(a,es16.8)')'CAPM beta: ',fit%beta
    write(*,'(a,es16.8)')'CAPM R-squared: ',fit%r_squared
    write(*,'(a,es16.8)')'tracking error: ',tracking_error(data(:,return_col),data(:,benchmark_col),scale)
    write(*,'(a,es16.8)')'information ratio: ',information_ratio(data(:,return_col),data(:,benchmark_col),scale)
  end if
end program analyze_csv
