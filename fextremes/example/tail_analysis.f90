! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program tail_analysis
  use fextremes_kinds, only: dp
  use fextremes_csv, only: read_numeric_series
  use fextremes_stats, only: quantile_type1
  use fextremes_fit, only: gpd_fit_result,fit_gpd
  use fextremes_diagnostics, only: curve_result,tail_index_result,mean_excess_curve,hill_estimator
  implicit none
  real(dp),allocatable::x(:)
  logical::ok
  real(dp)::u
  type(gpd_fit_result)::fit
  type(curve_result)::me
  type(tail_index_result)::hill
  call read_numeric_series('data/danishClaims.csv',x,ok)
  if(.not.ok) error stop 'Could not read data/danishClaims.csv'
  u=quantile_type1(x,0.90_dp); call fit_gpd(x,u,fit,'mle','observed')
  call mean_excess_curve(x,me); call hill_estimator(x,0.10_dp,hill)
  print '(a,f10.4)','threshold: ',u
  print '(a,2f12.6)','GPD xi beta: ',fit%xi,fit%beta
  print '(a,f12.6)','Hill mean: ',hill%mean
  print '(a,i0)','mean-excess curve points: ',size(me%x)
end program tail_analysis
