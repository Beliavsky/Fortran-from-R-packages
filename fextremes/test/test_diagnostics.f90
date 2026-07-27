! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program test_diagnostics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state,seed_rng,uniform_rng,exponential_rng
  use fextremes_diagnostics
  use fextremes_metrics, only: ema_filter,riskmetrics_volatility
  implicit none
  real(dp),allocatable::x(:),ratios(:,:),path(:),ema(:),vol(:)
  type(rng_state)::rng
  type(curve_result)::curve
  type(record_result)::rec
  type(exceedance_acf_result)::ac
  type(tail_index_result)::hill,pick,deh
  integer::i
  integer,allocatable::record_counts(:)
  call seed_rng(rng,12345); allocate(x(5000))
  do i=1,size(x); x(i)=exponential_rng(rng); end do
  call empirical_survival(x,curve)
  call pareto_qq(x,0.0_dp,curve)
  call assert_true(all(curve%y(2:)>=curve%y(:size(curve%y)-1)),'Pareto QQ')
  call assert_true(all(curve%x(2:)>=curve%x(:size(x)-1)),'empirical survival sorting')
  call mean_excess_curve(x,curve)
  call assert_true(abs(sum(curve%y(1:min(100,size(curve%y)))) / &
    real(min(100,size(curve%y)),dp)-1.0_dp) < 0.25_dp, &
    'exponential mean excess')
  call mean_residual_life(x,0.2_dp,2.0_dp,20,0.95_dp,curve)
  call assert_true(count(curve%y>0.0_dp)>15,'mean residual life')
  call records_development(x,rec)
  call subsample_record_counts(x,10,record_counts)
  call assert_true(all(record_counts>=1),'subsample records')
  call assert_true(all(rec%record(2:)>rec%record(:size(rec%record)-1)),'records increasing')
  call max_sum_ratios(x,[1.0_dp,2.0_dp],ratios)
  call assert_true(all(ratios>=0.0_dp .and. ratios<=1.0_dp),'max sum ratios')
  call slln_path(x,path); call assert_true(abs(path(size(path))-1.0_dp)<0.06_dp,'SLLN path')
  call lil_path(x,path)
  call assert_true(size(path)==size(x)-2,'LIL path')
  call exceedance_acf(x,3.0_dp,10,ac); call assert_true(abs(ac%height_acf(0)-1.0_dp)<1.0e-12_dp,'exceedance ACF')
  do i=1,size(x); x(i)=uniform_rng(rng)**(-0.25_dp); end do
  call hill_estimator(x,0.1_dp,hill); call pickands_estimator(x,0.2_dp,pick); call dehaan_estimator(x,0.1_dp,deh)
  call assert_true(abs(hill%mean-0.25_dp)<0.12_dp,'Hill estimate')
  call assert_true(ieee_is_finite(pick%mean) .and. ieee_is_finite(deh%mean), &
    'tail estimators finite')
  allocate(ema(size(x)),vol(size(x)))
  call ema_filter(x,0.1_dp,ema)
  call riskmetrics_volatility(x-1.0_dp,0.94_dp,vol)
  call assert_true(all(vol>=0.0_dp),'RiskMetrics volatility')
  call assert_true(abs(normal_mean_excess(0.0_dp,0.0_dp,1.0_dp)-0.7978845608028654_dp)<1.0e-12_dp,'normal mean excess')
  print '(a)','Exploratory and tail-index tests passed.'
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond; character(len=*),intent(in)::msg
    if(.not.cond) then; print *,trim(msg); error stop 1; end if
  end subroutine
end program test_diagnostics
