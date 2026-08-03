! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program example_ml_scores
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),mom(:,:),vol(:,:),labels(:,:),scores(:,:),ic(:)
  real(dp),allocatable::features(:,:,:),lag_mom(:,:),lag_vol(:,:)
  call generate_sample_prices(180,10,prices,8127_i8)
  call calc_momentum(prices,12,mom)
  call calc_rolling_volatility(prices,12,vol)
  call panel_lag(mom,1,lag_mom);call panel_lag(vol,1,lag_vol)
  call make_labels(prices,4,1,labels)
  allocate(features(180,10,2));features(:,:,1)=lag_mom;features(:,:,2)=lag_vol
  call rolling_fit_predict(features,labels,78,4,4,scores,lambda=1.0_dp)
  call ic_series(scores,labels,ic,.true.)
  print '(a,f10.4)','mean out-of-sample rank IC: ',finite_mean(ic)
  print '(a,i0)','finite forecasts: ',count(is_finite(scores))
end program example_ml_scores
