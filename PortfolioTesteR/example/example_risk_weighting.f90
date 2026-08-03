! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program example_risk_weighting
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),returns(:,:),cov(:,:),hrp(:),erc(:),maxdiv(:)
  integer::status
  call generate_sample_prices(260,6,prices,9182_i8)
  call panel_returns_log(prices,returns)
  call covariance_matrix(returns(2:,:),cov)
  call calculate_hrp_weights(returns(2:,:),hrp,status)
  call calculate_erc_weights(cov,erc,status)
  call calculate_max_div_weights(cov,maxdiv,status)
  print '(a,6f9.4)','HRP:     ',hrp
  print '(a,6f9.4)','ERC:     ',erc
  print '(a,6f9.4)','Max div: ',maxdiv
end program example_risk_weighting
