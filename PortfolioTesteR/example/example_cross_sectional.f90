! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program example_cross_sectional
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),mom(:,:),ranked(:,:),selection(:,:),breadth(:)
  call generate_sample_prices(100,8,prices,619_i8)
  call calc_momentum(prices,12,mom)
  call calc_relative_strength_rank(mom,ranked,'percentile')
  call filter_top_n(ranked,3,selection)
  call calc_market_breadth(selection,breadth)
  print '(a,f8.4)','last-period selected fraction: ',breadth(size(breadth))
  print '(a,8f8.3)','last relative-strength percentiles: ',ranked(size(ranked,1),:)
end program example_cross_sectional
