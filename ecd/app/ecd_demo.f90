! SPDX-License-Identifier: Artistic-2.0
program ecd_demo
  use ecd_api
  implicit none
  type(ecd_model) :: cusp
  type(ecld_model) :: ld
  type(ecd_stats_type) :: stats
  real(dp) :: price,iv,q95
  integer :: status

  cusp=ecd_new(with_stats=.true.,status=status)
  if(status/=ecd_ok) error stop 'failed to initialize standard cusp model'
  stats=ecd_statistics(cusp)
  q95=ecd_quantile(cusp,0.95_dp,status)

  ld=ecld_new(lambda=3.0_dp,sigma=0.4_dp)
  price=bs_call_price(0.25_dp,100.0_dp,100.0_dp,0.5_dp,0.03_dp,0.01_dp)
  iv=bs_implied_volatility(price,100.0_dp,100.0_dp,0.5_dp,0.03_dp,0.01_dp,'c',status)

  write(*,'(a)') 'ecd-fortran demonstration'
  write(*,'(a,f18.10)') 'standard cusp normalizing constant: ',cusp%norm_const
  write(*,'(a,f18.10)') 'standard cusp variance:             ',stats%variance
  write(*,'(a,f18.10)') 'standard cusp 95% quantile:         ',q95
  write(*,'(a,f18.10)') 'ECLD(lambda=3) CDF at 0.2:          ',ecld_cdf(ld,0.2_dp)
  write(*,'(a,f18.10)') 'Black-Scholes call price:           ',price
  write(*,'(a,f18.10)') 'recovered implied volatility:       ',iv
end program ecd_demo
