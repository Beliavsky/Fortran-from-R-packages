! SPDX-License-Identifier: GPL-3.0-or-later
program distributions_and_risk
  use nvmix
  implicit none
  type(nvmix_model) :: model
  real(dp) :: loc(1),scale(1,1)
  loc=0.0_dp; scale(1,1)=1.0_dp
  model=make_nvmix_model(loc,scale,mix_inverse_gamma,8.0_dp)
  print '(a,f12.6)','99% VaR: ',VaR_nvmix(0.99_dp,model)
  print '(a,f12.6)','99% ES:  ',ES_nvmix(0.99_dp,model)
  print '(a,f12.6)','P(X <= 1): ',nvmix_cdf_1d(1.0_dp,model)
end program
