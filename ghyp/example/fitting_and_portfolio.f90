! SPDX-License-Identifier: GPL-2.0-or-later
program fitting_and_portfolio
   use ghyp
   implicit none
   real(dp) :: data(12), scatter(2,2)
   type(fit_result) :: fit
   type(ghyp_model_type) :: market
   type(portfolio_result) :: portfolio

   data=[-1.1_dp,-0.7_dp,-0.4_dp,-0.2_dp,0.0_dp,0.1_dp, &
          0.2_dp,0.4_dp,0.7_dp,0.9_dp,1.3_dp,1.7_dp]
   fit=fit_gaussian_uv(data)
   if(.not.fit%ok)error stop trim(fit%message)
   print '(a,f10.5,a,f10.5)', 'Gaussian fit mean=',fit%model%mu(1), &
      ' sd=',sqrt(fit%model%scatter(1,1))

   scatter=reshape([1.0_dp,0.25_dp,0.25_dp,0.7_dp],[2,2])
   market=ghyp_mv(0.8_dp,1.5_dp,2.2_dp,[0.05_dp,0.02_dp],scatter,[0.1_dp,-0.05_dp])
   portfolio=portfolio_optimize(market,'sd','minimum.risk')
   if(.not.portfolio%ok)error stop trim(portfolio%message)
   print '(a,2f10.5)', 'minimum-risk weights: ',portfolio%weights
end program fitting_and_portfolio
