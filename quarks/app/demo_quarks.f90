program demo_quarks
   use iso_fortran_env, only : int64
   use quarks
   implicit none
   real(dp) :: returns(260), losses(40)
   type(rng_state) :: rng
   type(risk_result) :: current
   type(rollcast_result) :: rolling
   type(coverage_result) :: coverage
   type(loss_result) :: scores
   integer :: i

   do i = 1, size(returns)
      returns(i) = 0.009_dp * sin(0.23_dp * real(i, dp)) + &
         0.004_dp * cos(0.071_dp * real(i, dp))
   end do
   returns(232) = -0.045_dp
   returns(247) = -0.035_dp
   call seed_rng(rng, 20260802_int64)
   current = fhs(returns(1:220), p=0.95_dp, model=volatility_ewma, &
      nboot=10000, rng=rng)
   rolling = rollcast(returns, p=0.95_dp, method=method_age, lambda=0.98_dp, &
      nout=40, nwin=120)
   losses = -rolling%xout
   coverage = cvgtest(losses, rolling%var, rolling%p)
   scores = lossfun(losses, rolling%es)

   print '(a)', 'quarks-fortran demonstration'
   print '(a,f12.6)', 'current filtered VaR: ', current%var
   print '(a,f12.6)', 'current filtered ES : ', current%es
   print '(a,i0)', 'rolling violations   : ', coverage%violations
   print '(a,f10.6)', 'coverage p-value     : ', coverage%p_uc
   print '(a,f12.4)', 'ES loss function 1  : ', scores%lossfun1
end program demo_quarks
