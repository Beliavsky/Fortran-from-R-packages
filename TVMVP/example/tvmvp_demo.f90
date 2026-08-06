program tvmvp_demo
  use tvmvp, only : dp, portfolio_prediction_result, expanding_window_result, predict_portfolio, expanding_tvmvp
  implicit none
  integer, parameter :: n=70,p=6
  real(dp) :: returns(n,p),grid(8)
  integer :: i,j
  type(portfolio_prediction_result) :: prediction
  type(expanding_window_result) :: rolling
  do i=1,n
    do j=1,p
      returns(i,j)=0.0004_dp*real(j,dp)+0.015_dp*sin(0.10_dp*real(i,dp))*(0.4_dp+0.15_dp*j)+ &
                   0.005_dp*cos(0.19_dp*real(i*j,dp))
    end do
  end do
  do i=1,size(grid)
    grid(i)=0.01_dp+0.12_dp*real(i-1,dp)
  end do
  call predict_portfolio(returns,21,prediction,m=1,compute_max_sharpe=.true., &
                         min_return=0.03_dp,m0=5,rho_grid=grid)
  if (prediction%error%failed()) then
    print *, trim(prediction%error%message)
    error stop 1
  end if
  print '(a,*(f10.5,1x))','MVP weights: ',prediction%minimum_variance%weights
  print '(a,f10.6)','MVP expected return: ',prediction%minimum_variance%expected_return
  print '(a,f10.6)','MVP risk: ',prediction%minimum_variance%risk
  call expanding_tvmvp(returns,40,10,1,rolling,m0=4,rho_grid=grid,source_compatible_expanding=.false.)
  if (rolling%error%failed()) then
    print *, trim(rolling%error%message)
    error stop 1
  end if
  print '(a,f10.6)','Rolling TV-MVP Sharpe: ',rolling%tvmvp_metrics%sharpe
  print '(a,f10.6)','Rolling equal-weight Sharpe: ',rolling%equal_metrics%sharpe
end program tvmvp_demo
