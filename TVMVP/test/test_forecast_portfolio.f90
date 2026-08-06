program test_forecast_portfolio
  use tvmvp, only : dp, tvmvp_error, portfolio_result, portfolio_prediction_result, &
                    comp_expected_returns, minimum_variance_portfolio, maximum_sharpe_portfolio, &
                    constrained_minimum_variance, predict_portfolio
  implicit none
  integer, parameter :: n=50,p=4
  real(dp) :: sigma(2,2),mu(2),ret(n,p),grid(6)
  real(dp), allocatable :: expected(:)
  integer, allocatable :: orders(:,:)
  integer :: i,j
  type(tvmvp_error) :: err
  type(portfolio_result) :: port
  type(portfolio_prediction_result) :: pred
  sigma=reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2]); mu=[0.01_dp,0.02_dp]
  call minimum_variance_portfolio(sigma,mu,5,port)
  call check(.not.port%error%failed(),'MVP status')
  call check(abs(sum(port%weights)-1.0_dp)<1.0e-12_dp,'MVP budget')
  call check(maxval(abs(port%weights-[0.7272727272727273_dp,0.2727272727272727_dp]))<1.0e-10_dp,'MVP weights')
  call maximum_sharpe_portfolio(sigma,mu,5,port)
  call check(.not.port%error%failed() .and. abs(sum(port%weights)-1.0_dp)<1.0e-12_dp,'max Sharpe')
  call constrained_minimum_variance(sigma,mu,5,0.08_dp,port)
  call check(.not.port%error%failed(),'constrained status')
  call check(abs(dot_product(port%weights,mu)-0.016_dp)<1.0e-10_dp,'return constraint')
  do j=1,p
    ret(1,j)=0.001_dp*real(j,dp)
    do i=2,n
      ret(i,j)=0.0005_dp*real(j,dp)+0.72_dp*(ret(i-1,j)-0.0005_dp*real(j,dp))+ &
               0.002_dp*sin(real(i*(j+1),dp))
    end do
  end do
  call comp_expected_returns(ret,4,expected,err,orders)
  call check(.not.err%failed() .and. size(expected)==p,'expected-return forecast')
  call check(all(abs(expected)<0.1_dp),'forecast magnitude')
  do i=1,6
    grid(i)=0.01_dp+0.15_dp*real(i-1,dp)
  end do
  call predict_portfolio(ret,3,pred,m=1,compute_max_sharpe=.true.,min_return=0.004_dp,m0=4,rho_grid=grid)
  call check(.not.pred%error%failed(),'end-to-end prediction')
  call check(pred%minimum_variance%available .and. pred%maximum_sharpe%available .and. &
             pred%return_constrained%available,'portfolio availability')
  call check(abs(sum(pred%minimum_variance%weights)-1.0_dp)<1.0e-9_dp,'predicted MVP budget')
  print *, 'test_forecast_portfolio: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_forecast_portfolio
