! SPDX-License-Identifier: GPL-3.0-only
program filter_and_portfolio
    use rsdc, only: dp, rsdc_filter_result, rsdc_portfolio_result
    use rsdc, only: rsdc_hamilton, rsdc_minvar
    implicit none
    integer, parameter :: nobs=100, k=2, n=2
    real(dp) :: y(nobs,k), rho(n,1), p(n,n), sigma(nobs,k), path(nobs,1)
    type(rsdc_filter_result) :: filter
    type(rsdc_portfolio_result) :: portfolio
    integer :: t, s
    logical :: ok

    do t=1,nobs
        y(t,1)=0.01_dp*sin(0.15_dp*t)
        y(t,2)=0.01_dp*(0.4_dp*sin(0.15_dp*t)+cos(0.09_dp*t))
    end do
    rho(:,1)=[0.10_dp,0.70_dp]
    p=reshape([0.94_dp,0.10_dp,0.06_dp,0.90_dp],[2,2])
    call rsdc_hamilton(y,rho,filter,pmat=p)
    if(.not.filter%ok) error stop 'filter failed'
    path=0.0_dp
    do t=1,nobs
        do s=1,n
            path(t,1)=path(t,1)+filter%smoothed(s,t)*rho(s,1)
        end do
    end do
    sigma=0.20_dp
    call rsdc_minvar(sigma,path,y,portfolio,long_only=.true.,lagged=.true.,ok=ok)
    if(.not.ok) error stop 'portfolio failed'
    print '(a,f10.6)', 'realized volatility: ',portfolio%realized_volatility
end program filter_and_portfolio
