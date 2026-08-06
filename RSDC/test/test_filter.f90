! SPDX-License-Identifier: GPL-3.0-only
program test_filter
    use rsdc, only: dp, rsdc_filter_result, rsdc_hamilton
    use rsdc, only: rsdc_negative_log_likelihood, rsdc_nox, rsdc_tvtp
    implicit none
    integer, parameter :: t=80, k=2, n=2
    real(dp) :: y(t,k), rho(n,1), p(n,n), params(4), x(t,1), beta(n,1)
    type(rsdc_filter_result) :: out
    integer :: i
    do i=1,t
        y(i,1)=0.6_dp*sin(0.17_dp*i)
        y(i,2)=0.5_dp*cos(0.13_dp*i)+0.2_dp*y(i,1)
        x(i,1)=1.0_dp
    end do
    rho(:,1)=[0.1_dp,0.75_dp]
    p=reshape([0.94_dp,0.08_dp,0.06_dp,0.92_dp],[2,2])
    call rsdc_hamilton(y,rho,out,pmat=p)
    call check(out%ok,'fixed filter')
    call check(maxval(abs(sum(out%filtered,dim=1)-1.0_dp))<1.0e-12_dp,'filtered sums')
    call check(maxval(abs(sum(out%smoothed,dim=1)-1.0_dp))<1.0e-12_dp,'smoothed sums')
    params=[0.94_dp,0.92_dp,0.1_dp,0.75_dp]
    call check(abs(rsdc_negative_log_likelihood(params,y,rsdc_nox,n)+out%log_likelihood)<1.0e-10_dp, &
        'likelihood parity')
    beta(:,1)=[2.2_dp,1.8_dp]
    call rsdc_hamilton(y,rho,out,x=x,beta=beta)
    call check(out%ok,'tvtp filter')
    params=[2.2_dp,1.8_dp,0.1_dp,0.75_dp]
    call check(rsdc_negative_log_likelihood(params,y,rsdc_tvtp,n,x)<1.0e9_dp,'tvtp likelihood')
    print '(a)', 'test_filter: PASS'
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine check
end program test_filter
