! SPDX-License-Identifier: GPL-3.0-only
program test_portfolio_inference
    use rsdc, only: dp, rsdc_portfolio_result, rsdc_starts_result
    use rsdc, only: rsdc_model, rsdc_diagnostics_result
    use rsdc, only: rsdc_minvar, rsdc_maxdiv, rsdc_make_starts, rsdc_diagnostics
    use rsdc, only: rsdc_nox
    implicit none
    integer, parameter :: t=60,k=3,n=2
    real(dp) :: sigma(t,k), rho(t,3), y(t,k)
    type(rsdc_portfolio_result) :: mv,md
    type(rsdc_starts_result) :: st
    type(rsdc_model) :: model
    type(rsdc_diagnostics_result) :: dg
    logical :: ok
    integer :: i
    do i=1,t
        sigma(i,:)=[0.15_dp,0.20_dp,0.25_dp]
        rho(i,:)=[0.2_dp,0.1_dp,0.3_dp]
        y(i,:)=[sin(0.1_dp*i),cos(0.13_dp*i),sin(0.07_dp*i)]*0.01_dp
    end do
    call rsdc_minvar(sigma,rho,y,mv,long_only=.true.,lagged=.true.,ok=ok)
    call check(ok .and. maxval(abs(sum(mv%weights,dim=2)-1.0_dp))<1.0e-12_dp,'minvar')
    call rsdc_maxdiv(sigma,rho,y,md,long_only=.true.,ok=ok)
    call check(ok .and. all(md%weights>=0.0_dp) .and. md%mean_diversification>0.0_dp,'maxdiv')
    call rsdc_maxdiv(sigma,rho,y,md,long_only=.false.,ok=ok)
    call check(ok .and. maxval(abs(sum(md%weights,dim=2)-1.0_dp))<1.0e-8_dp,'bounded maxdiv')
    call rsdc_make_starts(y,n,rsdc_nox,st,window=15,n_starts=3,ok=ok)
    call check(ok .and. all(shape(st%starts)==[3,8]),'starts')
    model%n_regimes=2; allocate(model%transition_matrix(2,2))
    model%transition_matrix=reshape([0.9_dp,0.2_dp,0.1_dp,0.8_dp],[2,2])
    call rsdc_diagnostics(model,dg)
    call check(abs(sum(dg%ergodic_probability)-1.0_dp)<1.0e-10_dp,'ergodic')
    call check(abs(dg%expected_duration(1)-10.0_dp)<1.0e-10_dp,'duration')
    print '(a)', 'test_portfolio_inference: PASS'
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine check
end program test_portfolio_inference
