! SPDX-License-Identifier: GPL-3.0-only
program test_estimation
    use rsdc, only: dp, rsdc_model, rsdc_control, rsdc_estimate, rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc, only: rsdc_simulation_result, rsdc_simulate_fixed
    implicit none
    integer, parameter :: t=180,k=2,n=2
    real(dp) :: p(n,n), mu(n,k), cov(k,k,n), x(t,1)
    type(rsdc_simulation_result) :: sim
    type(rsdc_model) :: fit
    type(rsdc_control) :: ctl
    logical :: ok
    p=reshape([0.95_dp,0.08_dp,0.05_dp,0.92_dp],[2,2])
    mu=0.0_dp
    cov(:,:,1)=reshape([1.0_dp,0.10_dp,0.10_dp,1.0_dp],[2,2])
    cov(:,:,2)=reshape([1.0_dp,0.75_dp,0.75_dp,1.0_dp],[2,2])
    call rsdc_simulate_fixed(t,p,mu,cov,sim,9,ok)
    call check(ok,'simulation')
    ctl%population_size=35; ctl%max_global_iterations=90; ctl%max_local_iterations=60
    ctl%initial_step=0.06_dp; ctl%seed=11
    call rsdc_estimate(rsdc_nox,sim%observations,n,fit,control=ctl,ok=ok)
    call check(ok,'noX estimation')
    call check(fit%log_likelihood>-1.0e8_dp,'finite fit')
    call check(fit%correlations(1,1)<fit%correlations(2,1),'ordered regimes')
    call check(maxval(abs(sum(fit%transition_matrix,dim=2)-1.0_dp))<1.0e-10_dp,'transition rows')
    if (allocated(ctl%start)) deallocate(ctl%start)
    allocate(ctl%start(size(fit%parameters))); ctl%start=fit%parameters
    ctl%max_local_iterations=12
    call rsdc_estimate(rsdc_nox,sim%observations,n,fit,control=ctl,ok=ok)
    call check(ok,'warm-start estimation')
    deallocate(ctl%start)
    x = 1.0_dp
    ctl%max_global_iterations = 55; ctl%max_local_iterations = 35
    call rsdc_estimate(rsdc_tvtp,sim%observations,n,fit,x=x,control=ctl,ok=ok)
    call check(ok .and. allocated(fit%beta),'tvtp estimation')
    call rsdc_estimate(rsdc_const,sim%observations,1,fit,control=ctl,ok=ok)
    call check(ok .and. abs(fit%correlations(1,1))<1.0_dp,'constant estimation')
    print '(a)', 'test_estimation: PASS'
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine check
end program test_estimation
