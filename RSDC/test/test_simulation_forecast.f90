! SPDX-License-Identifier: GPL-3.0-only
program test_simulation_forecast
    use rsdc, only: dp, rsdc_model, rsdc_simulation_result, rsdc_forecast_result
    use rsdc, only: rsdc_simulate, rsdc_forecast_path, rsdc_forecast_ahead, rsdc_viterbi_path
    use rsdc, only: rsdc_tvtp
    implicit none
    integer, parameter :: t=160,k=2,n=2
    real(dp) :: x(t,1), beta(n,1), mu(n,k), cov(k,k,n), sigma(t,k)
    real(dp), allocatable :: reg(:,:), pc(:,:)
    integer, allocatable :: path(:)
    type(rsdc_simulation_result) :: sim
    type(rsdc_forecast_result) :: fc
    type(rsdc_model) :: model
    logical :: ok
    integer :: correct
    x=1.0_dp; beta(:,1)=[3.0_dp,2.6_dp]; mu=0.0_dp
    cov(:,:,1)=reshape([1.0_dp,0.05_dp,0.05_dp,1.0_dp],[2,2])
    cov(:,:,2)=reshape([1.0_dp,0.80_dp,0.80_dp,1.0_dp],[2,2])
    call rsdc_simulate(t,x,beta,mu,cov,sim,321,ok)
    call check(ok,'simulation')
    model%method=rsdc_tvtp; model%n_regimes=n; model%n_series=k; model%n_covariates=1
    allocate(model%beta(n,1),model%correlations(n,1),model%transition_matrix(n,n),model%covariance(k,k,n))
    model%beta=beta; model%correlations(:,1)=[0.05_dp,0.80_dp]; model%covariance=cov
    model%transition_matrix=reshape([0.952574_dp,0.069138_dp,0.047426_dp,0.930862_dp],[2,2])
    model%log_likelihood=-100.0_dp
    sigma=1.0_dp
    call rsdc_forecast_path(model,sim%observations,sigma,fc,x=x,ok=ok)
    call check(ok,'forecast path')
    call check(all(shape(fc%predicted_correlations)==[t,1]),'forecast dimensions')
    call check(maxval(abs(sum(fc%regime_probabilities,dim=2)-1.0_dp))<1.0e-12_dp,'forecast probabilities')
    call rsdc_forecast_ahead(model,sim%observations,5,reg,pc,x=x,ok=ok)
    call check(ok .and. all(shape(reg)==[5,n]),'ahead forecast')
    call rsdc_viterbi_path(model,sim%observations,path,x=x,ok=ok)
    call check(ok,'viterbi')
    correct=count(path==sim%states)
    call check(correct>t/2,'viterbi recovery')
    print '(a)', 'test_simulation_forecast: PASS'
contains
    subroutine check(condition,message)
        logical,intent(in)::condition
        character(len=*),intent(in)::message
        if(.not.condition) error stop message
    end subroutine check
end program test_simulation_forecast
