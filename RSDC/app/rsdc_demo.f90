! SPDX-License-Identifier: GPL-3.0-only
program rsdc_demo
    use rsdc, only: dp, rsdc_model, rsdc_control, rsdc_simulation_result
    use rsdc, only: rsdc_simulate_fixed, rsdc_estimate, rsdc_forecast_ahead
    use rsdc, only: rsdc_nox
    implicit none
    integer, parameter :: nobs = 220, k = 2, n = 2
    real(dp) :: p(n,n), mu(n,k), covariance(k,k,n)
    real(dp), allocatable :: regime_probs(:,:), correlations(:,:)
    type(rsdc_simulation_result) :: sim
    type(rsdc_model) :: fit
    type(rsdc_control) :: control
    logical :: ok

    p = reshape([0.96_dp, 0.08_dp, 0.04_dp, 0.92_dp], [2,2])
    mu = 0.0_dp
    covariance(:,:,1) = reshape([1.0_dp,0.10_dp,0.10_dp,1.0_dp],[2,2])
    covariance(:,:,2) = reshape([1.0_dp,0.78_dp,0.78_dp,1.0_dp],[2,2])
    call rsdc_simulate_fixed(nobs,p,mu,covariance,sim,seed=2026,ok=ok)
    if (.not. ok) error stop 'simulation failed'

    control%population_size = 40
    control%max_global_iterations = 100
    control%max_local_iterations = 70
    control%seed = 17
    call rsdc_estimate(rsdc_nox,sim%observations,n,fit,control=control,ok=ok)
    if (.not. ok) error stop 'estimation failed'

    call rsdc_forecast_ahead(fit,sim%observations,5,regime_probs,correlations,ok=ok)
    if (.not. ok) error stop 'forecast failed'

    print '(a,f12.4)', 'log likelihood: ', fit%log_likelihood
    print '(a,2f10.4)', 'regime correlations: ', fit%correlations(:,1)
    print '(a)', 'five-step predicted correlations:'
    print '(5f10.4)', correlations(:,1)
end program rsdc_demo
