program test_nfcp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use nfcp
  implicit none

  integer :: failures
  failures = 0
  call test_parameters(failures)
  call test_core_formulas(failures)
  call test_forecast_and_filter(failures)
  call test_options(failures)
  call test_stitching(failures)
  call test_mle(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  write(*,'(a)') 'All NFCP tests passed.'

contains

  subroutine test_parameters(failures)
    integer, intent(inout) :: failures
    type(nfcp_model_t) :: model, recovered
    type(nfcp_parameterization_t) :: spec
    real(dp), allocatable :: theta(:), x0_out(:), lower(:), upper(:)
    real(dp) :: x0(2)
    character(len=32), allocatable :: names(:)
    integer :: status

    call make_two_factor(model)
    spec = nfcp_parameterization_t(2,0,1,.true.,.true.)
    x0 = [log(50.0_dp), -0.2_dp]
    names = spec%names()
    call assert_true(size(names) == 10, 'parameter count', failures)
    call assert_true(trim(names(1)) == 'x_0_1' .and. trim(names(4)) == 'mu_rn', &
                     'parameter names', failures)
    call pack_parameters(model,spec,x0,theta,status)
    call assert_true(status == nfcp_ok, 'pack parameters', failures)
    call unpack_parameters(theta,spec,recovered,x0_out,status)
    call assert_true(status == nfcp_ok, 'unpack parameters', failures)
    call assert_close(maxval(abs(recovered%sigma-model%sigma)),0.0_dp,1.0e-14_dp, &
                      'parameter round trip',failures)
    call assert_close(maxval(abs(x0_out-x0)),0.0_dp,1.0e-14_dp,'initial state round trip',failures)
    call default_parameter_bounds(spec,lower,upper)
    call assert_true(all(lower < upper), 'default bounds', failures)
  end subroutine test_parameters

  subroutine test_core_formulas(failures)
    integer, intent(inout) :: failures
    type(nfcp_model_t) :: model
    real(dp) :: covariance(2,2), value
    call make_two_factor(model)
    value = nfcp_a_t(model,0.75_dp)
    call assert_close(value,-0.017236185911698815_dp,2.0e-14_dp,'A(T)',failures)
    call nfcp_covariance(model,0.1_dp,covariance)
    call assert_close(covariance(1,1),0.004_dp,2.0e-14_dp,'covariance 11',failures)
    call assert_close(covariance(1,2),-0.001130795632828425_dp,2.0e-14_dp,'covariance 12',failures)
    call assert_close(covariance(2,2),0.008001455210004244_dp,2.0e-14_dp,'covariance 22',failures)
    call assert_close(normal_quantile(0.975_dp),1.959963986120195_dp,2.0e-8_dp, &
                      'normal quantile',failures)
  end subroutine test_core_formulas

  subroutine test_forecast_and_filter(failures)
    integer, intent(inout) :: failures
    type(nfcp_model_t) :: model
    type(nfcp_futures_simulation_result_t) :: sim
    type(nfcp_filter_result_t) :: filter
    real(dp), allocatable :: forecast(:,:), theoretical(:), empirical(:)
    real(dp), allocatable :: ttm(:,:), log_futures(:,:)
    real(dp) :: x0(2), nanv
    integer :: t, status

    call make_two_factor(model)
    x0 = [log(50.0_dp),-0.2_dp]
    call futures_price_forecast(model,x0,0.0_dp,[0.25_dp,1.0_dp],forecast,status=status)
    call assert_true(status == nfcp_ok .and. all(forecast > 0.0_dp), 'futures forecast',failures)
    allocate(ttm(80,3))
    do t=1,80
      ttm(t,:)=[0.25_dp,0.75_dp,1.5_dp]
    end do
    call futures_price_simulate(model,x0,1.0_dp/52.0_dp,ttm,sim,seed=2468)
    call assert_true(sim%status == nfcp_ok, 'futures simulation', failures)
    log_futures=log(sim%futures_prices)
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    log_futures(10,2)=nanv
    call nfcp_kalman_filter(model,log_futures,1.0_dp/52.0_dp,ttm,filter,initial_state=x0,n_parameters=8)
    call assert_true(filter%status == nfcp_ok .and. ieee_is_finite(filter%log_likelihood), &
                     'Kalman filter', failures)
    call assert_true(all(ieee_is_finite(filter%final_state)), 'filtered state finite', failures)
    call tsfit_volatility(model,sim%futures_prices,[0.25_dp,0.75_dp,1.5_dp],1.0_dp/52.0_dp, &
                          theoretical,empirical,status)
    call assert_true(status == nfcp_ok .and. all(theoretical > 0.0_dp), &
                     'volatility term structure',failures)
  end subroutine test_forecast_and_filter

  subroutine test_options(failures)
    integer, intent(inout) :: failures
    type(nfcp_model_t) :: model
    type(nfcp_option_result_t) :: european, american
    real(dp) :: x0(1)

    call initialize_model(model,1,.true.,0,0)
    model%mu=0.02_dp;model%mu_rn=0.01_dp;model%sigma(1)=0.2_dp
    model%lambda(1)=model%mu-model%mu_rn
    x0=[log(40.0_dp)]
    call european_option_value(model,x0,1.0_dp,0.5_dp,40.0_dp,0.03_dp,.false.,european)
    call assert_true(european%status == nfcp_ok .and. european%value > 0.0_dp, &
                     'European option',failures)
    call american_option_value(model,x0,1.0_dp,0.5_dp,40.0_dp,0.03_dp,.false.,1200, &
                               1.0_dp/12.0_dp,american,degree=2,seed=1001)
    call assert_true(american%status == nfcp_ok .and. american%value > 0.0_dp, &
                     'American option',failures)
    call assert_true(american%value + 4.0_dp*american%standard_error >= european%value, &
                     'American value versus European value',failures)
  end subroutine test_options

  subroutine test_stitching(failures)
    integer, intent(inout) :: failures
    real(dp) :: futures(3,4), maturity(3,4), nanv
    real(dp), allocatable :: stitched(:,:), stitched_maturity(:,:)
    integer, allocatable :: selected(:,:)
    integer :: status

    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    futures=reshape([10.0_dp,11.0_dp,12.0_dp, 20.0_dp,21.0_dp,22.0_dp, &
                     30.0_dp,31.0_dp,32.0_dp, 40.0_dp,41.0_dp,42.0_dp],[3,4])
    futures(1,1)=nanv
    call stitch_contract_numbers(futures,[1,2],stitched,selected,status)
    call assert_true(status == nfcp_ok .and. selected(1,1)==2 .and. selected(2,1)==1, &
                     'contract-number stitching',failures)
    maturity=reshape([0.10_dp,0.08_dp,0.06_dp, 0.35_dp,0.33_dp,0.31_dp, &
                      0.60_dp,0.58_dp,0.56_dp, 0.85_dp,0.83_dp,0.81_dp],[3,4])
    call stitch_by_maturity(futures,maturity,[0.35_dp],0.10_dp,stitched,stitched_maturity, &
                            selected,status)
    call assert_true(status == nfcp_ok .and. size(stitched,2)==1, 'maturity stitching',failures)
  end subroutine test_stitching

  subroutine test_mle(failures)
    integer, intent(inout) :: failures
    type(nfcp_model_t) :: model
    type(nfcp_parameterization_t) :: spec
    type(nfcp_futures_simulation_result_t) :: sim
    type(nfcp_mle_control_t) :: control
    type(nfcp_mle_result_t) :: fit
    real(dp), allocatable :: ttm(:,:), theta(:), lower(:), upper(:)
    real(dp) :: x0(1)
    integer :: t,status

    call initialize_model(model,1,.true.,1,0)
    model%mu=0.025_dp;model%mu_rn=0.01_dp;model%sigma(1)=0.16_dp
    model%lambda(1)=model%mu-model%mu_rn;model%measurement_error(1)=0.008_dp
    x0=[log(45.0_dp)]
    allocate(ttm(35,2))
    do t=1,35;ttm(t,:)=[0.25_dp,0.75_dp];end do
    call futures_price_simulate(model,x0,1.0_dp/52.0_dp,ttm,sim,seed=777)
    spec=nfcp_parameterization_t(1,0,1,.true.,.false.)
    call pack_parameters(model,spec,theta=theta,status=status)
    allocate(lower(4),upper(4))
    lower=[-0.10_dp,-0.10_dp,0.05_dp,0.001_dp]
    upper=[ 0.10_dp, 0.10_dp,0.40_dp,0.05_dp]
    theta=theta+[0.005_dp,-0.004_dp,0.02_dp,0.003_dp]
    control=nfcp_mle_control_t(population_size=8,generations=1,wait_generations=1, &
      local_max_evaluations=300,seed=99,local_refinement=.true.,trace=0)
    call nfcp_fit_mle(spec,log(sim%futures_prices),1.0_dp/52.0_dp,ttm,lower,upper,fit, &
                      initial_state=x0,starting_values=theta,control=control)
    call assert_true(fit%status == nfcp_ok .and. ieee_is_finite(fit%log_likelihood), &
                     'maximum likelihood estimation',failures)
    call assert_true(all(fit%parameters>=lower) .and. all(fit%parameters<=upper), &
                     'MLE bounds',failures)
  end subroutine test_mle

  subroutine make_two_factor(model)
    type(nfcp_model_t),intent(out)::model
    call initialize_model(model,2,.true.,1,0)
    model%mu=0.03_dp;model%mu_rn=0.01_dp;model%equilibrium=0.0_dp
    model%lambda=[model%mu-model%mu_rn,0.1_dp]
    model%kappa=[0.0_dp,1.2_dp]
    model%sigma=[0.2_dp,0.3_dp]
    model%rho=reshape([1.0_dp,-0.2_dp,-0.2_dp,1.0_dp],[2,2])
    model%measurement_error=[0.01_dp]
  end subroutine make_two_factor

  subroutine assert_true(condition,name,failures)
    logical,intent(in)::condition
    character(len=*),intent(in)::name
    integer,intent(inout)::failures
    if(.not.condition)then
      failures=failures+1;write(*,'(a)')'FAIL: '//trim(name)
    end if
  end subroutine assert_true

  subroutine assert_close(actual,expected,tolerance,name,failures)
    real(dp),intent(in)::actual,expected,tolerance
    character(len=*),intent(in)::name
    integer,intent(inout)::failures
    if(abs(actual-expected)>tolerance)then
      failures=failures+1
      write(*,'(a,2es24.15)')'FAIL: '//trim(name)//' actual/expected: ',actual,expected
    end if
  end subroutine assert_close

end program test_nfcp
