! SPDX-License-Identifier: MIT
module bekks
  use bekks_kinds, only: dp
  use bekks_types
  use bekks_model, only: parameter_count, unpack_parameters, pack_parameters, &
    indicator_function, expected_indicator_value, valid_parameters, initial_parameters, &
    random_initial_parameters, log_likelihood, filter_bekk, simulate_core => simulate_bekk, &
    unconditional_covariance, covariance_to_volatility
  use bekks_estimation, only: score_core => score_bekk, hessian_core => hessian_bekk, &
    bhhh_fit, bekk_fit, qml_covariance, rmse_parameters
  use bekks_forecast, only: forecast_bekk, virf_bekk
  use bekks_risk, only: var_bekk_fit, var_bekk_forecast, coverage_tests, &
    backtest_forecasts, rolling_backtest
  use bekks_diagnostics, only: portmanteau_test, bekk_mc_eval
  use bekks_matrix
  use bekks_linalg, only: general_inverse, symmetric_sqrt
  use bekks_rng, only: rng_state, rng_seed
  implicit none
  private

  public :: dp, rng_state, rng_seed
  public :: bekk_full, bekk_diagonal, bekk_scalar
  public :: bekk_ok, bekk_invalid_input, bekk_invalid_parameters, &
    bekk_linalg_failure, bekk_no_convergence
  public :: bekk_spec_type, bekk_parameters, bekk_fit_result, bekk_filter_result
  public :: bekk_forecast_result, bekk_var_result, bekk_backtest_result
  public :: bekk_virf_result, bekk_portmanteau_result
  public :: bekk_spec, bekk_fit, bhhh_fit, qml_covariance, forecast_bekk, virf_bekk
  public :: var_bekk_fit, var_bekk_forecast, rolling_backtest, backtest_forecasts
  public :: coverage_tests, portmanteau_test, bekk_mc_eval
  public :: parameter_count, unpack_parameters, pack_parameters
  public :: indicator_function, expected_indicator_value
  public :: initial_parameters, random_initial_parameters
  public :: grid_search_bekk, grid_search_asymmetric_bekk
  public :: grid_search_dbekk, grid_search_asymmetric_dbekk
  public :: grid_search_sbekk, grid_search_asymmetric_sbekk
  public :: random_grid_search_bekk, random_grid_search_asymmetric_bekk
  public :: random_grid_search_dbekk, random_grid_search_asymmetric_dbekk
  public :: random_grid_search_sbekk, random_grid_search_asymmetric_sbekk
  public :: unconditional_covariance, covariance_to_volatility
  public :: elimination_mat, commutation_mat, duplication_mat, diag_selection_mat
  public :: cut_mat_symmetric, cut_mat_asymmetric, vech_lower, unvech_lower
  public :: y_lag_cr, extract_csd, general_inverse, symmetric_sqrt
  public :: loglike_bekk, loglike_asymm_bekk, loglike_dbekk, loglike_asymm_dbekk
  public :: loglike_sbekk, loglike_asymm_sbekk
  public :: score_bekk, score_asymm_bekk, score_dbekk, score_asymm_dbekk
  public :: score_sbekk, score_asymm_sbekk
  public :: hesse_bekk, hesse_asymm_bekk, hesse_dbekk, hesse_asymm_dbekk
  public :: hesse_sbekk, hesse_asymm_sbekk
  public :: sigma_bekk, sigma_bekk_asymm, sigma_dbekk, sigma_dbekk_asymm
  public :: sigma_sbekk, sigma_sbekk_asymm
  public :: valid_bekk, valid_asymm_bekk, valid_dbekk, valid_asymm_dbekk
  public :: valid_sbekk, valid_asymm_sbekk
  public :: simulate_bekk_model
  public :: simulate_bekk_full, simulate_bekk_asymm, simulate_dbekk
  public :: simulate_dbekk_asymm, simulate_sbekk, simulate_sbekk_asymm
  public :: rmse_parameters

contains

  function bekk_spec(model_type,asymmetric,signs,initial_theta) result(spec)
    integer, intent(in), optional :: model_type
    logical, intent(in), optional :: asymmetric
    real(dp), intent(in), optional :: signs(:),initial_theta(:)
    type(bekk_spec_type) :: spec
    spec%model_type=bekk_full
    spec%asymmetric=.false.
    if(present(model_type))spec%model_type=model_type
    if(present(asymmetric))spec%asymmetric=asymmetric
    if(present(signs))then
      allocate(spec%signs(size(signs)))
      spec%signs=signs
    end if
    if(present(initial_theta))then
      allocate(spec%initial_theta(size(initial_theta)))
      spec%initial_theta=initial_theta
    end if
  end function bekk_spec

  real(dp) function loglike_bekk(theta,data) result(v)
    real(dp), intent(in) :: theta(:),data(:,:)
    v=log_likelihood(theta,data,bekk_full,.false.)
  end function loglike_bekk
  real(dp) function loglike_asymm_bekk(theta,data,signs) result(v)
    real(dp), intent(in) :: theta(:),data(:,:),signs(:)
    v=log_likelihood(theta,data,bekk_full,.true.,signs)
  end function loglike_asymm_bekk
  real(dp) function loglike_dbekk(theta,data) result(v)
    real(dp), intent(in) :: theta(:),data(:,:)
    v=log_likelihood(theta,data,bekk_diagonal,.false.)
  end function loglike_dbekk
  real(dp) function loglike_asymm_dbekk(theta,data,signs) result(v)
    real(dp), intent(in) :: theta(:),data(:,:),signs(:)
    v=log_likelihood(theta,data,bekk_diagonal,.true.,signs)
  end function loglike_asymm_dbekk
  real(dp) function loglike_sbekk(theta,data) result(v)
    real(dp), intent(in) :: theta(:),data(:,:)
    v=log_likelihood(theta,data,bekk_scalar,.false.)
  end function loglike_sbekk
  real(dp) function loglike_asymm_sbekk(theta,data,signs) result(v)
    real(dp), intent(in) :: theta(:),data(:,:),signs(:)
    v=log_likelihood(theta,data,bekk_scalar,.true.,signs)
  end function loglike_asymm_sbekk

  subroutine score_bekk(theta,data,score,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_full,.false.,score=score,status=status)
  end subroutine score_bekk
  subroutine score_asymm_bekk(theta,data,signs,score,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_full,.true.,signs,score,status)
  end subroutine score_asymm_bekk
  subroutine score_dbekk(theta,data,score,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_diagonal,.false.,score=score,status=status)
  end subroutine score_dbekk
  subroutine score_asymm_dbekk(theta,data,signs,score,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_diagonal,.true.,signs,score,status)
  end subroutine score_asymm_dbekk
  subroutine score_sbekk(theta,data,score,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_scalar,.false.,score=score,status=status)
  end subroutine score_sbekk
  subroutine score_asymm_sbekk(theta,data,signs,score,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::score(:,:)
    integer,intent(out)::status
    call score_core(theta,data,bekk_scalar,.true.,signs,score,status)
  end subroutine score_asymm_sbekk

  subroutine hesse_bekk(theta,data,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_full,.false.,hessian=hessian,status=status)
  end subroutine hesse_bekk
  subroutine hesse_asymm_bekk(theta,data,signs,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_full,.true.,signs,hessian,status)
  end subroutine hesse_asymm_bekk
  subroutine hesse_dbekk(theta,data,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_diagonal,.false.,hessian=hessian,status=status)
  end subroutine hesse_dbekk
  subroutine hesse_asymm_dbekk(theta,data,signs,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_diagonal,.true.,signs,hessian,status)
  end subroutine hesse_asymm_dbekk
  subroutine hesse_sbekk(theta,data,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_scalar,.false.,hessian=hessian,status=status)
  end subroutine hesse_sbekk
  subroutine hesse_asymm_sbekk(theta,data,signs,hessian,status)
    real(dp),intent(in)::theta(:),data(:,:),signs(:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::status
    call hessian_core(theta,data,bekk_scalar,.true.,signs,hessian,status)
  end subroutine hesse_asymm_sbekk

  subroutine sigma_generic(theta,data,model_type,asymmetric,signs,result)
    real(dp),intent(in)::theta(:),data(:,:)
    integer,intent(in)::model_type
    logical,intent(in)::asymmetric
    real(dp),intent(in),optional::signs(:)
    type(bekk_filter_result),intent(out)::result
    call filter_bekk(theta,data,model_type,asymmetric,signs,result%h,result%residuals,result%status)
  end subroutine sigma_generic
  subroutine sigma_bekk(theta,data,result)
    real(dp),intent(in)::theta(:),data(:,:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_full,.false.,result=result)
  end subroutine sigma_bekk
  subroutine sigma_bekk_asymm(theta,data,signs,result)
    real(dp),intent(in)::theta(:),data(:,:),signs(:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_full,.true.,signs,result)
  end subroutine sigma_bekk_asymm
  subroutine sigma_dbekk(theta,data,result)
    real(dp),intent(in)::theta(:),data(:,:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_diagonal,.false.,result=result)
  end subroutine sigma_dbekk
  subroutine sigma_dbekk_asymm(theta,data,signs,result)
    real(dp),intent(in)::theta(:),data(:,:),signs(:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_diagonal,.true.,signs,result)
  end subroutine sigma_dbekk_asymm
  subroutine sigma_sbekk(theta,data,result)
    real(dp),intent(in)::theta(:),data(:,:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_scalar,.false.,result=result)
  end subroutine sigma_sbekk
  subroutine sigma_sbekk_asymm(theta,data,signs,result)
    real(dp),intent(in)::theta(:),data(:,:),signs(:);type(bekk_filter_result),intent(out)::result
    call sigma_generic(theta,data,bekk_scalar,.true.,signs,result)
  end subroutine sigma_sbekk_asymm

  logical function valid_generic(theta,n,model_type,asymmetric,e) result(ok)
    real(dp),intent(in)::theta(:),e
    integer,intent(in)::n,model_type
    logical,intent(in)::asymmetric
    type(bekk_parameters)::p
    integer::st
    call unpack_parameters(theta,n,model_type,asymmetric,p,st)
    ok=st==bekk_ok
    if(ok)ok=valid_parameters(p,e)
  end function valid_generic
  logical function valid_bekk(theta,n) result(ok)
    real(dp),intent(in)::theta(:);integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_full,.false.,0.0_dp)
  end function valid_bekk
  logical function valid_asymm_bekk(theta,n,expected_indicator) result(ok)
    real(dp),intent(in)::theta(:),expected_indicator;integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_full,.true.,expected_indicator)
  end function valid_asymm_bekk
  logical function valid_dbekk(theta,n) result(ok)
    real(dp),intent(in)::theta(:);integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_diagonal,.false.,0.0_dp)
  end function valid_dbekk
  logical function valid_asymm_dbekk(theta,n,expected_indicator) result(ok)
    real(dp),intent(in)::theta(:),expected_indicator;integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_diagonal,.true.,expected_indicator)
  end function valid_asymm_dbekk
  logical function valid_sbekk(theta,n) result(ok)
    real(dp),intent(in)::theta(:);integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_scalar,.false.,0.0_dp)
  end function valid_sbekk
  logical function valid_asymm_sbekk(theta,n,expected_indicator) result(ok)
    real(dp),intent(in)::theta(:),expected_indicator;integer,intent(in)::n
    ok=valid_generic(theta,n,bekk_scalar,.true.,expected_indicator)
  end function valid_asymm_sbekk

  subroutine simulate_generic(theta,nobs,n,model_type,asymmetric,state,data,h,status,signs,expected_indicator)
    real(dp),intent(in)::theta(:)
    integer,intent(in)::nobs,n,model_type
    logical,intent(in)::asymmetric
    type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::data(:,:),h(:,:,:)
    integer,intent(out)::status
    real(dp),intent(in),optional::signs(:),expected_indicator
    call simulate_core(theta,nobs,n,model_type,asymmetric,signs,expected_indicator,state,data,h,status)
  end subroutine simulate_generic
  subroutine simulate_bekk_full(theta,nobs,n,state,data,h,status)
    real(dp),intent(in)::theta(:);integer,intent(in)::nobs,n;type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_full,.false.,state,data,h,status)
  end subroutine simulate_bekk_full
  subroutine simulate_bekk_asymm(theta,nobs,n,state,signs,expected_indicator,data,h,status)
    real(dp),intent(in)::theta(:),signs(:),expected_indicator;integer,intent(in)::nobs,n
    type(rng_state),intent(inout)::state;real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_full,.true.,state,data,h,status,signs,expected_indicator)
  end subroutine simulate_bekk_asymm
  subroutine simulate_dbekk(theta,nobs,n,state,data,h,status)
    real(dp),intent(in)::theta(:);integer,intent(in)::nobs,n;type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_diagonal,.false.,state,data,h,status)
  end subroutine simulate_dbekk
  subroutine simulate_dbekk_asymm(theta,nobs,n,state,signs,expected_indicator,data,h,status)
    real(dp),intent(in)::theta(:),signs(:),expected_indicator;integer,intent(in)::nobs,n
    type(rng_state),intent(inout)::state;real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_diagonal,.true.,state,data,h,status,signs,expected_indicator)
  end subroutine simulate_dbekk_asymm
  subroutine simulate_sbekk(theta,nobs,n,state,data,h,status)
    real(dp),intent(in)::theta(:);integer,intent(in)::nobs,n;type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_scalar,.false.,state,data,h,status)
  end subroutine simulate_sbekk
  subroutine simulate_sbekk_asymm(theta,nobs,n,state,signs,expected_indicator,data,h,status)
    real(dp),intent(in)::theta(:),signs(:),expected_indicator;integer,intent(in)::nobs,n
    type(rng_state),intent(inout)::state;real(dp),allocatable,intent(out)::data(:,:),h(:,:,:);integer,intent(out)::status
    call simulate_generic(theta,nobs,n,bekk_scalar,.true.,state,data,h,status,signs,expected_indicator)
  end subroutine simulate_sbekk_asymm

  subroutine simulate_bekk_model(theta,nobs,n,model_type,asymmetric,state,data,h,status, &
      signs,expected_indicator,innovations)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs,n,model_type
    logical, intent(in) :: asymmetric
    type(rng_state), intent(inout) :: state
    real(dp), allocatable, intent(out) :: data(:,:),h(:,:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: signs(:),expected_indicator,innovations(:,:)
    call simulate_core(theta,nobs,n,model_type,asymmetric,signs,expected_indicator, &
      state,data,h,status,innovations)
  end subroutine simulate_bekk_model

  subroutine grid_search_generic(data,model_type,asymmetric,theta,likelihood,status,signs)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), allocatable, intent(out) :: theta(:)
    real(dp), intent(out) :: likelihood
    integer, intent(out) :: status
    real(dp), intent(in), optional :: signs(:)
    type(bekk_parameters) :: par
    real(dp), allocatable :: s(:)
    call initial_parameters(data,model_type,asymmetric,par,status)
    if(status/=bekk_ok)then
      likelihood=-huge(1.0_dp)
      return
    end if
    theta=pack_parameters(par)
    if(asymmetric)then
      allocate(s(size(data,2)));s=-1.0_dp
      if(present(signs))s=signs
      likelihood=log_likelihood(theta,data,model_type,.true.,s,status)
    else
      likelihood=log_likelihood(theta,data,model_type,.false.,status=status)
    end if
  end subroutine grid_search_generic

  subroutine random_grid_search_generic(data,model_type,asymmetric,state,theta,likelihood,status,signs,n_trials)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    type(rng_state), intent(inout) :: state
    real(dp), allocatable, intent(out) :: theta(:)
    real(dp), intent(out) :: likelihood
    integer, intent(out) :: status
    real(dp), intent(in), optional :: signs(:)
    integer, intent(in), optional :: n_trials
    type(bekk_parameters) :: par
    real(dp), allocatable :: s(:)
    call random_initial_parameters(data,model_type,asymmetric,state,par,status,n_trials)
    if(status/=bekk_ok)then
      likelihood=-huge(1.0_dp)
      return
    end if
    theta=pack_parameters(par)
    if(asymmetric)then
      allocate(s(size(data,2)));s=-1.0_dp
      if(present(signs))s=signs
      likelihood=log_likelihood(theta,data,model_type,.true.,s,status)
    else
      likelihood=log_likelihood(theta,data,model_type,.false.,status=status)
    end if
  end subroutine random_grid_search_generic

  subroutine grid_search_bekk(data,theta,likelihood,status)
    real(dp),intent(in)::data(:,:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_full,.false.,theta,likelihood,status)
  end subroutine grid_search_bekk

  subroutine grid_search_asymmetric_bekk(data,signs,theta,likelihood,status)
    real(dp),intent(in)::data(:,:),signs(:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_full,.true.,theta,likelihood,status,signs)
  end subroutine grid_search_asymmetric_bekk

  subroutine grid_search_dbekk(data,theta,likelihood,status)
    real(dp),intent(in)::data(:,:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_diagonal,.false.,theta,likelihood,status)
  end subroutine grid_search_dbekk

  subroutine grid_search_asymmetric_dbekk(data,signs,theta,likelihood,status)
    real(dp),intent(in)::data(:,:),signs(:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_diagonal,.true.,theta,likelihood,status,signs)
  end subroutine grid_search_asymmetric_dbekk

  subroutine grid_search_sbekk(data,theta,likelihood,status)
    real(dp),intent(in)::data(:,:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_scalar,.false.,theta,likelihood,status)
  end subroutine grid_search_sbekk

  subroutine grid_search_asymmetric_sbekk(data,signs,theta,likelihood,status)
    real(dp),intent(in)::data(:,:),signs(:);real(dp),allocatable,intent(out)::theta(:)
    real(dp),intent(out)::likelihood;integer,intent(out)::status
    call grid_search_generic(data,bekk_scalar,.true.,theta,likelihood,status,signs)
  end subroutine grid_search_asymmetric_sbekk

  subroutine random_grid_search_bekk(data,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_full,.false.,state,theta,likelihood,status,n_trials=n_trials)
  end subroutine random_grid_search_bekk

  subroutine random_grid_search_asymmetric_bekk(data,signs,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:),signs(:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_full,.true.,state,theta,likelihood,status,signs,n_trials)
  end subroutine random_grid_search_asymmetric_bekk

  subroutine random_grid_search_dbekk(data,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_diagonal,.false.,state,theta,likelihood,status,n_trials=n_trials)
  end subroutine random_grid_search_dbekk

  subroutine random_grid_search_asymmetric_dbekk(data,signs,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:),signs(:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_diagonal,.true.,state,theta,likelihood,status,signs,n_trials)
  end subroutine random_grid_search_asymmetric_dbekk

  subroutine random_grid_search_sbekk(data,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_scalar,.false.,state,theta,likelihood,status,n_trials=n_trials)
  end subroutine random_grid_search_sbekk

  subroutine random_grid_search_asymmetric_sbekk(data,signs,state,theta,likelihood,status,n_trials)
    real(dp),intent(in)::data(:,:),signs(:);type(rng_state),intent(inout)::state
    real(dp),allocatable,intent(out)::theta(:);real(dp),intent(out)::likelihood
    integer,intent(out)::status;integer,intent(in),optional::n_trials
    call random_grid_search_generic(data,bekk_scalar,.true.,state,theta,likelihood,status,signs,n_trials)
  end subroutine random_grid_search_asymmetric_sbekk

end module bekks
