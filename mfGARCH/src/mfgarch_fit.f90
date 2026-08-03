! SPDX-License-Identifier: MIT
module mfgarch_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use mfgarch_kinds, only : dp
  use mfgarch_components, only : beta_weights, likelihood_contributions, log_likelihood, &
    model_parameter_count, model_to_raw, raw_to_model, model_parameters, forecast_tau, &
    variance_ratio, valid_model
  use mfgarch_math, only : finite_value, invert_matrix, sample_mean
  use mfgarch_optimization, only : optimize_bfgs, optimize_nelder_mead
  use mfgarch_status, only : mfgarch_success, mfgarch_invalid_argument, &
    mfgarch_dimension_error, mfgarch_not_converged, mfgarch_singular_matrix, mfgarch_message
  use mfgarch_types, only : mfgarch_model, mfgarch_fit_control, mfgarch_fit_result
  implicit none
  private

  public :: fit_mfgarch, compute_mfgarch_inference

contains

  subroutine fit_mfgarch(returns, period, start_model, result, status, covariate, &
      period_two, covariate_two, variance_period, control)
    real(dp), intent(in) :: returns(:)
    integer, intent(in) :: period(:)
    type(mfgarch_model), intent(in) :: start_model
    type(mfgarch_fit_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: covariate(:), covariate_two(:)
    integer, intent(in), optional :: period_two(:), variance_period(:)
    type(mfgarch_fit_control), intent(in), optional :: control
    type(mfgarch_fit_control) :: ctrl
    type(mfgarch_model) :: fitted_model
    real(dp), allocatable :: raw_start(:), raw_solution(:), nm_solution(:), refined_solution(:)
    real(dp), allocatable :: contributions(:)
    real(dp) :: objective_value, nm_value, refined_value
    integer :: optimizer_status, nm_status, refined_status
    integer :: iterations, evaluations, nm_iterations, nm_evaluations
    integer :: refined_iterations, refined_evaluations, component_status, nvalid, i
    character(len=:), allocatable :: status_text

    ctrl = mfgarch_fit_control()
    if (present(control)) ctrl = control
    result%model = start_model
    status = mfgarch_success

    if (size(returns) /= size(period) .or. size(returns) < 5 .or. any(period < 1)) then
      status = mfgarch_dimension_error
      result%status = status
      result%message = mfgarch_message(status)
      return
    end if
    if (.not. valid_model(start_model)) then
      status = mfgarch_invalid_argument
      result%status = status
      result%message = mfgarch_message(status)
      return
    end if
    if (start_model%k > 0 .and. .not. present(covariate)) then
      status = mfgarch_invalid_argument
      result%status = status
      result%message = 'covariate is required when K > 0'
      return
    end if
    if (start_model%has_second .and. (.not. present(period_two) .or. &
        .not. present(covariate_two))) then
      status = mfgarch_invalid_argument
      result%status = status
      result%message = 'second period index and covariate are required'
      return
    end if

    call model_to_raw(start_model, raw_start, status)
    if (status /= mfgarch_success) then
      result%status = status
      result%message = mfgarch_message(status)
      return
    end if

    select case (trim(adjustl(ctrl%method)))
    case ('nelder-mead', 'nm')
      call optimize_nelder_mead(objective_function, raw_start, ctrl%initial_simplex_step, &
        ctrl%tolerance, ctrl%max_iterations, ctrl%max_function_evaluations, ctrl%trace, &
        raw_solution, objective_value, optimizer_status, iterations, evaluations)
    case default
      call optimize_bfgs(objective_function, raw_start, ctrl%tolerance, ctrl%gradient_step, &
        ctrl%max_iterations, ctrl%max_function_evaluations, ctrl%trace, raw_solution, &
        objective_value, optimizer_status, iterations, evaluations)
      if (ctrl%multi_start) then
        call optimize_nelder_mead(objective_function, raw_start, ctrl%initial_simplex_step, &
          max(10.0_dp*ctrl%tolerance,1.0e-6_dp), max(100,ctrl%max_iterations/2), &
          ctrl%max_function_evaluations, ctrl%trace, nm_solution, nm_value, nm_status, &
          nm_iterations, nm_evaluations)
        evaluations = evaluations + nm_evaluations
        if (finite_value(nm_value) .and. (nm_value < objective_value .or. &
            .not. finite_value(objective_value))) then
          raw_solution = nm_solution
          objective_value = nm_value
          optimizer_status = nm_status
          iterations = nm_iterations
        end if
        call optimize_bfgs(objective_function, raw_solution, ctrl%tolerance, ctrl%gradient_step, &
          ctrl%max_iterations, ctrl%max_function_evaluations, ctrl%trace, refined_solution, &
          refined_value, refined_status, refined_iterations, refined_evaluations)
        evaluations = evaluations + refined_evaluations
        if (finite_value(refined_value) .and. refined_value < objective_value) then
          raw_solution = refined_solution
          objective_value = refined_value
          optimizer_status = refined_status
          iterations = iterations + refined_iterations
        end if
      end if
    end select

    call raw_to_model(raw_solution, start_model, fitted_model, status)
    if (status /= mfgarch_success) then
      result%status = status
      result%message = mfgarch_message(status)
      return
    end if
    result%model = fitted_model
    result%log_likelihood = -objective_value
    result%iterations = iterations
    result%function_evaluations = evaluations
    result%status = optimizer_status
    status = optimizer_status
    status_text = mfgarch_message(optimizer_status)
    result%message = status_text

    if (fitted_model%k == 0) then
      call likelihood_contributions(fitted_model, returns, period, contributions, result%tau, &
        result%g, result%residuals, component_status)
    else if (fitted_model%has_second) then
      call likelihood_contributions(fitted_model, returns, period, contributions, result%tau, &
        result%g, result%residuals, component_status, covariate, period_two, covariate_two)
    else
      call likelihood_contributions(fitted_model, returns, period, contributions, result%tau, &
        result%g, result%residuals, component_status, covariate)
    end if
    if (component_status /= mfgarch_success) then
      status = component_status
      result%status = status
      result%message = mfgarch_message(status)
      return
    end if
    nvalid = count([(finite_value(contributions(i)), i=1,size(contributions))])
    result%bic = log(real(max(1,nvalid),dp)) * real(model_parameter_count(fitted_model),dp) - &
      2.0_dp * result%log_likelihood

    if (fitted_model%k == 0) then
      result%tau_forecast = exp(fitted_model%m)
      result%variance_ratio = 0.0_dp
    else
      if (fitted_model%has_second) then
        result%tau_forecast = forecast_tau(fitted_model, covariate, component_status, covariate_two)
      else
        result%tau_forecast = forecast_tau(fitted_model, covariate, component_status)
      end if
      if (component_status /= mfgarch_success) result%tau_forecast = 0.0_dp
      if (present(variance_period)) then
        result%variance_ratio = variance_ratio(result%tau, result%g, variance_period, component_status)
      else
        result%variance_ratio = variance_ratio(result%tau, result%g, period, component_status)
      end if
    end if

    call beta_weights(fitted_model%k, fitted_model%w1, fitted_model%w2, result%weights, component_status)
    if (fitted_model%has_second) then
      call beta_weights(fitted_model%k_two, fitted_model%w1_two, fitted_model%w2_two, &
        result%weights_two, component_status)
    else
      allocate(result%weights_two(0))
    end if

    if (ctrl%compute_inference) then
      if (fitted_model%k == 0) then
        call compute_mfgarch_inference(returns, period, fitted_model, raw_solution, result, &
          component_status, gradient_step=ctrl%gradient_step)
      else if (fitted_model%has_second) then
        call compute_mfgarch_inference(returns, period, fitted_model, raw_solution, result, &
          component_status, covariate, period_two, covariate_two, ctrl%gradient_step)
      else
        call compute_mfgarch_inference(returns, period, fitted_model, raw_solution, result, &
          component_status, covariate=covariate, gradient_step=ctrl%gradient_step)
      end if
      if (component_status /= mfgarch_success .and. result%status == mfgarch_success) then
        result%message = 'fit converged; covariance calculation failed'
      end if
    else
      allocate(result%covariance(0,0), result%robust_covariance(0,0), result%opg_covariance(0,0))
      allocate(result%standard_error(0), result%robust_standard_error(0), result%opg_standard_error(0))
    end if

  contains

    function objective_function(raw) result(value)
      real(dp), intent(in) :: raw(:)
      real(dp) :: value, llh
      type(mfgarch_model) :: trial_model
      integer :: local_status

      call raw_to_model(raw, start_model, trial_model, local_status)
      if (local_status /= mfgarch_success) then
        value = huge(1.0_dp) / 100.0_dp
        return
      end if
      if (trial_model%k == 0) then
        llh = log_likelihood(trial_model, returns, period, local_status)
      else if (trial_model%has_second) then
        llh = log_likelihood(trial_model, returns, period, local_status, covariate, &
          period_two, covariate_two)
      else
        llh = log_likelihood(trial_model, returns, period, local_status, covariate)
      end if
      if (local_status /= mfgarch_success .or. .not. finite_value(llh)) then
        value = huge(1.0_dp) / 100.0_dp
      else
        value = -llh
      end if
    end function objective_function

  end subroutine fit_mfgarch

  subroutine compute_mfgarch_inference(returns, period, template, raw, result, status, &
      covariate, period_two, covariate_two, gradient_step)
    real(dp), intent(in) :: returns(:), raw(:)
    integer, intent(in) :: period(:)
    type(mfgarch_model), intent(in) :: template
    type(mfgarch_fit_result), intent(inout) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: covariate(:), covariate_two(:), gradient_step
    integer, intent(in), optional :: period_two(:)
    real(dp), allocatable :: hessian(:,:), score(:,:), variance_score(:,:), meat(:,:), opg_meat(:,:), inv_hessian(:,:)
    real(dp), allocatable :: raw_cov(:,:), raw_robust(:,:), raw_opg(:,:), inv_meat(:,:)
    real(dp), allocatable :: jacobian(:,:), physical_plus(:), physical_minus(:)
    real(dp), allocatable :: c0(:), cp(:), cm(:), vp(:), vm(:)
    real(dp), allocatable :: xplus(:), xminus(:), xpp(:), xpm(:), xmp(:), xmm(:)
    real(dp) :: f0, fp, fm, fpp, fpm, fmp, fmm, hi, hj, step, factor, nan_value
    type(mfgarch_model) :: model_plus, model_minus
    integer :: npar, nobs, i, j, local_status, inv_status

    status = mfgarch_success
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    step = 1.0e-4_dp
    if (present(gradient_step)) step = max(gradient_step, 1.0e-6_dp)
    npar = size(raw)
    nobs = size(returns)
    allocate(hessian(npar,npar), score(nobs,npar), variance_score(nobs,npar), &
      meat(npar,npar), opg_meat(npar,npar), jacobian(npar,npar))
    allocate(xplus(npar), xminus(npar), xpp(npar), xpm(npar), xmp(npar), xmm(npar))
    hessian = 0.0_dp
    score = 0.0_dp
    variance_score = 0.0_dp

    call get_contributions(raw, c0, f0, local_status)
    if (local_status /= mfgarch_success) then
      status = local_status
      call allocate_failed_inference(result, npar, nan_value)
      return
    end if

    do j = 1, npar
      hi = step * max(1.0_dp, abs(raw(j)))
      xplus = raw
      xminus = raw
      xplus(j) = xplus(j) + hi
      xminus(j) = xminus(j) - hi
      call get_contributions(xplus, cp, fp, local_status)
      if (local_status /= mfgarch_success) then
        status = local_status
        call allocate_failed_inference(result, npar, nan_value)
        return
      end if
      call get_contributions(xminus, cm, fm, local_status)
      if (local_status /= mfgarch_success) then
        status = local_status
        call allocate_failed_inference(result, npar, nan_value)
        return
      end if
      call get_variance_contributions(xplus, vp, local_status)
      if (local_status /= mfgarch_success) then
        status = local_status
        call allocate_failed_inference(result, npar, nan_value)
        return
      end if
      call get_variance_contributions(xminus, vm, local_status)
      if (local_status /= mfgarch_success) then
        status = local_status
        call allocate_failed_inference(result, npar, nan_value)
        return
      end if
      do i = 1, nobs
        if (finite_value(cp(i)) .and. finite_value(cm(i))) score(i,j) = (cp(i) - cm(i)) / (2.0_dp * hi)
        if (finite_value(vp(i)) .and. finite_value(vm(i))) &
          variance_score(i,j) = (vp(i) - vm(i)) / (2.0_dp * hi)
      end do
      hessian(j,j) = (fp - 2.0_dp*f0 + fm) / (hi*hi)

      call raw_to_model(xplus, template, model_plus, local_status)
      call model_parameters(model_plus, physical_plus)
      call raw_to_model(xminus, template, model_minus, local_status)
      call model_parameters(model_minus, physical_minus)
      jacobian(:,j) = (physical_plus - physical_minus) / (2.0_dp * hi)
    end do

    do i = 1, npar - 1
      hi = step * max(1.0_dp, abs(raw(i)))
      do j = i + 1, npar
        hj = step * max(1.0_dp, abs(raw(j)))
        xpp = raw
        xpm = raw
        xmp = raw
        xmm = raw
        xpp(i) = xpp(i) + hi
        xpp(j) = xpp(j) + hj
        xpm(i) = xpm(i) + hi
        xpm(j) = xpm(j) - hj
        xmp(i) = xmp(i) - hi
        xmp(j) = xmp(j) + hj
        xmm(i) = xmm(i) - hi
        xmm(j) = xmm(j) - hj
        fpp = objective_only(xpp)
        fpm = objective_only(xpm)
        fmp = objective_only(xmp)
        fmm = objective_only(xmm)
        hessian(i,j) = (fpp - fpm - fmp + fmm) / (4.0_dp * hi * hj)
        hessian(j,i) = hessian(i,j)
      end do
    end do

    call invert_matrix(hessian, inv_hessian, inv_status, 1.0e-10_dp)
    if (inv_status /= mfgarch_success) then
      status = mfgarch_singular_matrix
      call allocate_failed_inference(result, npar, nan_value)
      return
    end if
    meat = matmul(transpose(score), score)
    opg_meat = matmul(transpose(variance_score), variance_score)
    raw_cov = inv_hessian
    raw_robust = matmul(matmul(inv_hessian, meat), inv_hessian)
    call invert_matrix(opg_meat, inv_meat, inv_status, 1.0e-12_dp)
    if (inv_status == mfgarch_success) then
      factor = residual_kurtosis_factor(result%residuals)
      raw_opg = factor * inv_meat
    else
      allocate(raw_opg(npar,npar))
      raw_opg = nan_value
    end if

    result%covariance = matmul(matmul(jacobian, raw_cov), transpose(jacobian))
    result%robust_covariance = matmul(matmul(jacobian, raw_robust), transpose(jacobian))
    if (matrix_is_finite(raw_opg)) then
      result%opg_covariance = matmul(matmul(jacobian, raw_opg), transpose(jacobian))
    else
      allocate(result%opg_covariance(npar,npar))
      result%opg_covariance = nan_value
    end if
    allocate(result%standard_error(npar), result%robust_standard_error(npar), &
      result%opg_standard_error(npar))
    do i = 1, npar
      result%standard_error(i) = sqrt(max(0.0_dp, result%covariance(i,i)))
      result%robust_standard_error(i) = sqrt(max(0.0_dp, result%robust_covariance(i,i)))
      if (finite_value(result%opg_covariance(i,i)) .and. result%opg_covariance(i,i) >= 0.0_dp) then
        result%opg_standard_error(i) = sqrt(result%opg_covariance(i,i))
      else
        result%opg_standard_error(i) = nan_value
      end if
    end do

  contains

    subroutine get_contributions(x, contributions, objective_value, local_status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: contributions(:)
      real(dp), intent(out) :: objective_value
      integer, intent(out) :: local_status
      type(mfgarch_model) :: model
      real(dp), allocatable :: local_tau(:), local_g(:), local_residuals(:)
      integer :: ii

      call raw_to_model(x, template, model, local_status)
      if (local_status /= mfgarch_success) return
      if (model%k == 0) then
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status)
      else if (model%has_second) then
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status, covariate, period_two, covariate_two)
      else
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status, covariate)
      end if
      if (local_status /= mfgarch_success) return
      objective_value = 0.0_dp
      do ii = 1, size(contributions)
        if (finite_value(contributions(ii))) objective_value = objective_value - contributions(ii)
      end do
    end subroutine get_contributions

    subroutine get_variance_contributions(x, contributions, local_status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: contributions(:)
      integer, intent(out) :: local_status
      type(mfgarch_model) :: model
      real(dp), allocatable :: local_tau(:), local_g(:), local_residuals(:)

      call raw_to_model(x, template, model, local_status)
      if (local_status /= mfgarch_success) return
      if (model%k == 0) then
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status, log_variance_only=.true.)
      else if (model%has_second) then
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status, covariate, period_two, covariate_two, .true.)
      else
        call likelihood_contributions(model, returns, period, contributions, local_tau, &
          local_g, local_residuals, local_status, covariate, log_variance_only=.true.)
      end if
    end subroutine get_variance_contributions

    function objective_only(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      real(dp), allocatable :: local_contributions(:)
      integer :: local_status

      call get_contributions(x, local_contributions, value, local_status)
      if (local_status /= mfgarch_success) value = huge(1.0_dp) / 100.0_dp
    end function objective_only

  end subroutine compute_mfgarch_inference

  subroutine allocate_failed_inference(result, npar, nan_value)
    type(mfgarch_fit_result), intent(inout) :: result
    integer, intent(in) :: npar
    real(dp), intent(in) :: nan_value

    allocate(result%covariance(npar,npar), result%robust_covariance(npar,npar), &
      result%opg_covariance(npar,npar), result%standard_error(npar), &
      result%robust_standard_error(npar), result%opg_standard_error(npar))
    result%covariance = nan_value
    result%robust_covariance = nan_value
    result%opg_covariance = nan_value
    result%standard_error = nan_value
    result%robust_standard_error = nan_value
    result%opg_standard_error = nan_value
  end subroutine allocate_failed_inference

  pure logical function matrix_is_finite(matrix) result(ok)
    real(dp), intent(in) :: matrix(:,:)
    integer :: i, j

    ok = .true.
    do j = 1, size(matrix,2)
      do i = 1, size(matrix,1)
        if (.not. finite_value(matrix(i,j))) then
          ok = .false.
          return
        end if
      end do
    end do
  end function matrix_is_finite

  function residual_kurtosis_factor(residuals) result(value)
    real(dp), intent(in) :: residuals(:)
    real(dp) :: value, fourth
    integer :: i, n

    fourth = 0.0_dp
    n = 0
    do i = 1, size(residuals)
      if (finite_value(residuals(i))) then
        fourth = fourth + residuals(i)**4
        n = n + 1
      end if
    end do
    if (n == 0) then
      value = 1.0_dp
    else
      value = max(0.0_dp, (fourth / real(n,dp) - 1.0_dp) / 2.0_dp)
    end if
  end function residual_kurtosis_factor

end module mfgarch_fit
