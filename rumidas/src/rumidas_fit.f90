! SPDX-License-Identifier: GPL-3.0-only
module rumidas_fit
  use rumidas_kinds, only: dp
  use rumidas_status
  use rumidas_types
  use rumidas_garch_midas, only: garch_midas_parameter_count, garch_midas_evaluate
  use rumidas_mem_models, only: mem_parameter_count, mem_evaluate
  use rumidas_statistics, only: sample_mean, sample_variance, information_criteria
  use maxlik, only: maxlik_problem, maxlik_control, maxlik_result, initialize_problem, set_bounds, &
    set_inequality_constraints, max_lik, robust_covariance_matrix, standard_errors
  implicit none
  private

  type :: garch_fit_context
    type(garch_midas_spec) :: spec
    real(dp), allocatable :: returns(:)
    real(dp), allocatable :: mv1(:, :)
    real(dp), allocatable :: mv2(:, :)
    real(dp), allocatable :: x(:)
  end type garch_fit_context

  type :: mem_fit_context
    type(mem_spec) :: spec
    real(dp), allocatable :: observations(:)
    real(dp), allocatable :: returns(:)
    real(dp), allocatable :: mv(:, :)
    real(dp), allocatable :: z(:)
  end type mem_fit_context

  type(garch_fit_context), save :: current_garch
  type(mem_fit_context), save :: current_mem

  public :: fit_garch_midas, fit_mem, ugmfit, umemfit

contains

  subroutine fit_garch_midas(spec, daily_ret, mv_m_1, result, status, mv_m_2, x_variable, start, control)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: daily_ret(:), mv_m_1(:, :)
    type(rumidas_fit_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: mv_m_2(:, :), x_variable(:), start(:)
    type(rumidas_fit_control), intent(in), optional :: control

    type(rumidas_fit_control) :: ctrl
    type(maxlik_problem) :: problem
    type(maxlik_control) :: max_control
    type(maxlik_result) :: max_result
    real(dp), allocatable :: lower(:), upper(:), a(:, :), b(:), initial(:), candidate(:)
    real(dp), allocatable :: ll(:), cv(:), lr(:), sr(:)
    real(dp) :: candidate_value, best_value
    integer :: npar, local_status, j, nstarts

    ctrl = rumidas_fit_control()
    if (present(control)) ctrl = control
    status = RUMIDAS_SUCCESS
    result = rumidas_fit_result()
    npar = garch_midas_parameter_count(spec)
    if (npar <= 0 .or. size(daily_ret) <= npar .or. size(mv_m_1, 2) /= size(daily_ret)) then
      status = RUMIDAS_INVALID_INPUT
      result%status = status
      result%message = 'invalid model specification or data dimensions'
      return
    end if
    if (present(start)) then
      if (size(start) /= npar) then
        status = RUMIDAS_DIMENSION_ERROR
        result%status = status
        result%message = 'starting vector has the wrong length'
        return
      end if
    end if

    call set_garch_context(spec, daily_ret, mv_m_1, mv_m_2, x_variable, local_status)
    if (local_status /= RUMIDAS_SUCCESS) then
      status = local_status
      result%status = status
      result%message = 'invalid optional data for model specification'
      return
    end if

    call initialize_problem(problem, npar, garch_objective, nobs=size(daily_ret))
    problem%scores => garch_scores
    allocate(lower(npar), upper(npar))
    call garch_bounds(spec, lower, upper)
    call set_bounds(problem, lower, upper, local_status)
    if (local_status /= 0) then
      status = RUMIDAS_INVALID_INPUT
      result%status = status
      result%message = 'failed to construct optimization bounds'
      return
    end if
    call garch_stationarity_constraint(spec, a, b)
    call set_inequality_constraints(problem, a, b, local_status)
    if (local_status /= 0) then
      status = RUMIDAS_INVALID_INPUT
      result%status = status
      result%message = 'failed to construct stationarity constraint'
      return
    end if

    allocate(initial(npar), candidate(npar))
    if (present(start)) then
      initial = start
    else
      call default_garch_start(spec, daily_ret, initial)
      best_value = -huge(1.0_dp)
      nstarts = max(1, ctrl%random_starts)
      do j = 1, nstarts
        call random_garch_start(spec, daily_ret, ctrl%random_seed + 7919 * j, candidate)
        call garch_objective(candidate, candidate_value, local_status)
        if (local_status == 0 .and. candidate_value > best_value) then
          best_value = candidate_value
          initial = candidate
        end if
      end do
    end if

    max_control = maxlik_control()
    max_control%iterlim = ctrl%max_iterations
    max_control%gradtol = ctrl%gradient_tolerance
    max_control%reltol = ctrl%relative_tolerance
    max_control%random_seed = ctrl%random_seed
    max_control%constraint_max_outer = 8
    max_control%constraint_tol = 1.0e-7_dp
    max_control%final_hessian = .true.
    call max_lik(problem, initial, max_result, trim(ctrl%method), max_control)
    if (.not. allocated(max_result%estimate)) then
      status = RUMIDAS_OPTIMIZATION_ERROR
      result%status = status
      result%message = trim(max_result%message)
      return
    end if

    result%coefficients = max_result%estimate
    result%loglik = max_result%maximum
    result%iterations = max_result%iterations
    result%function_count = max_result%function_count
    result%converged = max_result%converged
    result%message = trim(max_result%message)
    if (allocated(max_result%covariance)) then
      result%covariance = max_result%covariance
      if (allocated(max_result%std_error)) result%standard_errors = max_result%std_error
    end if
    if (ctrl%compute_robust_covariance .and. allocated(max_result%hessian) .and. &
        allocated(max_result%gradient_obs)) then
      call robust_covariance_matrix(max_result%hessian, max_result%gradient_obs, max_result%active, &
        result%robust_covariance, local_status)
      if (local_status == 0) then
        allocate(result%robust_standard_errors(npar))
        call standard_errors(result%robust_covariance, result%robust_standard_errors)
      end if
    end if

    call garch_midas_evaluate(result%coefficients, spec, daily_ret, mv_m_1, ll, cv, lr, sr, &
      local_status, mv_m_2, x_variable)
    if (local_status == RUMIDAS_SUCCESS) then
      result%loglik_obs = ll
      result%conditional = cv ** 2
      result%long_run = lr ** 2
      result%short_run = sr
      call information_criteria(result%loglik, npar, size(daily_ret), result%aic, result%bic)
    end if
    if (max_result%converged .and. local_status == RUMIDAS_SUCCESS) then
      status = RUMIDAS_SUCCESS
    else
      status = RUMIDAS_OPTIMIZATION_ERROR
    end if
    result%status = status
  end subroutine fit_garch_midas

  subroutine fit_mem(spec, x, result, status, daily_ret, mv_m, z_variable, start, control)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(in) :: x(:)
    type(rumidas_fit_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: daily_ret(:), mv_m(:, :), z_variable(:), start(:)
    type(rumidas_fit_control), intent(in), optional :: control

    type(rumidas_fit_control) :: ctrl
    type(maxlik_problem) :: problem
    type(maxlik_control) :: max_control
    type(maxlik_result) :: max_result
    real(dp), allocatable :: lower(:), upper(:), a(:, :), b(:), initial(:), candidate(:)
    real(dp), allocatable :: ll(:), pred(:), lr(:), sr(:)
    real(dp) :: candidate_value, best_value, variance_epsilon
    integer :: npar, local_status, j, nstarts

    ctrl = rumidas_fit_control()
    if (present(control)) ctrl = control
    status = RUMIDAS_SUCCESS
    result = rumidas_fit_result()
    npar = mem_parameter_count(spec)
    if (npar <= 0 .or. size(x) <= npar .or. any(x <= 0.0_dp)) then
      status = RUMIDAS_INVALID_INPUT
      result%status = status
      result%message = 'invalid model specification or nonpositive observations'
      return
    end if
    if (present(start)) then
      if (size(start) /= npar) then
        status = RUMIDAS_DIMENSION_ERROR
        result%status = status
        result%message = 'starting vector has the wrong length'
        return
      end if
    end if
    call set_mem_context(spec, x, daily_ret, mv_m, z_variable, local_status)
    if (local_status /= RUMIDAS_SUCCESS) then
      status = local_status
      result%status = status
      result%message = 'invalid optional data for MEM specification'
      return
    end if

    call initialize_problem(problem, npar, mem_objective, nobs=size(x))
    problem%scores => mem_scores
    allocate(lower(npar), upper(npar))
    call mem_bounds(spec, lower, upper)
    call set_bounds(problem, lower, upper, local_status)
    call mem_stationarity_constraint(spec, a, b)
    call set_inequality_constraints(problem, a, b, local_status)

    allocate(initial(npar), candidate(npar))
    if (present(start)) then
      initial = start
    else
      call default_mem_start(spec, x, initial)
      best_value = -huge(1.0_dp)
      nstarts = max(1, ctrl%random_starts)
      do j = 1, nstarts
        call random_mem_start(spec, x, ctrl%random_seed + 6151 * j, candidate)
        call mem_objective(candidate, candidate_value, local_status)
        if (local_status == 0 .and. candidate_value > best_value) then
          best_value = candidate_value
          initial = candidate
        end if
      end do
    end if

    max_control = maxlik_control()
    max_control%iterlim = ctrl%max_iterations
    max_control%gradtol = ctrl%gradient_tolerance
    max_control%reltol = ctrl%relative_tolerance
    max_control%constraint_max_outer = 8
    max_control%constraint_tol = 1.0e-7_dp
    max_control%final_hessian = .true.
    call max_lik(problem, initial, max_result, trim(ctrl%method), max_control)
    if (.not. allocated(max_result%estimate)) then
      status = RUMIDAS_OPTIMIZATION_ERROR
      result%status = status
      result%message = trim(max_result%message)
      return
    end if

    result%coefficients = max_result%estimate
    result%loglik = max_result%maximum
    result%iterations = max_result%iterations
    result%function_count = max_result%function_count
    result%converged = max_result%converged
    result%message = trim(max_result%message)
    if (allocated(max_result%covariance)) then
      result%covariance = max_result%covariance
      if (allocated(max_result%std_error)) result%standard_errors = max_result%std_error
    end if
    if (ctrl%compute_robust_covariance .and. allocated(max_result%hessian) .and. &
        allocated(max_result%gradient_obs)) then
      call robust_covariance_matrix(max_result%hessian, max_result%gradient_obs, max_result%active, &
        result%robust_covariance, local_status)
      if (local_status == 0) then
        allocate(result%robust_standard_errors(npar))
        call standard_errors(result%robust_covariance, result%robust_standard_errors)
      end if
    end if

    call mem_evaluate(result%coefficients, spec, x, ll, pred, lr, sr, local_status, &
      daily_ret, mv_m, z_variable)
    if (local_status == RUMIDAS_SUCCESS) then
      result%loglik_obs = ll
      result%conditional = pred
      result%long_run = lr
      result%short_run = sr
      call information_criteria(result%loglik, npar, size(x), result%aic, result%bic)
      if (.not. allocated(result%robust_standard_errors)) then
        variance_epsilon = sample_variance(x / pred)
        if (allocated(result%covariance)) then
          result%robust_covariance = 0.5_dp * variance_epsilon * result%covariance
          allocate(result%robust_standard_errors(npar))
          call standard_errors(result%robust_covariance, result%robust_standard_errors)
        end if
      end if
    end if
    if (max_result%converged .and. local_status == RUMIDAS_SUCCESS) then
      status = RUMIDAS_SUCCESS
    else
      status = RUMIDAS_OPTIMIZATION_ERROR
    end if
    result%status = status
  end subroutine fit_mem

  subroutine set_garch_context(spec, daily_ret, mv1, mv2, x, status)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: daily_ret(:), mv1(:, :)
    real(dp), intent(in), optional :: mv2(:, :), x(:)
    integer, intent(out) :: status
    status = RUMIDAS_SUCCESS
    current_garch%spec = spec
    current_garch%returns = daily_ret
    current_garch%mv1 = mv1
    if (allocated(current_garch%mv2)) deallocate(current_garch%mv2)
    if (allocated(current_garch%x)) deallocate(current_garch%x)
    if (present(mv2)) current_garch%mv2 = mv2
    if (present(x)) current_garch%x = x
    if ((spec%model == RUMIDAS_GM2M .or. spec%model == RUMIDAS_DAGM2M) .and. .not. allocated(current_garch%mv2)) &
      status = RUMIDAS_INVALID_INPUT
    if ((spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX) .and. .not. allocated(current_garch%x)) &
      status = RUMIDAS_INVALID_INPUT
  end subroutine set_garch_context

  subroutine set_mem_context(spec, x, daily_ret, mv, z, status)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: daily_ret(:), mv(:, :), z(:)
    integer, intent(out) :: status
    status = RUMIDAS_SUCCESS
    current_mem%spec = spec
    current_mem%observations = x
    if (allocated(current_mem%returns)) deallocate(current_mem%returns)
    if (allocated(current_mem%mv)) deallocate(current_mem%mv)
    if (allocated(current_mem%z)) deallocate(current_mem%z)
    if (present(daily_ret)) current_mem%returns = daily_ret
    if (present(mv)) current_mem%mv = mv
    if (present(z)) current_mem%z = z
    if (spec%skew .and. .not. allocated(current_mem%returns)) status = RUMIDAS_INVALID_INPUT
    if ((spec%model == RUMIDAS_MEM_MIDAS .or. spec%model == RUMIDAS_MEM_MIDAS_X) .and. &
        .not. allocated(current_mem%mv)) status = RUMIDAS_INVALID_INPUT
    if ((spec%model == RUMIDAS_MEM_X .or. spec%model == RUMIDAS_MEM_MIDAS_X) .and. &
        .not. allocated(current_mem%z)) status = RUMIDAS_INVALID_INPUT
  end subroutine set_mem_context

  subroutine garch_objective(parameters, value, callback_status)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: callback_status
    real(dp), allocatable :: ll(:), cv(:), lr(:), sr(:)
    integer :: status
    call evaluate_current_garch(parameters, ll, cv, lr, sr, status)
    callback_status = 0
    if (status == RUMIDAS_SUCCESS .and. all(ll > -huge(1.0_dp) / 10.0_dp)) then
      value = sum(ll)
    else
      value = -1.0e100_dp
    end if
  end subroutine garch_objective

  subroutine garch_scores(parameters, scores, callback_status)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: scores(:, :)
    integer, intent(out) :: callback_status
    real(dp), allocatable :: plus(:), minus(:), base(:), cv(:), lr(:), sr(:), xp(:), xm(:)
    real(dp) :: step
    integer :: j, status_plus, status_minus, status_base
    allocate(xp(size(parameters)), xm(size(parameters)))
    call evaluate_current_garch(parameters, base, cv, lr, sr, status_base)
    if (status_base /= RUMIDAS_SUCCESS) then
      scores = 0.0_dp
      callback_status = 1
      return
    end if
    do j = 1, size(parameters)
      step = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp) * (abs(parameters(j)) + 1.0_dp)
      xp = parameters; xm = parameters
      xp(j) = xp(j) + step; xm(j) = xm(j) - step
      call evaluate_current_garch(xp, plus, cv, lr, sr, status_plus)
      call evaluate_current_garch(xm, minus, cv, lr, sr, status_minus)
      if (status_plus == RUMIDAS_SUCCESS .and. status_minus == RUMIDAS_SUCCESS) then
        scores(:, j) = (plus - minus) / (2.0_dp * step)
      else if (status_plus == RUMIDAS_SUCCESS) then
        scores(:, j) = (plus - base) / step
      else if (status_minus == RUMIDAS_SUCCESS) then
        scores(:, j) = (base - minus) / step
      else
        scores(:, j) = 0.0_dp
      end if
    end do
    callback_status = 0
  end subroutine garch_scores

  subroutine evaluate_current_garch(parameters, ll, cv, lr, sr, status)
    real(dp), intent(in) :: parameters(:)
    real(dp), allocatable, intent(out) :: ll(:), cv(:), lr(:), sr(:)
    integer, intent(out) :: status
    if (allocated(current_garch%mv2) .and. allocated(current_garch%x)) then
      call garch_midas_evaluate(parameters, current_garch%spec, current_garch%returns, current_garch%mv1, &
        ll, cv, lr, sr, status, current_garch%mv2, current_garch%x)
    else if (allocated(current_garch%mv2)) then
      call garch_midas_evaluate(parameters, current_garch%spec, current_garch%returns, current_garch%mv1, &
        ll, cv, lr, sr, status, mv_m_2=current_garch%mv2)
    else if (allocated(current_garch%x)) then
      call garch_midas_evaluate(parameters, current_garch%spec, current_garch%returns, current_garch%mv1, &
        ll, cv, lr, sr, status, x_variable=current_garch%x)
    else
      call garch_midas_evaluate(parameters, current_garch%spec, current_garch%returns, current_garch%mv1, &
        ll, cv, lr, sr, status)
    end if
  end subroutine evaluate_current_garch

  subroutine mem_objective(parameters, value, callback_status)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: callback_status
    real(dp), allocatable :: ll(:), pred(:), lr(:), sr(:)
    integer :: status
    call evaluate_current_mem(parameters, ll, pred, lr, sr, status)
    callback_status = 0
    if (status == RUMIDAS_SUCCESS .and. all(ll > -huge(1.0_dp) / 10.0_dp)) then
      value = sum(ll)
    else
      value = -1.0e100_dp
    end if
  end subroutine mem_objective

  subroutine mem_scores(parameters, scores, callback_status)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(out) :: scores(:, :)
    integer, intent(out) :: callback_status
    real(dp), allocatable :: plus(:), minus(:), base(:), pred(:), lr(:), sr(:), xp(:), xm(:)
    real(dp) :: step
    integer :: j, status_plus, status_minus, status_base
    allocate(xp(size(parameters)), xm(size(parameters)))
    call evaluate_current_mem(parameters, base, pred, lr, sr, status_base)
    if (status_base /= RUMIDAS_SUCCESS) then
      scores = 0.0_dp
      callback_status = 1
      return
    end if
    do j = 1, size(parameters)
      step = epsilon(1.0_dp) ** (1.0_dp / 3.0_dp) * (abs(parameters(j)) + 1.0_dp)
      xp = parameters; xm = parameters
      xp(j) = xp(j) + step; xm(j) = xm(j) - step
      call evaluate_current_mem(xp, plus, pred, lr, sr, status_plus)
      call evaluate_current_mem(xm, minus, pred, lr, sr, status_minus)
      if (status_plus == RUMIDAS_SUCCESS .and. status_minus == RUMIDAS_SUCCESS) then
        scores(:, j) = (plus - minus) / (2.0_dp * step)
      else if (status_plus == RUMIDAS_SUCCESS) then
        scores(:, j) = (plus - base) / step
      else if (status_minus == RUMIDAS_SUCCESS) then
        scores(:, j) = (base - minus) / step
      else
        scores(:, j) = 0.0_dp
      end if
    end do
    callback_status = 0
  end subroutine mem_scores

  subroutine evaluate_current_mem(parameters, ll, pred, lr, sr, status)
    real(dp), intent(in) :: parameters(:)
    real(dp), allocatable, intent(out) :: ll(:), pred(:), lr(:), sr(:)
    integer, intent(out) :: status
    if (allocated(current_mem%returns) .and. allocated(current_mem%mv) .and. allocated(current_mem%z)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        current_mem%returns, current_mem%mv, current_mem%z)
    else if (allocated(current_mem%returns) .and. allocated(current_mem%mv)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        daily_ret=current_mem%returns, mv_m=current_mem%mv)
    else if (allocated(current_mem%returns) .and. allocated(current_mem%z)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        daily_ret=current_mem%returns, z_variable=current_mem%z)
    else if (allocated(current_mem%mv) .and. allocated(current_mem%z)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        mv_m=current_mem%mv, z_variable=current_mem%z)
    else if (allocated(current_mem%returns)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        daily_ret=current_mem%returns)
    else if (allocated(current_mem%mv)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        mv_m=current_mem%mv)
    else if (allocated(current_mem%z)) then
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status, &
        z_variable=current_mem%z)
    else
      call mem_evaluate(parameters, current_mem%spec, current_mem%observations, ll, pred, lr, sr, status)
    end if
  end subroutine evaluate_current_mem

  subroutine garch_bounds(spec, lower, upper)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(out) :: lower(:), upper(:)
    integer :: index, j, number_long
    number_long = 0
    lower = -huge(1.0_dp)
    upper = huge(1.0_dp)
    lower(1) = 1.0e-5_dp; upper(1) = 0.999_dp
    lower(2) = 1.0e-5_dp; upper(2) = 0.999_dp
    index = 3
    if (spec%skew) index = index + 1
    if (spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX) index = index + 1
    index = index + 1
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      number_long = 1
    case (RUMIDAS_GM2M, RUMIDAS_DAGM, RUMIDAS_DAGMX)
      number_long = 2
    case (RUMIDAS_DAGM2M)
      number_long = 4
    end select
    do j = 1, number_long
      index = index + 1
      if (spec%lag_function == RUMIDAS_BETA_LAG) lower(index) = 1.001_dp
      if (spec%lag_function == RUMIDAS_ALMON_LAG) upper(index) = -1.0e-6_dp
      index = index + 1
    end do
    if (spec%distribution == RUMIDAS_STUDENT_T) lower(size(lower)) = 2.001_dp
  end subroutine garch_bounds

  subroutine garch_stationarity_constraint(spec, a, b)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: a(:, :), b(:)
    allocate(a(1, garch_midas_parameter_count(spec)), b(1))
    a = 0.0_dp
    a(1, 1) = -1.0_dp
    a(1, 2) = -1.0_dp
    if (spec%skew) a(1, 3) = -0.5_dp
    b(1) = 0.999_dp
  end subroutine garch_stationarity_constraint

  subroutine mem_bounds(spec, lower, upper)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(out) :: lower(:), upper(:)
    integer :: index
    lower = -huge(1.0_dp)
    upper = huge(1.0_dp)
    lower(1) = 1.0e-5_dp; upper(1) = 0.999_dp
    lower(2) = 1.0e-5_dp; upper(2) = 0.999_dp
    index = 3
    if (spec%skew) index = index + 1
    if (spec%model == RUMIDAS_MEM_MIDAS .or. spec%model == RUMIDAS_MEM_MIDAS_X) then
      index = index + 2
      lower(index) = 1.001_dp
    end if
  end subroutine mem_bounds

  subroutine mem_stationarity_constraint(spec, a, b)
    type(mem_spec), intent(in) :: spec
    real(dp), allocatable, intent(out) :: a(:, :), b(:)
    allocate(a(1, mem_parameter_count(spec)), b(1))
    a = 0.0_dp
    a(1, 1) = -1.0_dp
    a(1, 2) = -1.0_dp
    if (spec%skew) a(1, 3) = -0.5_dp
    b(1) = 0.999_dp
  end subroutine mem_stationarity_constraint

  subroutine default_garch_start(spec, daily_ret, start)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: daily_ret(:)
    real(dp), intent(out) :: start(:)
    integer :: index, j, number_long
    number_long = 0
    start = 0.0_dp
    start(1) = 0.05_dp
    start(2) = 0.85_dp
    index = 3
    if (spec%skew) then
      start(index) = 0.05_dp
      index = index + 1
    end if
    if (spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX) then
      start(index) = 0.01_dp
      index = index + 1
    end if
    start(index) = log(max(sample_variance(daily_ret), 1.0e-6_dp))
    index = index + 1
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      number_long = 1
    case (RUMIDAS_GM2M, RUMIDAS_DAGM, RUMIDAS_DAGMX)
      number_long = 2
    case (RUMIDAS_DAGM2M)
      number_long = 4
    end select
    do j = 1, number_long
      start(index) = 0.05_dp
      start(index + 1) = merge(2.0_dp, -0.1_dp, spec%lag_function == RUMIDAS_BETA_LAG)
      index = index + 2
    end do
    if (spec%distribution == RUMIDAS_STUDENT_T) start(index) = 8.0_dp
  end subroutine default_garch_start

  subroutine random_garch_start(spec, daily_ret, seed, start)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: daily_ret(:)
    integer, intent(in) :: seed
    real(dp), intent(out) :: start(:)
    integer :: state, index, j, number_long
    real(dp) :: alpha, beta, gamma_value
    number_long = 0
    state = max(1, seed)
    alpha = 0.01_dp + 0.12_dp * lcg_uniform(state)
    beta = 0.60_dp + 0.25_dp * lcg_uniform(state)
    gamma_value = 0.10_dp * lcg_uniform(state)
    if (alpha + beta + 0.5_dp * gamma_value >= 0.97_dp) beta = 0.80_dp
    start = 0.0_dp
    start(1) = alpha; start(2) = beta
    index = 3
    if (spec%skew) then
      start(index) = gamma_value
      index = index + 1
    end if
    if (spec%model == RUMIDAS_GMX .or. spec%model == RUMIDAS_DAGMX) then
      start(index) = -0.05_dp + 0.10_dp * lcg_uniform(state)
      index = index + 1
    end if
    start(index) = log(max(sample_variance(daily_ret), 1.0e-6_dp)) - 0.5_dp + lcg_uniform(state)
    index = index + 1
    select case (spec%model)
    case (RUMIDAS_GM, RUMIDAS_GMX)
      number_long = 1
    case (RUMIDAS_GM2M, RUMIDAS_DAGM, RUMIDAS_DAGMX)
      number_long = 2
    case (RUMIDAS_DAGM2M)
      number_long = 4
    end select
    do j = 1, number_long
      start(index) = -0.5_dp + lcg_uniform(state)
      if (spec%lag_function == RUMIDAS_BETA_LAG) then
        start(index + 1) = 1.05_dp + 8.0_dp * lcg_uniform(state)
      else
        start(index + 1) = -0.01_dp - 0.5_dp * lcg_uniform(state)
      end if
      index = index + 2
    end do
    if (spec%distribution == RUMIDAS_STUDENT_T) start(index) = 2.5_dp + 12.0_dp * lcg_uniform(state)
  end subroutine random_garch_start

  subroutine default_mem_start(spec, x, start)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: start(:)
    integer :: index
    start = 0.0_dp
    start(1) = 0.25_dp; start(2) = 0.65_dp
    index = 3
    if (spec%skew) then
      start(index) = 0.05_dp
      index = index + 1
    end if
    if (spec%model == RUMIDAS_MEM_MIDAS .or. spec%model == RUMIDAS_MEM_MIDAS_X) then
      start(index) = log(max(sample_mean(x), 1.0e-6_dp))
      start(index + 1) = 0.05_dp
      start(index + 2) = 2.0_dp
      index = index + 3
    end if
    if (spec%model == RUMIDAS_MEM_X .or. spec%model == RUMIDAS_MEM_MIDAS_X) start(index) = 0.01_dp
  end subroutine default_mem_start

  subroutine random_mem_start(spec, x, seed, start)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: seed
    real(dp), intent(out) :: start(:)
    integer :: state, index
    real(dp) :: alpha, beta, gamma_value
    state = max(1, seed)
    alpha = 0.05_dp + 0.35_dp * lcg_uniform(state)
    beta = 0.40_dp + 0.45_dp * lcg_uniform(state)
    gamma_value = 0.10_dp * lcg_uniform(state)
    if (alpha + beta + 0.5_dp * gamma_value >= 0.97_dp) beta = 0.55_dp
    start = 0.0_dp
    start(1) = alpha; start(2) = beta
    index = 3
    if (spec%skew) then
      start(index) = gamma_value
      index = index + 1
    end if
    if (spec%model == RUMIDAS_MEM_MIDAS .or. spec%model == RUMIDAS_MEM_MIDAS_X) then
      start(index) = log(max(sample_mean(x), 1.0e-6_dp)) - 0.5_dp + lcg_uniform(state)
      start(index + 1) = -0.5_dp + lcg_uniform(state)
      start(index + 2) = 1.05_dp + 8.0_dp * lcg_uniform(state)
      index = index + 3
    end if
    if (spec%model == RUMIDAS_MEM_X .or. spec%model == RUMIDAS_MEM_MIDAS_X) &
      start(index) = -0.05_dp + 0.10_dp * lcg_uniform(state)
  end subroutine random_mem_start

  real(dp) function lcg_uniform(state) result(value)
    integer, intent(inout) :: state
    integer(kind=8) :: next
    next = modulo(1103515245_8 * int(state, 8) + 12345_8, 2147483647_8)
    state = int(max(next, 1_8))
    value = real(state, dp) / 2147483647.0_dp
  end function lcg_uniform

  subroutine ugmfit(spec, daily_ret, mv_m_1, result, status, mv_m_2, x_variable, start, control)
    type(garch_midas_spec), intent(in) :: spec
    real(dp), intent(in) :: daily_ret(:), mv_m_1(:, :)
    type(rumidas_fit_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: mv_m_2(:, :), x_variable(:), start(:)
    type(rumidas_fit_control), intent(in), optional :: control
    call fit_garch_midas(spec, daily_ret, mv_m_1, result, status, mv_m_2, x_variable, start, control)
  end subroutine ugmfit

  subroutine umemfit(spec, x, result, status, daily_ret, mv_m, z_variable, start, control)
    type(mem_spec), intent(in) :: spec
    real(dp), intent(in) :: x(:)
    type(rumidas_fit_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), intent(in), optional :: daily_ret(:), mv_m(:, :), z_variable(:), start(:)
    type(rumidas_fit_control), intent(in), optional :: control
    call fit_mem(spec, x, result, status, daily_ret, mv_m, z_variable, start, control)
  end subroutine umemfit

end module rumidas_fit
