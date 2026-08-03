! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_solvers
  use maxlik_kinds, only: dp
  use maxlik_types, only: maxlik_problem, maxlik_control, maxlik_result
  use maxlik_status
  use maxlik_linalg, only: vector_norm, outer_product, identity_matrix, solve_linear
  use maxlik_evaluation, only: valid_problem, project_parameters, evaluate_value, evaluate_gradient, &
    evaluate_hessian, evaluate_scores
  use maxlik_random, only: maxlik_rng
  implicit none
  private

  public :: solve_newton, solve_bfgs, solve_bfgsr, solve_bhhh
  public :: solve_cg, solve_nelder_mead, solve_sann, solve_sga, solve_adam

contains

  subroutine prepare_result(problem, start, method, control, result, status)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    character(len=*), intent(in) :: method
    type(maxlik_control), intent(in) :: control
    type(maxlik_result), intent(out) :: result
    integer, intent(out) :: status

    status = MAXLIK_INVALID_INPUT
    if (.not. valid_problem(problem, start)) then
      result%code = MAXLIK_INVALID_START
      result%message = maxlik_message(result%code)
      return
    end if
    if (control%iterlim < 0 .or. control%gradtol < 0.0_dp .or. control%steptol < 0.0_dp) then
      result%code = MAXLIK_INVALID_INPUT
      result%message = maxlik_message(result%code)
      return
    end if
    allocate(result%estimate(problem%npar), result%gradient(problem%npar), &
      result%hessian(problem%npar, problem%npar), result%active(problem%npar))
    result%estimate = start
    result%gradient = 0.0_dp
    result%hessian = 0.0_dp
    result%active = problem%active
    result%method = method
    if (control%store_values) then
      allocate(result%stored_values(control%iterlim + 1))
      result%stored_values = 0.0_dp
    end if
    if (control%store_parameters) then
      allocate(result%stored_parameters(problem%npar, control%iterlim + 1))
      result%stored_parameters = 0.0_dp
    end if
    status = 0
  end subroutine prepare_result

  subroutine save_iteration(result, control, iteration, x, value)
    type(maxlik_result), intent(inout) :: result
    type(maxlik_control), intent(in) :: control
    integer, intent(in) :: iteration
    real(dp), intent(in) :: x(:), value
    integer :: position

    position = min(iteration + 1, control%iterlim + 1)
    if (allocated(result%stored_values)) result%stored_values(position) = value
    if (allocated(result%stored_parameters)) result%stored_parameters(:, position) = x
  end subroutine save_iteration

  subroutine set_convergence(result, code, iterations, x, value)
    type(maxlik_result), intent(inout) :: result
    integer, intent(in) :: code, iterations
    real(dp), intent(in) :: x(:), value

    result%code = code
    result%iterations = iterations
    result%estimate = x
    result%maximum = value
    result%message = maxlik_message(code)
    result%converged = code == MAXLIK_SUCCESS_GRADIENT .or. code == MAXLIK_SUCCESS_VALUE
  end subroutine set_convergence

  pure real(dp) function active_norm(problem, x) result(value)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    value = sqrt(sum(merge(x * x, 0.0_dp, problem%active)))
  end function active_norm


  pure subroutine projected_gradient(problem, x, gradient, projected)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:), gradient(:)
    real(dp), intent(out) :: projected(:)
    real(dp) :: tolerance
    integer :: i

    projected = gradient
    do i = 1, size(x)
      if (.not. problem%active(i)) then
        projected(i) = 0.0_dp
        cycle
      end if
      tolerance = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)))
      if (x(i) <= problem%lower(i) + tolerance .and. gradient(i) < 0.0_dp) projected(i) = 0.0_dp
      if (x(i) >= problem%upper(i) - tolerance .and. gradient(i) > 0.0_dp) projected(i) = 0.0_dp
    end do
  end subroutine projected_gradient

  pure subroutine make_direction_feasible(problem, x, direction)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: direction(:)
    real(dp) :: tolerance
    integer :: i

    do i = 1, size(x)
      if (.not. problem%active(i)) then
        direction(i) = 0.0_dp
        cycle
      end if
      tolerance = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)))
      if (x(i) <= problem%lower(i) + tolerance .and. direction(i) < 0.0_dp) direction(i) = 0.0_dp
      if (x(i) >= problem%upper(i) - tolerance .and. direction(i) > 0.0_dp) direction(i) = 0.0_dp
    end do
  end subroutine make_direction_feasible

  subroutine line_search(problem, control, x, value, gradient, direction, x_new, value_new, alpha, &
      function_count, status)
    type(maxlik_problem), intent(in) :: problem
    type(maxlik_control), intent(in) :: control
    real(dp), intent(in) :: x(:), value, gradient(:), direction(:)
    real(dp), intent(out) :: x_new(:), value_new, alpha
    integer, intent(inout) :: function_count
    integer, intent(out) :: status

    real(dp) :: slope, scale
    integer :: trial

    slope = dot_product(gradient, direction)
    if (slope <= 0.0_dp) then
      status = MAXLIK_STEP_FAILURE
      return
    end if
    alpha = 1.0_dp
    status = MAXLIK_STEP_FAILURE
    do trial = 1, 80
      x_new = x + alpha * direction
      where (.not. problem%active) x_new = x
      call project_parameters(problem, x_new)
      call evaluate_value(problem, x_new, value_new, function_count, status)
      if (status == 0) then
        if (value_new >= value + 1.0e-4_dp * alpha * slope) then
          status = 0
          return
        end if
      end if
      alpha = 0.5_dp * alpha
      scale = alpha * max(1.0_dp, active_norm(problem, direction))
      if (scale <= control%steptol) exit
    end do
    status = MAXLIK_STEP_FAILURE
  end subroutine line_search

  subroutine solve_bfgs(problem, start, control, result, method_name)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result
    character(len=*), intent(in), optional :: method_name

    type(maxlik_control) :: ctrl
    real(dp), allocatable :: x(:), x_new(:), gradient(:), gradient_new(:), projected(:), projected_new(:), direction(:)
    real(dp), allocatable :: s(:), y(:), inverse_info(:, :), ident(:, :)
    real(dp) :: value, value_new, alpha, ys, rho
    integer :: iteration, status
    character(len=40) :: label

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    label = 'BFGS maximization'
    if (present(method_name)) label = method_name
    call prepare_result(problem, start, label, ctrl, result, status)
    if (status /= 0) return

    allocate(x(problem%npar), x_new(problem%npar), gradient(problem%npar), gradient_new(problem%npar), &
      projected(problem%npar), projected_new(problem%npar), direction(problem%npar), s(problem%npar), y(problem%npar), &
      inverse_info(problem%npar, problem%npar), ident(problem%npar, problem%npar))
    x = start
    ident = identity_matrix(problem%npar)
    inverse_info = ident
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    call evaluate_gradient(problem, x, gradient, result%function_count, result%gradient_count, status, &
      ctrl%use_central_differences)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_EVALUATION_ERROR, 0, x, value)
      return
    end if
    call projected_gradient(problem, x, gradient, projected)
    call save_iteration(result, ctrl, 0, x, value)

    do iteration = 1, ctrl%iterlim
      if (active_norm(problem, projected) <= ctrl%gradtol) then
        call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, iteration - 1, x, value)
        return
      end if
      direction = matmul(inverse_info, projected)
      call make_direction_feasible(problem, x, direction)
      if (dot_product(projected, direction) <= epsilon(1.0_dp) * max(1.0_dp, active_norm(problem, projected))) then
        inverse_info = ident
        direction = projected
        call make_direction_feasible(problem, x, direction)
      end if
      call line_search(problem, ctrl, x, value, projected, direction, x_new, value_new, alpha, &
        result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_STEP_FAILURE, iteration - 1, x, value)
        return
      end if
      call evaluate_gradient(problem, x_new, gradient_new, result%function_count, result%gradient_count, status, &
        ctrl%use_central_differences)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration, x, value)
        return
      end if
      call projected_gradient(problem, x_new, gradient_new, projected_new)
      s = x_new - x
      y = projected - projected_new
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, active_norm(problem, y) * active_norm(problem, s))) then
        rho = 1.0_dp / ys
        inverse_info = matmul(ident - rho * outer_product(s, y), &
          matmul(inverse_info, ident - rho * outer_product(y, s))) + rho * outer_product(s, s)
      else
        inverse_info = ident
      end if
      call save_iteration(result, ctrl, iteration, x_new, value_new)
      if (abs(value_new - value) <= ctrl%reltol * (abs(value) + ctrl%reltol) .or. &
          active_norm(problem, s) <= ctrl%tol) then
        call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, x_new, value_new)
        return
      end if
      x = x_new
      value = value_new
      gradient = gradient_new
      projected = projected_new
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, x, value)
  end subroutine solve_bfgs

  subroutine solve_bfgsr(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    if (present(control)) then
      call solve_bfgs(problem, start, control, result, 'BFGSR maximization')
    else
      call solve_bfgs(problem, start, result=result, method_name='BFGSR maximization')
    end if
  end subroutine solve_bfgsr

  subroutine solve_cg(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    real(dp), allocatable :: x(:), x_new(:), gradient(:), gradient_new(:), projected(:), projected_new(:), direction(:), step(:)
    real(dp) :: value, value_new, alpha, beta, denominator
    integer :: iteration, status

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call prepare_result(problem, start, 'Conjugate-gradient maximization', ctrl, result, status)
    if (status /= 0) return
    allocate(x(problem%npar), x_new(problem%npar), gradient(problem%npar), gradient_new(problem%npar), &
      projected(problem%npar), projected_new(problem%npar), direction(problem%npar), step(problem%npar))
    x = start
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    call evaluate_gradient(problem, x, gradient, result%function_count, result%gradient_count, status, &
      ctrl%use_central_differences)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_EVALUATION_ERROR, 0, x, value)
      return
    end if
    call projected_gradient(problem, x, gradient, projected)
    direction = projected
    call make_direction_feasible(problem, x, direction)
    call save_iteration(result, ctrl, 0, x, value)

    do iteration = 1, ctrl%iterlim
      if (active_norm(problem, projected) <= ctrl%gradtol) then
        call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, iteration - 1, x, value)
        return
      end if
      if (dot_product(projected, direction) <= 0.0_dp) direction = projected
      call make_direction_feasible(problem, x, direction)
      call line_search(problem, ctrl, x, value, projected, direction, x_new, value_new, alpha, &
        result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_STEP_FAILURE, iteration - 1, x, value)
        return
      end if
      call evaluate_gradient(problem, x_new, gradient_new, result%function_count, result%gradient_count, status, &
        ctrl%use_central_differences)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration, x, value)
        return
      end if
      call projected_gradient(problem, x_new, gradient_new, projected_new)
      denominator = max(dot_product(projected, projected), tiny(1.0_dp))
      beta = max(0.0_dp, dot_product(projected_new, projected_new - projected) / denominator)
      direction = projected_new + beta * direction
      call make_direction_feasible(problem, x_new, direction)
      step = x_new - x
      call save_iteration(result, ctrl, iteration, x_new, value_new)
      if (abs(value_new - value) <= ctrl%reltol * (abs(value) + ctrl%reltol) .or. &
          active_norm(problem, step) <= ctrl%tol) then
        call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, x_new, value_new)
        return
      end if
      x = x_new
      value = value_new
      gradient = gradient_new
      projected = projected_new
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, x, value)
  end subroutine solve_cg

  subroutine solve_newton(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    real(dp), allocatable :: x(:), x_new(:), gradient(:), projected(:), hessian(:, :), information(:, :), direction(:), step(:)
    real(dp) :: value, value_new, alpha, lambda
    integer :: i, iteration, status, solve_status, trial

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call prepare_result(problem, start, 'Newton-Raphson maximization', ctrl, result, status)
    if (status /= 0) return
    allocate(x(problem%npar), x_new(problem%npar), gradient(problem%npar), projected(problem%npar), &
      hessian(problem%npar, problem%npar), information(problem%npar, problem%npar), &
      direction(problem%npar), step(problem%npar))
    x = start
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    call save_iteration(result, ctrl, 0, x, value)

    do iteration = 1, ctrl%iterlim
      call evaluate_gradient(problem, x, gradient, result%function_count, result%gradient_count, status, &
        ctrl%use_central_differences)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration - 1, x, value)
        return
      end if
      call projected_gradient(problem, x, gradient, projected)
      if (active_norm(problem, projected) <= ctrl%gradtol) then
        call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, iteration - 1, x, value)
        return
      end if
      call evaluate_hessian(problem, x, hessian, result%function_count, result%gradient_count, &
        result%hessian_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration - 1, x, value)
        return
      end if

      lambda = 0.0_dp
      direction = 0.0_dp
      do trial = 1, 60
        information = -hessian
        do i = 1, problem%npar
          if (problem%active(i)) then
            information(i, i) = information(i, i) + lambda
          else
            information(i, :) = 0.0_dp
            information(:, i) = 0.0_dp
            information(i, i) = 1.0_dp
          end if
        end do
        call solve_linear(information, projected, direction, solve_status, ctrl%qrtol)
        call make_direction_feasible(problem, x, direction)
        if (solve_status == 0 .and. dot_product(projected, direction) > 0.0_dp) exit
        if (lambda <= 0.0_dp) then
          lambda = ctrl%marquardt_lambda0
        else
          lambda = lambda * ctrl%marquardt_lambda_step
        end if
        if (lambda > ctrl%marquardt_max_lambda) exit
      end do
      if (solve_status /= 0 .or. dot_product(projected, direction) <= 0.0_dp) then
        call set_convergence(result, MAXLIK_SINGULAR_HESSIAN, iteration - 1, x, value)
        return
      end if

      call line_search(problem, ctrl, x, value, projected, direction, x_new, value_new, alpha, &
        result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_STEP_FAILURE, iteration - 1, x, value)
        return
      end if
      step = x_new - x
      call save_iteration(result, ctrl, iteration, x_new, value_new)
      if (abs(value_new - value) <= ctrl%reltol * (abs(value) + ctrl%reltol) .or. &
          active_norm(problem, step) <= ctrl%tol) then
        call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, x_new, value_new)
        return
      end if
      x = x_new
      value = value_new
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, x, value)
  end subroutine solve_newton

  subroutine solve_bhhh(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    real(dp), allocatable :: x(:), x_new(:), gradient(:), projected(:), scores(:, :), information(:, :), direction(:), step(:)
    real(dp) :: value, value_new, alpha
    integer :: i, iteration, status, solve_status

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call prepare_result(problem, start, 'BHHH maximization', ctrl, result, status)
    if (status /= 0) return
    if (.not. associated(problem%scores) .or. problem%nobs <= 0) then
      call set_convergence(result, MAXLIK_INVALID_INPUT, 0, start, -huge(1.0_dp))
      result%message = 'BHHH requires observation-level score contributions'
      return
    end if
    allocate(x(problem%npar), x_new(problem%npar), gradient(problem%npar), projected(problem%npar), &
      scores(problem%nobs, problem%npar), information(problem%npar, problem%npar), &
      direction(problem%npar), step(problem%npar))
    x = start
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    call save_iteration(result, ctrl, 0, x, value)

    do iteration = 1, ctrl%iterlim
      call evaluate_scores(problem, x, scores, result%gradient_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration - 1, x, value)
        return
      end if
      gradient = sum(scores, dim=1)
      if (problem%penalty_rho > 0.0_dp) then
        call evaluate_gradient(problem, x, gradient, result%function_count, result%gradient_count, status)
        if (status /= 0) then
          call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration - 1, x, value)
          return
        end if
      end if
      call projected_gradient(problem, x, gradient, projected)
      if (active_norm(problem, projected) <= ctrl%gradtol) then
        call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, iteration - 1, x, value)
        return
      end if
      information = matmul(transpose(scores), scores)
      do i = 1, problem%npar
        if (.not. problem%active(i)) then
          information(i, :) = 0.0_dp
          information(:, i) = 0.0_dp
          information(i, i) = 1.0_dp
        else
          information(i, i) = information(i, i) + ctrl%qrtol
        end if
      end do
      call solve_linear(information, projected, direction, solve_status, ctrl%qrtol)
      call make_direction_feasible(problem, x, direction)
      if (solve_status /= 0 .or. dot_product(projected, direction) <= 0.0_dp) then
        call set_convergence(result, MAXLIK_SINGULAR_HESSIAN, iteration - 1, x, value)
        return
      end if
      call line_search(problem, ctrl, x, value, projected, direction, x_new, value_new, alpha, &
        result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_STEP_FAILURE, iteration - 1, x, value)
        return
      end if
      step = x_new - x
      call save_iteration(result, ctrl, iteration, x_new, value_new)
      if (abs(value_new - value) <= ctrl%reltol * (abs(value) + ctrl%reltol) .or. &
          active_norm(problem, step) <= ctrl%tol) then
        call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, x_new, value_new)
        return
      end if
      x = x_new
      value = value_new
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, x, value)
  end subroutine solve_bhhh

  subroutine solve_nelder_mead(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
    real(dp), allocatable :: temp_x(:)
    real(dp) :: reflected_value, expanded_value, contracted_value, spread_value, spread_x, delta
    integer, allocatable :: active_index(:)
    integer :: i, j, m, iteration, status

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call prepare_result(problem, start, 'Nelder-Mead maximization', ctrl, result, status)
    if (status /= 0) return
    m = count(problem%active)
    if (m == 0) then
      call evaluate_value(problem, start, reflected_value, result%function_count, status)
      call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, 0, start, reflected_value)
      return
    end if
    allocate(active_index(m), simplex(problem%npar, m + 1), values(m + 1), centroid(problem%npar), &
      reflected(problem%npar), expanded(problem%npar), contracted(problem%npar), temp_x(problem%npar))
    active_index = pack([(i, i=1,problem%npar)], problem%active)
    simplex(:, 1) = start
    call evaluate_value(problem, simplex(:, 1), values(1), result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, start, values(1))
      return
    end if
    do j = 1, m
      simplex(:, j + 1) = start
      i = active_index(j)
      delta = 0.05_dp * max(1.0_dp, abs(start(i)))
      simplex(i, j + 1) = min(problem%upper(i), start(i) + delta)
      if (abs(simplex(i, j + 1) - start(i)) <= tiny(1.0_dp)) &
        simplex(i, j + 1) = max(problem%lower(i), start(i) - delta)
      call evaluate_value(problem, simplex(:, j + 1), values(j + 1), result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, 0, start, values(1))
        return
      end if
    end do
    call sort_simplex(simplex, values)
    call save_iteration(result, ctrl, 0, simplex(:, 1), values(1))

    do iteration = 1, ctrl%iterlim
      call sort_simplex(simplex, values)
      spread_value = maxval(abs(values - values(1)))
      spread_x = 0.0_dp
      do j = 2, m + 1
        spread_x = max(spread_x, active_norm(problem, simplex(:, j) - simplex(:, 1)))
      end do
      if (spread_value <= ctrl%reltol * (abs(values(1)) + ctrl%reltol) .or. spread_x <= ctrl%tol) then
        call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration - 1, simplex(:, 1), values(1))
        return
      end if
      centroid = sum(simplex(:, 1:m), dim=2) / real(m, dp)
      reflected = centroid + ctrl%nm_alpha * (centroid - simplex(:, m + 1))
      where (.not. problem%active) reflected = start
      call project_parameters(problem, reflected)
      call evaluate_value(problem, reflected, reflected_value, result%function_count, status)
      if (status /= 0) reflected_value = -huge(1.0_dp)

      if (reflected_value > values(1)) then
        expanded = centroid + ctrl%nm_gamma * (reflected - centroid)
        where (.not. problem%active) expanded = start
        call project_parameters(problem, expanded)
        call evaluate_value(problem, expanded, expanded_value, result%function_count, status)
        if (status == 0 .and. expanded_value > reflected_value) then
          simplex(:, m + 1) = expanded
          values(m + 1) = expanded_value
        else
          simplex(:, m + 1) = reflected
          values(m + 1) = reflected_value
        end if
      else if (reflected_value > values(m)) then
        simplex(:, m + 1) = reflected
        values(m + 1) = reflected_value
      else
        if (reflected_value > values(m + 1)) then
          contracted = centroid + ctrl%nm_beta * (reflected - centroid)
        else
          contracted = centroid + ctrl%nm_beta * (simplex(:, m + 1) - centroid)
        end if
        where (.not. problem%active) contracted = start
        call project_parameters(problem, contracted)
        call evaluate_value(problem, contracted, contracted_value, result%function_count, status)
        if (status == 0 .and. contracted_value > values(m + 1)) then
          simplex(:, m + 1) = contracted
          values(m + 1) = contracted_value
        else
          do j = 2, m + 1
            simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
            call project_parameters(problem, simplex(:, j))
            call evaluate_value(problem, simplex(:, j), values(j), result%function_count, status)
            if (status /= 0) values(j) = -huge(1.0_dp)
          end do
        end if
      end if
      call sort_simplex(simplex, values)
      call save_iteration(result, ctrl, iteration, simplex(:, 1), values(1))
    end do
    call sort_simplex(simplex, values)
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, simplex(:, 1), values(1))
  end subroutine solve_nelder_mead

  subroutine sort_simplex(simplex, values)
    real(dp), intent(inout) :: simplex(:, :), values(:)
    real(dp) :: value_key
    real(dp) :: x_key(size(simplex, 1))
    integer :: i, j

    do i = 2, size(values)
      value_key = values(i)
      x_key = simplex(:, i)
      j = i - 1
      do while (j >= 1)
        if (values(j) >= value_key) exit
        values(j + 1) = values(j)
        simplex(:, j + 1) = simplex(:, j)
        j = j - 1
      end do
      values(j + 1) = value_key
      simplex(:, j + 1) = x_key
    end do
  end subroutine sort_simplex

  subroutine solve_sann(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    type(maxlik_rng) :: rng
    real(dp), allocatable :: x(:), candidate(:), best(:)
    real(dp) :: value, candidate_value, best_value, temperature, delta, scale, log_uniform
    integer :: i, iteration, status

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call prepare_result(problem, start, 'Simulated annealing maximization', ctrl, result, status)
    if (status /= 0) return
    allocate(x(problem%npar), candidate(problem%npar), best(problem%npar))
    call rng%seed(ctrl%random_seed)
    x = start
    best = start
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    best_value = value
    call save_iteration(result, ctrl, 0, best, best_value)

    do iteration = 1, ctrl%iterlim
      temperature = ctrl%sann_temp / log(real(iteration / max(1, ctrl%sann_tmax) + 2, dp))
      scale = max(temperature, ctrl%steptol)
      candidate = x
      do i = 1, problem%npar
        if (problem%active(i)) candidate(i) = x(i) + scale * rng%normal()
      end do
      call project_parameters(problem, candidate)
      call evaluate_value(problem, candidate, candidate_value, result%function_count, status)
      if (status == 0) then
        delta = candidate_value - value
        if (delta >= 0.0_dp) then
          x = candidate
          value = candidate_value
        else
          log_uniform = log(rng%uniform())
          if (log_uniform < delta / max(temperature, tiny(1.0_dp))) then
            x = candidate
            value = candidate_value
          end if
        end if
        if (candidate_value > best_value) then
          best = candidate
          best_value = candidate_value
        end if
      end if
      call save_iteration(result, ctrl, iteration, best, best_value)
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, best, best_value)
    result%message = 'simulated annealing completed the requested iterations'
  end subroutine solve_sann

  subroutine solve_sga(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call solve_stochastic(problem, start, ctrl, result, .false.)
  end subroutine solve_sga

  subroutine solve_adam(problem, start, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in), optional :: control
    type(maxlik_result), intent(out) :: result

    type(maxlik_control) :: ctrl
    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    call solve_stochastic(problem, start, ctrl, result, .true.)
  end subroutine solve_adam

  subroutine solve_stochastic(problem, start, ctrl, result, use_adam)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_control), intent(in) :: ctrl
    type(maxlik_result), intent(out) :: result
    logical, intent(in) :: use_adam

    type(maxlik_rng) :: rng
    real(dp), allocatable :: x(:), gradient(:), scores(:, :), velocity(:), second_moment(:), best(:)
    integer, allocatable :: indices(:)
    real(dp) :: value, old_value, best_value, gradient_norm, beta1_power, beta2_power
    integer :: iteration, status, batch, i, no_improvement
    character(len=40) :: label

    label = 'Stochastic-gradient maximization'
    if (use_adam) label = 'Adam maximization'
    call prepare_result(problem, start, label, ctrl, result, status)
    if (status /= 0) return
    if (.not. associated(problem%scores) .or. problem%nobs <= 0) then
      call set_convergence(result, MAXLIK_INVALID_INPUT, 0, start, -huge(1.0_dp))
      result%message = 'stochastic methods require observation-level scores'
      return
    end if
    allocate(x(problem%npar), gradient(problem%npar), scores(problem%nobs, problem%npar), &
      velocity(problem%npar), second_moment(problem%npar), best(problem%npar), indices(problem%nobs))
    call rng%seed(ctrl%random_seed)
    indices = [(i, i=1,problem%nobs)]
    x = start
    velocity = 0.0_dp
    second_moment = 0.0_dp
    beta1_power = 1.0_dp
    beta2_power = 1.0_dp
    batch = ctrl%batch_size
    if (batch <= 0 .or. batch > problem%nobs) batch = problem%nobs
    call evaluate_value(problem, x, value, result%function_count, status)
    if (status /= 0) then
      call set_convergence(result, MAXLIK_INVALID_START, 0, x, value)
      return
    end if
    best = x
    best_value = value
    old_value = value
    no_improvement = 0
    call save_iteration(result, ctrl, 0, x, value)

    do iteration = 1, ctrl%iterlim
      call rng%shuffle(indices)
      call evaluate_scores(problem, x, scores, result%gradient_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration - 1, x, value)
        return
      end if
      gradient = real(problem%nobs, dp) / real(batch, dp) * sum(scores(indices(1:batch), :), dim=1)
      where (.not. problem%active) gradient = 0.0_dp
      call make_direction_feasible(problem, x, gradient)
      gradient_norm = active_norm(problem, gradient)
      if (ctrl%gradient_clip > 0.0_dp .and. gradient_norm > ctrl%gradient_clip) then
        gradient = gradient * ctrl%gradient_clip / gradient_norm
        gradient_norm = ctrl%gradient_clip
      end if
      if (gradient_norm <= ctrl%gradtol) then
        call set_convergence(result, MAXLIK_SUCCESS_GRADIENT, iteration - 1, x, value)
        return
      end if
      if (use_adam) then
        velocity = ctrl%adam_momentum1 * velocity + (1.0_dp - ctrl%adam_momentum1) * gradient
        second_moment = ctrl%adam_momentum2 * second_moment + &
          (1.0_dp - ctrl%adam_momentum2) * gradient * gradient
        beta1_power = beta1_power * ctrl%adam_momentum1
        beta2_power = beta2_power * ctrl%adam_momentum2
        x = x + ctrl%learning_rate * (velocity / max(1.0_dp - beta1_power, tiny(1.0_dp))) / &
          (sqrt(second_moment / max(1.0_dp - beta2_power, tiny(1.0_dp))) + 1.0e-8_dp)
      else
        velocity = ctrl%sga_momentum * velocity + ctrl%learning_rate * gradient
        x = x + velocity
      end if
      where (.not. problem%active) x = start
      call project_parameters(problem, x)
      call evaluate_value(problem, x, value, result%function_count, status)
      if (status /= 0) then
        call set_convergence(result, MAXLIK_EVALUATION_ERROR, iteration, best, best_value)
        return
      end if
      if (value > best_value) then
        best = x
        best_value = value
        no_improvement = 0
      else
        no_improvement = no_improvement + 1
      end if
      call save_iteration(result, ctrl, iteration, x, value)
      if (abs(value - old_value) <= ctrl%reltol * (abs(old_value) + ctrl%reltol)) then
        if (ctrl%patience <= 0) then
          call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, best, best_value)
          return
        end if
      end if
      if (ctrl%patience > 0 .and. mod(iteration, max(1, ctrl%patience_step)) == 0) then
        if (no_improvement >= ctrl%patience) then
          call set_convergence(result, MAXLIK_SUCCESS_VALUE, iteration, best, best_value)
          result%message = 'stopped after stochastic-gradient patience limit'
          return
        end if
      end if
      old_value = value
    end do
    call set_convergence(result, MAXLIK_ITERATION_LIMIT, ctrl%iterlim, best, best_value)
  end subroutine solve_stochastic

end module maxlik_solvers
