! SPDX-License-Identifier: GPL-3.0-only
!
! Modern Fortran translation of the computational core of NlcOptim 0.6.
module nlcoptim
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp, qp_success, qp_inconsistent_constraints, &
    qp_not_positive_definite
  implicit none
  private

  integer, parameter, public :: nlc_success = 0
  integer, parameter, public :: nlc_max_iterations = 1
  integer, parameter, public :: nlc_max_evaluations = 2
  integer, parameter, public :: nlc_line_search_failed = 3
  integer, parameter, public :: nlc_invalid_input = 4
  integer, parameter, public :: nlc_qp_failed = 5
  integer, parameter, public :: nlc_nonfinite_evaluation = 6
  integer, parameter, public :: nlc_constraint_size_changed = 7

  type, public :: nlc_options
    real(dp) :: tol_x = 1.0e-5_dp
    real(dp) :: tol_fun = 1.0e-6_dp
    real(dp) :: tol_con = 1.0e-6_dp
    integer :: max_function_evaluations = 10000000
    integer :: max_iterations = 4000
    integer :: max_line_search = 30
    real(dp) :: armijo = 1.0e-4_dp
    real(dp) :: minimum_step = 1.0e-8_dp
    real(dp) :: finite_difference_step = sqrt(epsilon(1.0_dp))
    logical :: central_differences = .false.
    logical :: use_elastic_fallback = .true.
    real(dp) :: elastic_penalty = 1.0e4_dp
    real(dp) :: hessian_regularization = 1.0e-10_dp
  end type nlc_options

  type, public :: nlc_multipliers
    real(dp), allocatable :: equality_linear(:)
    real(dp), allocatable :: equality_nonlinear(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp), allocatable :: inequality_linear(:)
    real(dp), allocatable :: inequality_nonlinear(:)
  end type nlc_multipliers

  type, public :: nlc_result
    real(dp), allocatable :: x(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: hessian(:, :)
    type(nlc_multipliers) :: lambda
    integer :: function_evaluations = 0
    integer :: constraint_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: iterations = 0
    integer :: qp_iterations = 0
    real(dp) :: max_constraint_violation = huge(1.0_dp)
    real(dp) :: kkt_error = huge(1.0_dp)
    logical :: used_elastic_qp = .false.
    integer :: status = nlc_invalid_input
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => nlc_result_succeeded
  end type nlc_result

  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function

    subroutine constraint_function(x, equality, inequality)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: equality(:)
      real(dp), allocatable, intent(out) :: inequality(:)
    end subroutine constraint_function

    subroutine gradient_function(x, gradient)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
    end subroutine gradient_function

    subroutine jacobian_function(x, equality_jacobian, inequality_jacobian)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: equality_jacobian(:, :)
      real(dp), allocatable, intent(out) :: inequality_jacobian(:, :)
    end subroutine jacobian_function
  end interface

  public :: solnl
  public :: objective_function, constraint_function, gradient_function, jacobian_function

contains

  subroutine solnl(x0, objfun, fit, confun, a, b, aeq, beq, lb, ub, options, &
      gradfun, jacfun)
    real(dp), intent(in) :: x0(:)
    procedure(objective_function) :: objfun
    type(nlc_result), intent(out) :: fit
    procedure(constraint_function), optional :: confun
    real(dp), intent(in), optional :: a(:, :), b(:), aeq(:, :), beq(:), lb(:), ub(:)
    type(nlc_options), intent(in), optional :: options
    procedure(gradient_function), optional :: gradfun
    procedure(jacobian_function), optional :: jacfun

    type(nlc_options) :: opt
    real(dp), allocatable :: x(:), x_new(:), g(:), g_new(:), hess(:, :)
    real(dp), allocatable :: ceq(:), cineq(:), ceq_new(:), cineq_new(:)
    real(dp), allocatable :: jeq(:, :), jineq(:, :), jeq_new(:, :), jineq_new(:, :)
    real(dp), allocatable :: amat(:, :), avec(:), aeqmat(:, :), aeqvec(:)
    real(dp), allocatable :: lower(:), upper(:), step(:), lambda_all(:)
    real(dp), allocatable :: grad_l_old(:), grad_l_new(:), s(:), y(:)
    real(dp) :: f, f_new, alpha, violation, violation_new, merit, merit_new
    real(dp) :: rho, directional, step_norm, best_merit, best_f, reg
    real(dp), allocatable :: best_x(:), best_ceq(:), best_cineq(:)
    integer :: n, meq_nl, mineq_nl, iter, ls, qp_stat, qp_iters
    logical :: accepted, valid, elastic_used

    opt = nlc_options()
    if (present(options)) opt = options
    call initialize_result(fit, size(x0))

    n = size(x0)
    if (n < 1) then
      call set_failure(fit, nlc_invalid_input, 'x0 must contain at least one variable.')
      return
    end if
    if (.not. all(ieee_is_finite(x0))) then
      call set_failure(fit, nlc_invalid_input, 'x0 must be finite.')
      return
    end if
    if (opt%max_iterations < 1 .or. opt%max_function_evaluations < 1 .or. &
        opt%tol_x <= 0.0_dp .or. opt%tol_fun <= 0.0_dp .or. opt%tol_con <= 0.0_dp) then
      call set_failure(fit, nlc_invalid_input, 'Invalid optimization options.')
      return
    end if

    call prepare_linear_inputs(n, a, b, aeq, beq, lb, ub, amat, avec, aeqmat, &
      aeqvec, lower, upper, valid, fit%message)
    if (.not. valid) then
      fit%status = nlc_invalid_input
      return
    end if

    allocate(x(n), x_new(n), g(n), g_new(n), hess(n, n), step(n))
    allocate(grad_l_old(n), grad_l_new(n), s(n), y(n))
    x = min(max(x0, lower), upper)
    hess = 0.0_dp
    hess = identity_matrix(n)

    call evaluate_point(x, objfun, confun, f, ceq, cineq, fit, valid)
    if (.not. valid) then
      call set_failure(fit, nlc_nonfinite_evaluation, &
        'The objective or constraints were nonfinite at the starting point.')
      fit%x = x
      return
    end if
    meq_nl = size(ceq)
    mineq_nl = size(cineq)

    call evaluate_derivatives(x, f, ceq, cineq, objfun, confun, lower, upper, &
      opt, g, jeq, jineq, fit, valid, gradfun, jacfun)
    if (.not. valid) then
      call set_failure(fit, nlc_nonfinite_evaluation, &
        'Could not evaluate finite derivatives at the starting point.')
      fit%x = x
      return
    end if

    violation = constraint_violation(x, ceq, cineq, amat, avec, aeqmat, aeqvec, &
      lower, upper)
    allocate(lambda_all(size(aeqvec) + meq_nl + count(ieee_is_finite(lower)) + &
      count(ieee_is_finite(upper)) + size(avec) + mineq_nl))
    lambda_all = 0.0_dp
    elastic_used = .false.
    reg = max(opt%hessian_regularization, epsilon(1.0_dp))

    do iter = 1, opt%max_iterations
      fit%iterations = iter
      if (fit%function_evaluations >= opt%max_function_evaluations) then
        call set_failure(fit, nlc_max_evaluations, &
          'Maximum objective-function evaluations reached.')
        exit
      end if

      call solve_sqp_subproblem(x, g, hess, ceq, cineq, jeq, jineq, amat, avec, &
        aeqmat, aeqvec, lower, upper, opt, step, lambda_all, qp_stat, qp_iters, &
        elastic_used)
      fit%qp_iterations = fit%qp_iterations + qp_iters
      fit%used_elastic_qp = fit%used_elastic_qp .or. elastic_used
      if (qp_stat /= qp_success) then
        hess = identity_matrix(n)
        hess = hess * max(1.0_dp, reg * 100.0_dp)
        call solve_sqp_subproblem(x, g, hess, ceq, cineq, jeq, jineq, amat, avec, &
          aeqmat, aeqvec, lower, upper, opt, step, lambda_all, qp_stat, qp_iters, &
          elastic_used)
        fit%qp_iterations = fit%qp_iterations + qp_iters
        fit%used_elastic_qp = fit%used_elastic_qp .or. elastic_used
      end if
      if (qp_stat /= qp_success) then
        call set_failure(fit, nlc_qp_failed, 'The SQP quadratic subproblem failed.')
        exit
      end if

      call lagrangian_gradient(g, jeq, jineq, aeqmat, amat, lower, upper, &
        lambda_all, grad_l_old)
      fit%kkt_error = maxval(abs(grad_l_old))
      step_norm = maxval(abs(step))
      if (violation <= opt%tol_con .and. fit%kkt_error <= opt%tol_fun .and. &
          step_norm <= 2.0_dp * opt%tol_x) then
        fit%status = nlc_success
        fit%message = 'success'
        exit
      end if

      rho = max(10.0_dp, 1.0_dp + maxval(abs(lambda_all)))
      merit = f + rho * l1_constraint_violation(x, ceq, cineq, amat, avec, &
        aeqmat, aeqvec, lower, upper)
      directional = dot_product(g, step) - rho * &
        l1_constraint_violation(x, ceq, cineq, amat, avec, aeqmat, aeqvec, &
          lower, upper)
      if (directional >= -epsilon(1.0_dp)) directional = -max(1.0_dp, abs(merit))

      alpha = 1.0_dp
      accepted = .false.
      best_merit = huge(1.0_dp)
      best_f = huge(1.0_dp)
      allocate(best_x(n))
      best_x = x
      allocate(best_ceq(meq_nl), best_cineq(mineq_nl))
      best_ceq = ceq
      best_cineq = cineq

      do ls = 1, opt%max_line_search
        x_new = min(max(x + alpha * step, lower), upper)
        call evaluate_point(x_new, objfun, confun, f_new, ceq_new, cineq_new, fit, valid)
        if (valid) then
          if (size(ceq_new) /= meq_nl .or. size(cineq_new) /= mineq_nl) then
            call set_failure(fit, nlc_constraint_size_changed, &
              'The nonlinear constraint count changed during optimization.')
            deallocate(best_x, best_ceq, best_cineq)
            return
          end if
          merit_new = f_new + rho * l1_constraint_violation(x_new, ceq_new, &
            cineq_new, amat, avec, aeqmat, aeqvec, lower, upper)
          if (merit_new < best_merit) then
            best_merit = merit_new
            best_f = f_new
            best_x = x_new
            best_ceq = ceq_new
            best_cineq = cineq_new
          end if
          if (merit_new <= merit + opt%armijo * alpha * directional) then
            accepted = .true.
            exit
          end if
        end if
        alpha = 0.5_dp * alpha
        if (alpha < opt%minimum_step) exit
      end do

      if (.not. accepted) then
        if (best_merit < merit) then
          x_new = best_x
          f_new = best_f
          ceq_new = best_ceq
          cineq_new = best_cineq
          alpha = max(opt%minimum_step, alpha)
          accepted = .true.
        end if
      end if
      deallocate(best_x, best_ceq, best_cineq)

      if (.not. accepted) then
        if (violation <= opt%tol_con .and. step_norm <= opt%tol_x) then
          fit%status = nlc_success
          fit%message = 'success'
        else
          call set_failure(fit, nlc_line_search_failed, &
            'The merit-function line search failed.')
        end if
        exit
      end if

      call evaluate_derivatives(x_new, f_new, ceq_new, cineq_new, objfun, confun, &
        lower, upper, opt, g_new, jeq_new, jineq_new, fit, valid, gradfun, jacfun)
      if (.not. valid) then
        call set_failure(fit, nlc_nonfinite_evaluation, &
          'Could not evaluate finite derivatives after an SQP step.')
        exit
      end if

      call lagrangian_gradient(g_new, jeq_new, jineq_new, aeqmat, amat, lower, &
        upper, lambda_all, grad_l_new)
      s = x_new - x
      y = grad_l_new - grad_l_old
      call damped_bfgs_update(hess, s, y)

      x = x_new
      f = f_new
      g = g_new
      ceq = ceq_new
      cineq = cineq_new
      jeq = jeq_new
      jineq = jineq_new
      violation_new = constraint_violation(x, ceq, cineq, amat, avec, aeqmat, &
        aeqvec, lower, upper)

      if (maxval(abs(s)) <= opt%tol_x .and. abs(dot_product(g, s)) <= opt%tol_fun &
          .and. violation_new <= opt%tol_con) then
        violation = violation_new
        fit%status = nlc_success
        fit%message = 'success'
        exit
      end if
      violation = violation_new
    end do

    if (fit%status == nlc_invalid_input) then
      call set_failure(fit, nlc_max_iterations, 'Maximum SQP iterations reached.')
    end if

    fit%x = x
    fit%objective = f
    fit%gradient = g
    fit%hessian = hess
    fit%max_constraint_violation = violation
    call lagrangian_gradient(g, jeq, jineq, aeqmat, amat, lower, upper, lambda_all, &
      grad_l_old)
    fit%kkt_error = maxval(abs(grad_l_old))
    call unpack_multipliers(lambda_all, size(aeqvec), meq_nl, lower, upper, &
      size(avec), mineq_nl, fit%lambda)
  end subroutine solnl

  subroutine evaluate_point(x, objfun, confun, f, ceq, cineq, fit, valid)
    real(dp), intent(in) :: x(:)
    procedure(objective_function) :: objfun
    procedure(constraint_function), optional :: confun
    real(dp), intent(out) :: f
    real(dp), allocatable, intent(out) :: ceq(:), cineq(:)
    type(nlc_result), intent(inout) :: fit
    logical, intent(out) :: valid

    f = objfun(x)
    fit%function_evaluations = fit%function_evaluations + 1
    if (present(confun)) then
      call confun(x, ceq, cineq)
      fit%constraint_evaluations = fit%constraint_evaluations + 1
      if (.not. allocated(ceq)) allocate(ceq(0))
      if (.not. allocated(cineq)) allocate(cineq(0))
    else
      allocate(ceq(0), cineq(0))
    end if
    valid = ieee_is_finite(f) .and. all(ieee_is_finite(ceq)) .and. &
      all(ieee_is_finite(cineq))
  end subroutine evaluate_point

  subroutine evaluate_derivatives(x, f, ceq, cineq, objfun, confun, lower, upper, &
      opt, g, jeq, jineq, fit, valid, gradfun, jacfun)
    real(dp), intent(in) :: x(:), f, ceq(:), cineq(:), lower(:), upper(:)
    procedure(objective_function) :: objfun
    procedure(constraint_function), optional :: confun
    type(nlc_options), intent(in) :: opt
    real(dp), allocatable, intent(out) :: g(:), jeq(:, :), jineq(:, :)
    type(nlc_result), intent(inout) :: fit
    logical, intent(out) :: valid
    procedure(gradient_function), optional :: gradfun
    procedure(jacobian_function), optional :: jacfun

    integer :: n, i
    real(dp) :: h, fp, fm
    real(dp), allocatable :: xp(:), xm(:), ep(:), ip(:), em(:), im(:)
    logical :: can_forward, can_backward

    n = size(x)
    allocate(g(n), jeq(size(ceq), n), jineq(size(cineq), n))
    g = 0.0_dp
    jeq = 0.0_dp
    jineq = 0.0_dp

    if (present(gradfun)) then
      call gradfun(x, g)
    end if
    if (present(jacfun)) then
      call jacfun(x, jeq, jineq)
      if (.not. allocated(jeq)) allocate(jeq(0, n))
      if (.not. allocated(jineq)) allocate(jineq(0, n))
      if (size(jeq, 1) /= size(ceq) .or. size(jeq, 2) /= n .or. &
          size(jineq, 1) /= size(cineq) .or. size(jineq, 2) /= n) then
        valid = .false.
        return
      end if
    end if

    if (present(gradfun) .and. (present(jacfun) .or. .not. present(confun))) then
      fit%gradient_evaluations = fit%gradient_evaluations + 1
      valid = all(ieee_is_finite(g)) .and. all(ieee_is_finite(jeq)) .and. &
        all(ieee_is_finite(jineq))
      return
    end if

    allocate(xp(n), xm(n))
    do i = 1, n
      h = opt%finite_difference_step * max(1.0_dp, abs(x(i)))
      can_forward = x(i) + h <= upper(i)
      can_backward = x(i) - h >= lower(i)
      if (.not. can_forward .and. .not. can_backward) then
        h = min(upper(i) - x(i), x(i) - lower(i))
        if (h <= 0.0_dp) cycle
        can_forward = x(i) + h <= upper(i)
        can_backward = x(i) - h >= lower(i)
      end if

      if (opt%central_differences .and. can_forward .and. can_backward) then
        xp = x
        xm = x
        xp(i) = xp(i) + h
        xm(i) = xm(i) - h
        if (.not. present(gradfun)) then
          fp = objfun(xp)
          fm = objfun(xm)
          fit%function_evaluations = fit%function_evaluations + 2
          g(i) = (fp - fm) / (2.0_dp * h)
        end if
        if (present(confun) .and. .not. present(jacfun)) then
          call confun(xp, ep, ip)
          call confun(xm, em, im)
          if (.not. allocated(ep)) allocate(ep(0))
          if (.not. allocated(ip)) allocate(ip(0))
          if (.not. allocated(em)) allocate(em(0))
          if (.not. allocated(im)) allocate(im(0))
          fit%constraint_evaluations = fit%constraint_evaluations + 2
          if (size(ep) /= size(ceq) .or. size(em) /= size(ceq) .or. &
              size(ip) /= size(cineq) .or. size(im) /= size(cineq)) then
            valid = .false.
            return
          end if
          if (size(ceq) > 0) jeq(:, i) = (ep - em) / (2.0_dp * h)
          if (size(cineq) > 0) jineq(:, i) = (ip - im) / (2.0_dp * h)
        end if
      else
        xp = x
        if (can_forward) then
          xp(i) = xp(i) + h
          if (.not. present(gradfun)) then
            fp = objfun(xp)
            fit%function_evaluations = fit%function_evaluations + 1
            g(i) = (fp - f) / h
          end if
          if (present(confun) .and. .not. present(jacfun)) then
            call confun(xp, ep, ip)
            if (.not. allocated(ep)) allocate(ep(0))
            if (.not. allocated(ip)) allocate(ip(0))
            fit%constraint_evaluations = fit%constraint_evaluations + 1
            if (size(ep) /= size(ceq) .or. size(ip) /= size(cineq)) then
              valid = .false.
              return
            end if
            if (size(ceq) > 0) jeq(:, i) = (ep - ceq) / h
            if (size(cineq) > 0) jineq(:, i) = (ip - cineq) / h
          end if
        else
          xm = x
          xm(i) = xm(i) - h
          if (.not. present(gradfun)) then
            fm = objfun(xm)
            fit%function_evaluations = fit%function_evaluations + 1
            g(i) = (f - fm) / h
          end if
          if (present(confun) .and. .not. present(jacfun)) then
            call confun(xm, em, im)
            if (.not. allocated(em)) allocate(em(0))
            if (.not. allocated(im)) allocate(im(0))
            fit%constraint_evaluations = fit%constraint_evaluations + 1
            if (size(em) /= size(ceq) .or. size(im) /= size(cineq)) then
              valid = .false.
              return
            end if
            if (size(ceq) > 0) jeq(:, i) = (ceq - em) / h
            if (size(cineq) > 0) jineq(:, i) = (cineq - im) / h
          end if
        end if
      end if
    end do
    fit%gradient_evaluations = fit%gradient_evaluations + 1
    valid = all(ieee_is_finite(g)) .and. all(ieee_is_finite(jeq)) .and. &
      all(ieee_is_finite(jineq))
  end subroutine evaluate_derivatives

  subroutine solve_sqp_subproblem(x, g, hess, ceq, cineq, jeq, jineq, a, b, &
      aeq, beq, lower, upper, opt, step, lambda, status, iterations, elastic_used)
    real(dp), intent(in) :: x(:), g(:), hess(:, :), ceq(:), cineq(:)
    real(dp), intent(in) :: jeq(:, :), jineq(:, :), a(:, :), b(:), aeq(:, :)
    real(dp), intent(in) :: beq(:), lower(:), upper(:)
    type(nlc_options), intent(in) :: opt
    real(dp), intent(out) :: step(:), lambda(:)
    integer, intent(out) :: status, iterations
    logical, intent(out) :: elastic_used

    type(qp_result) :: qp
    real(dp), allocatable :: h(:, :), amat(:, :), bvec(:)
    integer :: n, neq, q, col, i, nf_lower, nf_upper

    n = size(x)
    neq = size(beq) + size(ceq)
    nf_lower = count(ieee_is_finite(lower))
    nf_upper = count(ieee_is_finite(upper))
    q = neq + nf_lower + nf_upper + size(b) + size(cineq)
    allocate(h(n, n), amat(n, q), bvec(q))
    h = 0.5_dp * (hess + transpose(hess))
    do i = 1, n
      h(i, i) = h(i, i) + opt%hessian_regularization
    end do
    amat = 0.0_dp
    bvec = 0.0_dp
    col = 0
    do i = 1, size(beq)
      col = col + 1
      amat(:, col) = aeq(i, :)
      bvec(col) = beq(i) - dot_product(aeq(i, :), x)
    end do
    do i = 1, size(ceq)
      col = col + 1
      amat(:, col) = jeq(i, :)
      bvec(col) = -ceq(i)
    end do
    do i = 1, n
      if (ieee_is_finite(lower(i))) then
        col = col + 1
        amat(i, col) = 1.0_dp
        bvec(col) = lower(i) - x(i)
      end if
    end do
    do i = 1, n
      if (ieee_is_finite(upper(i))) then
        col = col + 1
        amat(i, col) = -1.0_dp
        bvec(col) = x(i) - upper(i)
      end if
    end do
    do i = 1, size(b)
      col = col + 1
      amat(:, col) = -a(i, :)
      bvec(col) = dot_product(a(i, :), x) - b(i)
    end do
    do i = 1, size(cineq)
      col = col + 1
      amat(:, col) = -jineq(i, :)
      bvec(col) = cineq(i)
    end do

    qp = solve_qp(h, -g, amat, bvec, meq=neq)
    iterations = sum(qp%iterations)
    elastic_used = .false.
    if (qp%status == qp_success) then
      step = qp%solution
      lambda = qp%lagrangian
      status = qp_success
      return
    end if

    if (opt%use_elastic_fallback .and. &
        (qp%status == qp_inconsistent_constraints .or. &
         qp%status == qp_not_positive_definite)) then
      call solve_elastic_subproblem(x, g, h, ceq, cineq, jeq, jineq, a, b, aeq, &
        beq, lower, upper, opt, step, lambda, status, iterations)
      elastic_used = status == qp_success
    else
      step = 0.0_dp
      lambda = 0.0_dp
      status = qp%status
    end if
  end subroutine solve_sqp_subproblem

  subroutine solve_elastic_subproblem(x, g, h, ceq, cineq, jeq, jineq, a, b, &
      aeq, beq, lower, upper, opt, step, lambda, status, iterations)
    real(dp), intent(in) :: x(:), g(:), h(:, :), ceq(:), cineq(:)
    real(dp), intent(in) :: jeq(:, :), jineq(:, :), a(:, :), b(:), aeq(:, :)
    real(dp), intent(in) :: beq(:), lower(:), upper(:)
    type(nlc_options), intent(in) :: opt
    real(dp), intent(out) :: step(:), lambda(:)
    integer, intent(out) :: status, iterations

    type(qp_result) :: qp
    real(dp), allocatable :: he(:, :), de(:), ae(:, :), bv(:)
    integer :: n, q, col, i, nf_lower, nf_upper, neq

    n = size(x)
    neq = size(beq) + size(ceq)
    nf_lower = count(ieee_is_finite(lower))
    nf_upper = count(ieee_is_finite(upper))
    q = 2 * neq + nf_lower + nf_upper + size(b) + size(cineq) + 1
    allocate(he(n + 1, n + 1), de(n + 1), ae(n + 1, q), bv(q))
    he = 0.0_dp
    he(1:n, 1:n) = h
    he(n + 1, n + 1) = max(opt%hessian_regularization, 1.0e-10_dp)
    de(1:n) = -g
    de(n + 1) = -opt%elastic_penalty
    ae = 0.0_dp
    bv = 0.0_dp
    col = 0

    do i = 1, size(beq)
      col = col + 1
      ae(1:n, col) = -aeq(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = dot_product(aeq(i, :), x) - beq(i)
      col = col + 1
      ae(1:n, col) = aeq(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = beq(i) - dot_product(aeq(i, :), x)
    end do
    do i = 1, size(ceq)
      col = col + 1
      ae(1:n, col) = -jeq(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = ceq(i)
      col = col + 1
      ae(1:n, col) = jeq(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = -ceq(i)
    end do
    do i = 1, n
      if (ieee_is_finite(lower(i))) then
        col = col + 1
        ae(i, col) = 1.0_dp
        ae(n + 1, col) = 1.0_dp
        bv(col) = lower(i) - x(i)
      end if
    end do
    do i = 1, n
      if (ieee_is_finite(upper(i))) then
        col = col + 1
        ae(i, col) = -1.0_dp
        ae(n + 1, col) = 1.0_dp
        bv(col) = x(i) - upper(i)
      end if
    end do
    do i = 1, size(b)
      col = col + 1
      ae(1:n, col) = -a(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = dot_product(a(i, :), x) - b(i)
    end do
    do i = 1, size(cineq)
      col = col + 1
      ae(1:n, col) = -jineq(i, :)
      ae(n + 1, col) = 1.0_dp
      bv(col) = cineq(i)
    end do
    col = col + 1
    ae(n + 1, col) = 1.0_dp
    bv(col) = 0.0_dp

    qp = solve_qp(he, de, ae, bv, meq=0)
    iterations = iterations + sum(qp%iterations)
    status = qp%status
    if (status == qp_success) then
      step = qp%solution(1:n)
      lambda = 0.0_dp
    else
      step = 0.0_dp
      lambda = 0.0_dp
    end if
  end subroutine solve_elastic_subproblem

  subroutine damped_bfgs_update(h, s, y)
    real(dp), intent(inout) :: h(:, :)
    real(dp), intent(in) :: s(:), y(:)
    real(dp), allocatable :: hs(:), ybar(:)
    real(dp) :: sty, sths, theta
    integer :: i

    allocate(hs(size(s)), ybar(size(y)))
    hs = matmul(h, s)
    sths = dot_product(s, hs)
    sty = dot_product(s, y)
    if (sths <= epsilon(1.0_dp) * max(1.0_dp, dot_product(s, s))) then
      h = identity_matrix(size(s))
      return
    end if
    if (sty < 0.2_dp * sths) then
      theta = 0.8_dp * sths / max(sths - sty, epsilon(1.0_dp))
      ybar = theta * y + (1.0_dp - theta) * hs
    else
      ybar = y
    end if
    sty = dot_product(s, ybar)
    if (sty <= epsilon(1.0_dp) * max(1.0_dp, dot_product(s, s))) return
    h = h - outer_product(hs, hs) / sths + outer_product(ybar, ybar) / sty
    h = 0.5_dp * (h + transpose(h))
    do i = 1, size(s)
      h(i, i) = max(h(i, i), 1.0e-12_dp)
    end do
  end subroutine damped_bfgs_update

  subroutine lagrangian_gradient(g, jeq, jineq, aeq, a, lower, upper, lambda, gl)
    real(dp), intent(in) :: g(:), jeq(:, :), jineq(:, :), aeq(:, :), a(:, :)
    real(dp), intent(in) :: lower(:), upper(:), lambda(:)
    real(dp), intent(out) :: gl(:)
    integer :: k, i

    gl = g
    k = 0
    do i = 1, size(aeq, 1)
      k = k + 1
      gl = gl - lambda(k) * aeq(i, :)
    end do
    do i = 1, size(jeq, 1)
      k = k + 1
      gl = gl - lambda(k) * jeq(i, :)
    end do
    do i = 1, size(lower)
      if (ieee_is_finite(lower(i))) then
        k = k + 1
        gl(i) = gl(i) - lambda(k)
      end if
    end do
    do i = 1, size(upper)
      if (ieee_is_finite(upper(i))) then
        k = k + 1
        gl(i) = gl(i) + lambda(k)
      end if
    end do
    do i = 1, size(a, 1)
      k = k + 1
      gl = gl + lambda(k) * a(i, :)
    end do
    do i = 1, size(jineq, 1)
      k = k + 1
      gl = gl + lambda(k) * jineq(i, :)
    end do
  end subroutine lagrangian_gradient

  function constraint_violation(x, ceq, cineq, a, b, aeq, beq, lower, upper) result(v)
    real(dp), intent(in) :: x(:), ceq(:), cineq(:), a(:, :), b(:), aeq(:, :)
    real(dp), intent(in) :: beq(:), lower(:), upper(:)
    real(dp) :: v
    integer :: i

    v = 0.0_dp
    if (size(ceq) > 0) v = max(v, maxval(abs(ceq)))
    if (size(cineq) > 0) v = max(v, max(0.0_dp, maxval(cineq)))
    if (size(beq) > 0) v = max(v, maxval(abs(matmul(aeq, x) - beq)))
    if (size(b) > 0) v = max(v, max(0.0_dp, maxval(matmul(a, x) - b)))
    do i = 1, size(x)
      if (ieee_is_finite(lower(i))) v = max(v, lower(i) - x(i))
      if (ieee_is_finite(upper(i))) v = max(v, x(i) - upper(i))
    end do
  end function constraint_violation

  function l1_constraint_violation(x, ceq, cineq, a, b, aeq, beq, lower, upper) result(v)
    real(dp), intent(in) :: x(:), ceq(:), cineq(:), a(:, :), b(:), aeq(:, :)
    real(dp), intent(in) :: beq(:), lower(:), upper(:)
    real(dp) :: v
    integer :: i

    v = sum(abs(ceq)) + sum(max(0.0_dp, cineq))
    if (size(beq) > 0) v = v + sum(abs(matmul(aeq, x) - beq))
    if (size(b) > 0) v = v + sum(max(0.0_dp, matmul(a, x) - b))
    do i = 1, size(x)
      if (ieee_is_finite(lower(i))) v = v + max(0.0_dp, lower(i) - x(i))
      if (ieee_is_finite(upper(i))) v = v + max(0.0_dp, x(i) - upper(i))
    end do
  end function l1_constraint_violation

  subroutine prepare_linear_inputs(n, a_in, b_in, aeq_in, beq_in, lb_in, ub_in, &
      a, b, aeq, beq, lower, upper, valid, message)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: a_in(:, :), b_in(:), aeq_in(:, :), beq_in(:)
    real(dp), intent(in), optional :: lb_in(:), ub_in(:)
    real(dp), allocatable, intent(out) :: a(:, :), b(:), aeq(:, :), beq(:)
    real(dp), allocatable, intent(out) :: lower(:), upper(:)
    logical, intent(out) :: valid
    character(len=:), allocatable, intent(inout) :: message
    real(dp) :: inf

    inf = ieee_value(0.0_dp, ieee_positive_inf)
    valid = .false.
    allocate(lower(n), upper(n))
    lower = -inf
    upper = inf
    if (present(lb_in)) then
      if (size(lb_in) /= n) then
        message = 'lb must have the same length as x0.'
        return
      end if
      lower = lb_in
    end if
    if (present(ub_in)) then
      if (size(ub_in) /= n) then
        message = 'ub must have the same length as x0.'
        return
      end if
      upper = ub_in
    end if
    if (any(lower > upper) .or. any(ieee_is_nan(lower)) .or. &
        any(ieee_is_nan(upper))) then
      message = 'Invalid variable bounds.'
      return
    end if

    if (present(a_in)) then
      if (size(a_in, 2) /= n .or. .not. present(b_in)) then
        message = 'A must have n columns and requires B.'
        return
      end if
      if (size(b_in) /= size(a_in, 1)) then
        message = 'B must have one element per row of A.'
        return
      end if
      allocate(a(size(a_in, 1), n), b(size(b_in)))
      a = a_in
      b = b_in
    else
      if (present(b_in)) then
        message = 'B was supplied without A.'
        return
      end if
      allocate(a(0, n), b(0))
    end if

    if (present(aeq_in)) then
      if (size(aeq_in, 2) /= n .or. .not. present(beq_in)) then
        message = 'Aeq must have n columns and requires Beq.'
        return
      end if
      if (size(beq_in) /= size(aeq_in, 1)) then
        message = 'Beq must have one element per row of Aeq.'
        return
      end if
      allocate(aeq(size(aeq_in, 1), n), beq(size(beq_in)))
      aeq = aeq_in
      beq = beq_in
    else
      if (present(beq_in)) then
        message = 'Beq was supplied without Aeq.'
        return
      end if
      allocate(aeq(0, n), beq(0))
    end if

    if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b)) .or. &
        .not. all(ieee_is_finite(aeq)) .or. .not. all(ieee_is_finite(beq))) then
      message = 'Linear constraints must be finite.'
      return
    end if
    valid = .true.
  end subroutine prepare_linear_inputs

  subroutine unpack_multipliers(lambda, nleq, nne, lower, upper, nli, nni, out)
    real(dp), intent(in) :: lambda(:), lower(:), upper(:)
    integer, intent(in) :: nleq, nne, nli, nni
    type(nlc_multipliers), intent(out) :: out
    integer :: k, i

    allocate(out%equality_linear(nleq), out%equality_nonlinear(nne))
    allocate(out%lower(size(lower)), out%upper(size(upper)))
    allocate(out%inequality_linear(nli), out%inequality_nonlinear(nni))
    out%lower = 0.0_dp
    out%upper = 0.0_dp
    k = 0
    if (nleq > 0) then
      out%equality_linear = -lambda(1:nleq)
      k = nleq
    end if
    if (nne > 0) then
      out%equality_nonlinear = -lambda(k + 1:k + nne)
      k = k + nne
    end if
    do i = 1, size(lower)
      if (ieee_is_finite(lower(i))) then
        k = k + 1
        out%lower(i) = lambda(k)
      end if
    end do
    do i = 1, size(upper)
      if (ieee_is_finite(upper(i))) then
        k = k + 1
        out%upper(i) = lambda(k)
      end if
    end do
    if (nli > 0) then
      out%inequality_linear = lambda(k + 1:k + nli)
      k = k + nli
    end if
    if (nni > 0) out%inequality_nonlinear = lambda(k + 1:k + nni)
  end subroutine unpack_multipliers

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity_matrix

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: i
    do i = 1, size(x)
      a(i, :) = x(i) * y
    end do
  end function outer_product

  subroutine initialize_result(fit, n)
    type(nlc_result), intent(out) :: fit
    integer, intent(in) :: n
    allocate(fit%x(n), fit%gradient(n), fit%hessian(n, n))
    fit%x = ieee_value(0.0_dp, ieee_quiet_nan)
    fit%gradient = ieee_value(0.0_dp, ieee_quiet_nan)
    fit%hessian = ieee_value(0.0_dp, ieee_quiet_nan)
    fit%status = nlc_invalid_input
    fit%message = 'not started'
  end subroutine initialize_result

  subroutine set_failure(fit, status, message)
    type(nlc_result), intent(inout) :: fit
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    fit%status = status
    fit%message = message
  end subroutine set_failure

  pure logical function nlc_result_succeeded(this)
    class(nlc_result), intent(in) :: this
    nlc_result_succeeded = this%status == nlc_success
  end function nlc_result_succeeded

end module nlcoptim
