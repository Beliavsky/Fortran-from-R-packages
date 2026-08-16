! Modern Fortran computational port of roptim 0.1.7.
! Copyright (C) 2018 Yi Pan <ypan1988@gmail.com>
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-2.0-or-later
!
! The L-BFGS-B numerical kernel is separately licensed. See LICENSES.


module roptim_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
       ieee_positive_inf, ieee_negative_inf
  use, intrinsic :: iso_fortran_env, only : output_unit, int64
  use lbfgsb3_core_mod, only : setulb, set_core_output
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  character(len=*), parameter, public :: method_nelder_mead = 'Nelder-Mead'
  character(len=*), parameter, public :: method_bfgs = 'BFGS'
  character(len=*), parameter, public :: method_cg = 'CG'
  character(len=*), parameter, public :: method_lbfgsb = 'L-BFGS-B'
  character(len=*), parameter, public :: method_sann = 'SANN'

  integer, parameter, public :: roptim_success = 0
  integer, parameter, public :: roptim_iteration_limit = 1
  integer, parameter, public :: roptim_invalid_input = 10
  integer, parameter, public :: roptim_nonfinite = 11
  integer, parameter, public :: roptim_line_search_failure = 12
  integer, parameter, public :: roptim_user_stop = 13

  integer, parameter :: task_new_x = 1
  integer, parameter :: task_start = 2
  integer, parameter :: task_stop = 3
  integer, parameter :: task_fg = 4
  integer, parameter :: task_fg_line_search = 20
  integer, parameter :: task_fg_start = 21

  type, public :: roptim_control_t
    integer :: trace = 0
    real(dp) :: fnscale = 1.0_dp
    real(dp), allocatable :: parscale(:)
    real(dp), allocatable :: ndeps(:)
    integer :: max_iterations = 0
    real(dp) :: abstol = -huge(1.0_dp)
    real(dp) :: reltol = sqrt(epsilon(1.0_dp))
    real(dp) :: alpha = 1.0_dp
    real(dp) :: beta = 0.5_dp
    real(dp) :: gamma = 2.0_dp
    integer :: report = 10
    integer :: cg_type = 1
    integer :: lbfgsb_memory = 5
    real(dp) :: factr = 1.0e7_dp
    real(dp) :: pgtol = 0.0_dp
    real(dp) :: temperature = 10.0_dp
    integer :: tmax = 10
    integer :: seed = 0
    logical :: compute_hessian = .false.
    integer :: max_line_search = 40
  end type roptim_control_t

  type, public :: roptim_result_t
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: function_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: iterations = 0
    integer :: convergence = roptim_invalid_input
    logical :: success = .false.
    character(len=:), allocatable :: message
    real(dp), allocatable :: hessian(:, :)
  end type roptim_result_t

  abstract interface
    function objective_interface(x, user_data) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
    end function objective_interface

    subroutine gradient_interface(x, g, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: user_data
    end subroutine gradient_interface

    subroutine proposal_interface(current, candidate, scale, user_data)
      import dp
      real(dp), intent(in) :: current(:)
      real(dp), intent(out) :: candidate(:)
      real(dp), intent(in) :: scale
      class(*), intent(inout), optional :: user_data
    end subroutine proposal_interface

    subroutine monitor_interface(x, f, iteration, evaluations, stop, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
      class(*), intent(inout), optional :: user_data
    end subroutine monitor_interface

    function simple_objective_interface(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function simple_objective_interface

    subroutine simple_gradient_interface(x, g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine simple_gradient_interface

    subroutine simple_proposal_interface(current, candidate, scale)
      import dp
      real(dp), intent(in) :: current(:)
      real(dp), intent(out) :: candidate(:)
      real(dp), intent(in) :: scale
    end subroutine simple_proposal_interface

    subroutine simple_monitor_interface(x, f, iteration, evaluations, stop)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
    end subroutine simple_monitor_interface
  end interface

  public :: roptim_minimize
  public :: roptim_approximate_gradient
  public :: roptim_approximate_hessian
  public :: roptim_message

contains

  subroutine roptim_minimize(x, objective, result, method, gradient, lower, upper, &
                             control, user_data, proposal, monitor)
    real(dp), intent(inout) :: x(:)
    procedure(objective_interface) :: objective
    type(roptim_result_t), intent(out) :: result
    character(len=*), intent(in), optional :: method
    procedure(gradient_interface), optional :: gradient
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(roptim_control_t), intent(in), optional :: control
    class(*), intent(inout), optional :: user_data
    procedure(proposal_interface), optional :: proposal
    procedure(monitor_interface), optional :: monitor

    type(roptim_control_t) :: ctrl
    character(len=:), allocatable :: meth
    real(dp), allocatable :: scale(:), steps(:), z(:), zl(:), zu(:)
    real(dp) :: f
    integer :: n, status, iterations
    integer :: fncount, grcount
    character(len=:), allocatable :: message

    n = size(x)
    call initialize_result(result, n)
    if (n <= 0) then
      call set_result_failure(result, roptim_invalid_input, &
           'parameter vector must not be empty')
      return
    end if
    if (.not. all(ieee_is_finite(x))) then
      call set_result_failure(result, roptim_invalid_input, &
           'initial parameters must be finite')
      return
    end if

    ctrl = roptim_control_t()
    if (present(control)) ctrl = control
    meth = method_nelder_mead
    if (present(method)) meth = canonical_method(method)
    if (len(meth) == 0) then
      call set_result_failure(result, roptim_invalid_input, 'unknown optimization method')
      return
    end if

    if (.not. valid_control(ctrl, meth, n, result)) return

    allocate(scale(n), steps(n), z(n), zl(n), zu(n))
    scale = 1.0_dp
    steps = 1.0e-3_dp
    if (allocated(ctrl%parscale)) scale = ctrl%parscale
    if (allocated(ctrl%ndeps)) steps = ctrl%ndeps

    if (size(scale) /= n .or. any(.not. ieee_is_finite(scale)) .or. &
        any(scale <= 0.0_dp)) then
      call set_result_failure(result, roptim_invalid_input, &
           'parscale must have size(x), be finite, and be positive')
      return
    end if
    if (size(steps) /= n .or. any(.not. ieee_is_finite(steps)) .or. &
        any(steps <= 0.0_dp)) then
      call set_result_failure(result, roptim_invalid_input, &
           'ndeps must have size(x), be finite, and be positive')
      return
    end if

    call make_bounds(n, lower, upper, zl, zu, status, message)
    if (status /= roptim_success) then
      call set_result_failure(result, status, message)
      return
    end if
    zl = zl / scale
    zu = zu / scale
    z = x / scale

    if (meth == method_lbfgsb) then
      z = max(zl, min(zu, z))
    else if (present(lower) .or. present(upper)) then
      call set_result_failure(result, roptim_invalid_input, &
           'bounds are supported only by L-BFGS-B')
      return
    end if

    fncount = 0
    grcount = 0
    status = roptim_success
    iterations = 0
    message = 'converged'

    select case (meth)
    case (method_nelder_mead)
      call nelder_mead(z, f, objective_scaled, ctrl, fncount, iterations, status, &
                       message, monitor_scaled)
    case (method_bfgs)
      call bfgs_optimize(z, f, objective_scaled, gradient_scaled, ctrl, fncount, &
                         iterations, status, message, monitor_scaled)
    case (method_cg)
      call cg_optimize(z, f, objective_scaled, gradient_scaled, ctrl, fncount, &
                       iterations, status, message, monitor_scaled)
    case (method_lbfgsb)
      call lbfgsb_optimize(z, f, objective_scaled, gradient_scaled, zl, zu, ctrl, &
                           fncount, iterations, status, message, &
                           monitor_scaled)
    case (method_sann)
      if (ctrl%seed /= 0) call set_random_seed(ctrl%seed)
      call sann_optimize(z, f, objective_scaled, proposal_scaled, ctrl, fncount, &
                         iterations, status, message, monitor_scaled)
    case default
      status = roptim_invalid_input
      f = huge(1.0_dp)
      message = 'unknown optimization method'
    end select

    x = z * scale
    result%par = x
    result%value = f * ctrl%fnscale
    result%function_evaluations = fncount
    result%gradient_evaluations = grcount
    result%iterations = iterations
    result%convergence = status
    result%success = status == roptim_success
    result%message = message

    if (ctrl%compute_hessian .and. status /= roptim_invalid_input) then
      allocate(result%hessian(n, n))
      if (present(user_data)) then
        if (present(gradient)) then
          call roptim_approximate_hessian(x, objective, result%hessian, &
               gradient=gradient, ndeps=steps*scale, user_data=user_data)
        else
          call roptim_approximate_hessian(x, objective, result%hessian, &
               ndeps=steps*scale, user_data=user_data)
        end if
      else
        if (present(gradient)) then
          call roptim_approximate_hessian(x, objective, result%hessian, &
               gradient=gradient, ndeps=steps*scale)
        else
          call roptim_approximate_hessian(x, objective, result%hessian, &
               ndeps=steps*scale)
        end if
      end if
    end if

  contains

    function objective_scaled(v) result(value)
      real(dp), intent(in) :: v(:)
      real(dp) :: value
      real(dp), allocatable :: physical(:)

      physical = v * scale
      if (present(user_data)) then
        value = invoke_objective_with_data(objective, physical, user_data)
      else
        value = invoke_objective(objective, physical)
      end if
      fncount = fncount + 1
      value = value / ctrl%fnscale
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
    end function objective_scaled

    subroutine gradient_scaled(v, g)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: physical(:), gp(:)
      integer :: j
      real(dp) :: hp, hm, fp, fm
      real(dp), allocatable :: vp(:), vm(:)

      grcount = grcount + 1
      physical = v * scale
      allocate(gp(n))

      if (present(gradient)) then
        if (present(user_data)) then
          call invoke_gradient_with_data(gradient, physical, gp, user_data)
        else
          call invoke_gradient(gradient, physical, gp)
        end if
        g = gp * scale / ctrl%fnscale
      else
        allocate(vp(n), vm(n))
        do j = 1, n
          hp = min(steps(j), zu(j) - v(j))
          hm = min(steps(j), v(j) - zl(j))
          if (.not. ieee_is_finite(hp)) hp = steps(j)
          if (.not. ieee_is_finite(hm)) hm = steps(j)
          if (hp > 0.0_dp .and. hm > 0.0_dp) then
            vp = v
            vm = v
            vp(j) = v(j) + hp
            vm(j) = v(j) - hm
            fp = objective_scaled(vp)
            fm = objective_scaled(vm)
            g(j) = (fp - fm) / (hp + hm)
          else if (hp > 0.0_dp) then
            vp = v
            vp(j) = v(j) + hp
            fp = objective_scaled(vp)
            fm = objective_scaled(v)
            g(j) = (fp - fm) / hp
          else if (hm > 0.0_dp) then
            vm = v
            vm(j) = v(j) - hm
            fp = objective_scaled(v)
            fm = objective_scaled(vm)
            g(j) = (fp - fm) / hm
          else
            g(j) = 0.0_dp
          end if
        end do
      end if

      if (any(.not. ieee_is_finite(g))) g = huge(1.0_dp)
    end subroutine gradient_scaled

    subroutine proposal_scaled(current, candidate, proposal_scale)
      real(dp), intent(in) :: current(:)
      real(dp), intent(out) :: candidate(:)
      real(dp), intent(in) :: proposal_scale
      real(dp), allocatable :: pcurrent(:), pcandidate(:)

      if (present(proposal)) then
        pcurrent = current * scale
        allocate(pcandidate(n))
        if (present(user_data)) then
          call invoke_proposal_with_data(proposal, pcurrent, pcandidate, &
                                         proposal_scale, user_data)
        else
          call invoke_proposal(proposal, pcurrent, pcandidate, proposal_scale)
        end if
        candidate = pcandidate / scale
      else
        call gaussian_proposal(current, candidate, proposal_scale)
      end if
    end subroutine proposal_scaled

    subroutine monitor_scaled(v, value, iteration, evaluations, stop)
      real(dp), intent(in) :: v(:)
      real(dp), intent(in) :: value
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
      real(dp), allocatable :: physical(:)

      stop = .false.
      if (.not. present(monitor)) return
      physical = v * scale
      if (present(user_data)) then
        call invoke_monitor_with_data(monitor, physical, value*ctrl%fnscale, &
                                      iteration, evaluations, stop, user_data)
      else
        call invoke_monitor(monitor, physical, value*ctrl%fnscale, iteration, &
                            evaluations, stop)
      end if
    end subroutine monitor_scaled

  end subroutine roptim_minimize

  function invoke_objective(callback, x) result(f)
    procedure(objective_interface) :: callback
    real(dp), intent(in) :: x(:)
    real(dp) :: f

    f = callback(x)
  end function invoke_objective

  function invoke_objective_with_data(callback, x, user_data) result(f)
    procedure(objective_interface) :: callback
    real(dp), intent(in) :: x(:)
    class(*), intent(inout) :: user_data
    real(dp) :: f

    f = callback(x, user_data)
  end function invoke_objective_with_data

  subroutine invoke_gradient(callback, x, g)
    procedure(gradient_interface) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)

    call callback(x, g)
  end subroutine invoke_gradient

  subroutine invoke_gradient_with_data(callback, x, g, user_data)
    procedure(gradient_interface) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout) :: user_data

    call callback(x, g, user_data)
  end subroutine invoke_gradient_with_data

  subroutine invoke_proposal(callback, current, candidate, scale)
    procedure(proposal_interface) :: callback
    real(dp), intent(in) :: current(:)
    real(dp), intent(out) :: candidate(:)
    real(dp), intent(in) :: scale

    call callback(current, candidate, scale)
  end subroutine invoke_proposal

  subroutine invoke_proposal_with_data(callback, current, candidate, scale, &
                                       user_data)
    procedure(proposal_interface) :: callback
    real(dp), intent(in) :: current(:)
    real(dp), intent(out) :: candidate(:)
    real(dp), intent(in) :: scale
    class(*), intent(inout) :: user_data

    call callback(current, candidate, scale, user_data)
  end subroutine invoke_proposal_with_data

  subroutine invoke_monitor(callback, x, f, iteration, evaluations, stop)
    procedure(monitor_interface) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: f
    integer, intent(in) :: iteration, evaluations
    logical, intent(out) :: stop

    call callback(x, f, iteration, evaluations, stop)
  end subroutine invoke_monitor

  subroutine invoke_monitor_with_data(callback, x, f, iteration, evaluations, &
                                      stop, user_data)
    procedure(monitor_interface) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: f
    integer, intent(in) :: iteration, evaluations
    logical, intent(out) :: stop
    class(*), intent(inout) :: user_data

    call callback(x, f, iteration, evaluations, stop, user_data)
  end subroutine invoke_monitor_with_data

  subroutine nelder_mead(x, fbest, objective, ctrl, fncount, iterations, status, &
                         message, monitor)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fbest
    procedure(simple_objective_interface) :: objective
    type(roptim_control_t), intent(in) :: ctrl
    integer, intent(inout) :: fncount
    integer, intent(out) :: iterations, status
    character(len=:), allocatable, intent(out) :: message
    procedure(simple_monitor_interface) :: monitor

    integer :: n, j, maxit
    integer, allocatable :: order(:)
    real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: fr, fe, fc, threshold, fspread, xspread, step
    logical :: stop

    n = size(x)
    maxit = default_max_iterations(ctrl, method_nelder_mead)
    allocate(simplex(n, n+1), f(n+1), order(n+1), centroid(n), xr(n), xe(n), xc(n))

    simplex(:, 1) = x
    do j = 1, n
      simplex(:, j+1) = x
      step = 0.1_dp * max(abs(x(j)), 1.0_dp)
      simplex(j, j+1) = simplex(j, j+1) + step
    end do
    do j = 1, n+1
      f(j) = objective(simplex(:, j))
    end do

    status = roptim_iteration_limit
    message = 'maximum iterations reached'
    stop = .false.

    do iterations = 0, maxit
      call sort_indices(f, order)
      simplex = simplex(:, order)
      f = f(order)

      fbest = f(1)
      x = simplex(:, 1)
      fspread = maxval(abs(f - fbest))
      xspread = maxval(abs(simplex - spread(x, 2, n+1)))
      threshold = convergence_threshold(fbest, ctrl%abstol, ctrl%reltol)

      if (fspread <= threshold .and. xspread <= sqrt(max(threshold, epsilon(1.0_dp)))) then
        status = roptim_success
        message = 'converged'
        exit
      end if
      if (iterations >= maxit) exit

      centroid = sum(simplex(:, 1:n), dim=2) / real(n, dp)
      xr = centroid + ctrl%alpha * (centroid - simplex(:, n+1))
      fr = objective(xr)

      if (fr < f(1)) then
        xe = centroid + ctrl%gamma * (xr - centroid)
        fe = objective(xe)
        if (fe < fr) then
          simplex(:, n+1) = xe
          f(n+1) = fe
        else
          simplex(:, n+1) = xr
          f(n+1) = fr
        end if
      else if (fr < f(n)) then
        simplex(:, n+1) = xr
        f(n+1) = fr
      else
        if (fr < f(n+1)) then
          xc = centroid + ctrl%beta * (xr - centroid)
        else
          xc = centroid + ctrl%beta * (simplex(:, n+1) - centroid)
        end if
        fc = objective(xc)
        if (fc < min(fr, f(n+1))) then
          simplex(:, n+1) = xc
          f(n+1) = fc
        else
          do j = 2, n+1
            simplex(:, j) = simplex(:, 1) + ctrl%beta * &
                 (simplex(:, j) - simplex(:, 1))
            f(j) = objective(simplex(:, j))
          end do
        end if
      end if

      if (ctrl%trace > 0 .and. mod(iterations+1, max(1, ctrl%report)) == 0) then
        write(output_unit, '(a,i0,a,es16.8)') 'Nelder-Mead iteration ', &
             iterations+1, ': f=', minval(f)
      end if
      call monitor(x, fbest, iterations+1, fncount, stop)
      if (stop) then
        status = roptim_user_stop
        message = 'stopped by monitor callback'
        exit
      end if
    end do

    call sort_indices(f, order)
    x = simplex(:, order(1))
    fbest = f(order(1))
  end subroutine nelder_mead

  subroutine bfgs_optimize(x, f, objective, gradient, ctrl, fncount, &
                           iterations, status, message, monitor)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: f
    procedure(simple_objective_interface) :: objective
    procedure(simple_gradient_interface) :: gradient
    type(roptim_control_t), intent(in) :: ctrl
    integer, intent(inout) :: fncount
    integer, intent(out) :: iterations, status
    character(len=:), allocatable, intent(out) :: message
    procedure(simple_monitor_interface) :: monitor

    integer :: n, maxit, iter
    real(dp), allocatable :: h(:, :), g(:), gnew(:), p(:), xnew(:), s(:), y(:)
    real(dp), allocatable :: identity(:, :), a(:, :), tmp(:, :)
    real(dp) :: fnew, ys, rho, oldf
    logical :: accepted, stop

    n = size(x)
    maxit = default_max_iterations(ctrl, method_bfgs)
    allocate(h(n,n), g(n), gnew(n), p(n), xnew(n), s(n), y(n))
    allocate(identity(n,n), a(n,n), tmp(n,n))
    call identity_matrix(identity)
    h = identity

    f = objective(x)
    call gradient(x, g)
    if (.not. finite_state(f, g)) then
      status = roptim_nonfinite
      message = 'non-finite objective or gradient'
      iterations = 0
      return
    end if

    status = roptim_iteration_limit
    message = 'maximum iterations reached'
    stop = .false.

    iterations = 0
    do iter = 1, maxit
      if (maxval(abs(g)) <= max(ctrl%reltol, sqrt(epsilon(1.0_dp)))) then
        status = roptim_success
        message = 'converged'
        exit
      end if
      p = -matmul(h, g)
      if (dot_product(p, g) >= -epsilon(1.0_dp)*max(1.0_dp, dot_product(g,g))) then
        h = identity
        p = -g
      end if

      oldf = f
      call line_search(x, f, g, p, objective, gradient, ctrl, xnew, fnew, gnew, &
                       accepted)
      if (.not. accepted) then
        status = roptim_line_search_failure
        message = 'line search failed'
        exit
      end if

      s = xnew - x
      y = gnew - g
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp))*max(1.0_dp, norm2(s)*norm2(y))) then
        rho = 1.0_dp / ys
        a = identity - rho * outer_product(s, y)
        tmp = matmul(a, matmul(h, transpose(a)))
        h = tmp + rho * outer_product(s, s)
      else
        h = identity
      end if

      x = xnew
      f = fnew
      g = gnew

      iterations = iter
      if (ctrl%trace > 0 .and. mod(iter, max(1, ctrl%report)) == 0) then
        write(output_unit, '(a,i0,a,es16.8)') 'BFGS iteration ', iter, &
             ': f=', f
      end if
      call monitor(x, f, iter, fncount, stop)
      if (stop) then
        status = roptim_user_stop
        message = 'stopped by monitor callback'
        exit
      end if
      if ((abs(oldf-f) <= max(convergence_threshold(f, ctrl%abstol, ctrl%reltol), &
           epsilon(1.0_dp)) .or. norm2(s) <= sqrt(max(ctrl%reltol, epsilon(1.0_dp))) * &
           max(1.0_dp, norm2(x))) .and. maxval(abs(g)) <= &
           max(1.0e-6_dp, sqrt(max(ctrl%reltol, epsilon(1.0_dp))))) then
        status = roptim_success
        message = 'converged'
        exit
      end if
    end do
  end subroutine bfgs_optimize

  subroutine cg_optimize(x, f, objective, gradient, ctrl, fncount, &
                         iterations, status, message, monitor)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: f
    procedure(simple_objective_interface) :: objective
    procedure(simple_gradient_interface) :: gradient
    type(roptim_control_t), intent(in) :: ctrl
    integer, intent(inout) :: fncount
    integer, intent(out) :: iterations, status
    character(len=:), allocatable, intent(out) :: message
    procedure(simple_monitor_interface) :: monitor

    integer :: n, maxit, iter
    real(dp), allocatable :: g(:), gnew(:), d(:), xnew(:), delta(:), step_vector(:)
    real(dp) :: fnew, beta_cg, denominator, oldf
    logical :: accepted, stop

    n = size(x)
    maxit = default_max_iterations(ctrl, method_cg)
    allocate(g(n), gnew(n), d(n), xnew(n), delta(n), step_vector(n))

    f = objective(x)
    call gradient(x, g)
    d = -g
    status = roptim_iteration_limit
    message = 'maximum iterations reached'
    stop = .false.

    iterations = 0
    do iter = 1, maxit
      if (.not. finite_state(f, g)) then
        status = roptim_nonfinite
        message = 'non-finite objective or gradient'
        exit
      end if
      if (maxval(abs(g)) <= max(ctrl%reltol, sqrt(epsilon(1.0_dp)))) then
        status = roptim_success
        message = 'converged'
        exit
      end if
      if (dot_product(d, g) >= 0.0_dp) d = -g

      oldf = f
      call line_search(x, f, g, d, objective, gradient, ctrl, xnew, fnew, gnew, &
                       accepted)
      if (.not. accepted) then
        d = -g
        call line_search(x, f, g, d, objective, gradient, ctrl, xnew, fnew, &
                         gnew, accepted)
      end if
      if (.not. accepted) then
        status = roptim_line_search_failure
        message = 'line search failed'
        exit
      end if

      delta = gnew - g
      select case (ctrl%cg_type)
      case (1)
        denominator = max(dot_product(g, g), tiny(1.0_dp))
        beta_cg = dot_product(gnew, gnew) / denominator
      case (2)
        denominator = max(dot_product(g, g), tiny(1.0_dp))
        beta_cg = dot_product(gnew, delta) / denominator
      case (3)
        denominator = dot_product(d, delta)
        if (abs(denominator) <= tiny(1.0_dp)) then
          beta_cg = 0.0_dp
        else
          beta_cg = dot_product(gnew, delta) / denominator
        end if
      end select
      beta_cg = max(0.0_dp, beta_cg)
      if (mod(iter, n) == 0) beta_cg = 0.0_dp
      d = -gnew + beta_cg*d

      step_vector = xnew - x
      x = xnew
      f = fnew
      g = gnew
      iterations = iter

      if (ctrl%trace > 0 .and. mod(iter, max(1, ctrl%report)) == 0) then
        write(output_unit, '(a,i0,a,es16.8)') 'CG iteration ', iter, &
             ': f=', f
      end if
      call monitor(x, f, iter, fncount, stop)
      if (stop) then
        status = roptim_user_stop
        message = 'stopped by monitor callback'
        exit
      end if
      if ((abs(oldf-f) <= max(convergence_threshold(f, ctrl%abstol, ctrl%reltol), &
           epsilon(1.0_dp)) .or. norm2(step_vector) <= &
           sqrt(max(ctrl%reltol, epsilon(1.0_dp))) * max(1.0_dp, norm2(x))) .and. &
           maxval(abs(g)) <= max(1.0e-6_dp, &
           sqrt(max(ctrl%reltol, epsilon(1.0_dp))))) then
        status = roptim_success
        message = 'converged'
        exit
      end if
    end do
  end subroutine cg_optimize

  subroutine line_search(x, f, g, direction, objective, gradient, ctrl, xnew, &
                         fnew, gnew, accepted)
    real(dp), intent(in) :: x(:), f, g(:), direction(:)
    procedure(simple_objective_interface) :: objective
    procedure(simple_gradient_interface) :: gradient
    type(roptim_control_t), intent(in) :: ctrl
    real(dp), intent(out) :: xnew(:), fnew, gnew(:)
    logical, intent(out) :: accepted

    real(dp) :: step, slope, c1, c2, previous_f
    integer :: k

    c1 = 1.0e-4_dp
    c2 = 0.9_dp
    slope = dot_product(g, direction)
    accepted = .false.
    xnew = x
    fnew = f
    gnew = g
    if (slope >= 0.0_dp) return

    step = 1.0_dp
    previous_f = huge(1.0_dp)
    do k = 1, max(1, ctrl%max_line_search)
      xnew = x + step*direction
      fnew = objective(xnew)
      if (ieee_is_finite(fnew) .and. fnew <= f + c1*step*slope) then
        call gradient(xnew, gnew)
        if (all(ieee_is_finite(gnew))) then
          if (abs(dot_product(gnew, direction)) <= c2*abs(slope) .or. &
              step <= 1.0e-8_dp) then
            accepted = .true.
            return
          end if
          if (dot_product(gnew, direction) < 0.0_dp .and. fnew < previous_f) then
            previous_f = fnew
            step = min(2.0_dp*step, 8.0_dp)
            cycle
          end if
        end if
      end if
      previous_f = fnew
      step = 0.5_dp*step
      if (step < 1.0e-16_dp) exit
    end do

    if (ieee_is_finite(fnew) .and. fnew < f) then
      call gradient(xnew, gnew)
      accepted = all(ieee_is_finite(gnew))
    end if
  end subroutine line_search

  subroutine lbfgsb_optimize(x, f, objective, gradient, lower, upper, ctrl, &
                             fncount, iterations, status, message, monitor)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: f
    procedure(simple_objective_interface) :: objective
    procedure(simple_gradient_interface) :: gradient
    real(dp), intent(in) :: lower(:), upper(:)
    type(roptim_control_t), intent(in) :: ctrl
    integer, intent(inout) :: fncount
    integer, intent(out) :: iterations, status
    character(len=:), allocatable, intent(out) :: message
    procedure(simple_monitor_interface) :: monitor

    integer :: n, m, nwa, itask, iprint, icsave, maxit, i
    integer, allocatable :: nbd(:), iwa(:)
    integer :: lsave(4), isave(44)
    real(dp), allocatable :: wa(:), g(:), last_x(:)
    real(dp) :: dsave(29)
    logical :: stop

    n = size(x)
    m = ctrl%lbfgsb_memory
    maxit = default_max_iterations(ctrl, method_lbfgsb)
    allocate(nbd(n), iwa(3*n), g(n), last_x(n))
    nwa = 2*m*n + 11*m*m + 5*n + 8*m
    allocate(wa(nwa))

    do i = 1, n
      if (.not. ieee_is_finite(lower(i))) then
        if (.not. ieee_is_finite(upper(i))) then
          nbd(i) = 0
        else
          nbd(i) = 3
        end if
      else
        if (.not. ieee_is_finite(upper(i))) then
          nbd(i) = 1
        else
          nbd(i) = 2
        end if
      end if
    end do

    wa = 0.0_dp
    iwa = 0
    g = 0.0_dp
    lsave = 0
    isave = 0
    dsave = 0.0_dp
    icsave = 0
    itask = task_start
    iprint = -1
    call set_core_output(.false.)
    f = huge(1.0_dp)
    iterations = 0
    last_x = x
    status = roptim_iteration_limit
    message = 'maximum iterations reached'
    stop = .false.

    do
      call setulb(n, m, x, lower, upper, nbd, f, g, ctrl%factr, ctrl%pgtol, &
                  wa, iwa, itask, iprint, icsave, lsave, isave, dsave)
      select case (itask)
      case (task_fg, task_fg_line_search, task_fg_start)
        f = objective(x)
        call gradient(x, g)
        if (.not. finite_state(f, g)) then
          status = roptim_nonfinite
          message = 'non-finite objective or gradient'
          exit
        end if
      case (task_new_x)
        iterations = isave(30)
        if (ctrl%trace > 0 .and. mod(iterations, max(1, ctrl%report)) == 0) then
          write(output_unit, '(a,i0,a,es16.8,a,es16.8)') &
               'L-BFGS-B iteration ', iterations, ': f=', f, &
               ', projected gradient=', dsave(13)
        end if
        call monitor(x, f, iterations, fncount, stop)
        if (stop) then
          status = roptim_user_stop
          message = 'stopped by monitor callback'
          exit
        end if
        if (iterations >= maxit) then
          status = roptim_iteration_limit
          message = 'maximum iterations reached'
          exit
        end if
        if (maxval(abs(x-last_x)) <= ctrl%reltol*max(1.0_dp, maxval(abs(x)))) then
          status = roptim_success
          message = 'converged by parameter change'
          exit
        end if
        last_x = x
      case default
        if (itask == task_stop .or. itask >= 5) then
          if (itask == task_stop .or. (itask >= 6 .and. itask <= 9) .or. &
              itask == 22 .or. itask == 23 .or. itask == 24 .or. &
              itask == 25 .or. itask == 26) then
            status = roptim_success
            message = lbfgsb_task_message(itask)
          else
            status = roptim_line_search_failure
            message = lbfgsb_task_message(itask)
          end if
          exit
        end if
      end select
    end do

    iterations = max(iterations, isave(30))
  end subroutine lbfgsb_optimize

  subroutine sann_optimize(x, fbest, objective, proposal, ctrl, fncount, &
                           iterations, status, message, monitor)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fbest
    procedure(simple_objective_interface) :: objective
    procedure(simple_proposal_interface) :: proposal
    type(roptim_control_t), intent(in) :: ctrl
    integer, intent(inout) :: fncount
    integer, intent(out) :: iterations, status
    character(len=:), allocatable, intent(out) :: message
    procedure(simple_monitor_interface) :: monitor

    real(dp), allocatable :: current(:), candidate(:), best(:)
    real(dp) :: fcurrent, fcandidate, temperature, delta, u, proposal_scale
    integer :: k, stage, maxit
    logical :: stop

    maxit = default_max_iterations(ctrl, method_sann)
    allocate(current(size(x)), candidate(size(x)), best(size(x)))
    current = x
    best = x
    fcurrent = objective(current)
    fbest = fcurrent
    iterations = 1
    status = roptim_success
    message = 'completed simulated annealing schedule'
    stop = .false.
    stage = 1

    do while (iterations < maxit)
      temperature = ctrl%temperature / &
           log(real(stage, dp) + 1.7182818_dp)
      proposal_scale = temperature / ctrl%temperature
      do k = 1, ctrl%tmax
        if (iterations >= maxit) exit
        call proposal(current, candidate, proposal_scale)
        if (.not. all(ieee_is_finite(candidate))) then
          iterations = iterations + 1
          cycle
        end if
        fcandidate = objective(candidate)
        delta = fcandidate - fcurrent
        call random_number(u)
        if (delta <= 0.0_dp .or. u < exp(-delta/max(temperature, tiny(1.0_dp)))) then
          current = candidate
          fcurrent = fcandidate
          if (fcurrent <= fbest) then
            best = current
            fbest = fcurrent
          end if
        end if
        iterations = iterations + 1
      end do

      if (ctrl%trace > 0 .and. mod(stage, max(1, ctrl%report)) == 0) then
        write(output_unit, '(a,i0,a,es16.8)') 'SANN evaluation ', iterations, &
             ': best f=', fbest
      end if
      call monitor(best, fbest, iterations, fncount, stop)
      if (stop) then
        status = roptim_user_stop
        message = 'stopped by monitor callback'
        exit
      end if
      stage = stage + 1
    end do
    x = best
  end subroutine sann_optimize

  subroutine roptim_approximate_gradient(x, objective, gradient, ndeps, lower, &
                                         upper, user_data)
    real(dp), intent(in) :: x(:)
    procedure(objective_interface) :: objective
    real(dp), intent(out) :: gradient(:)
    real(dp), intent(in), optional :: ndeps(:), lower(:), upper(:)
    class(*), intent(inout), optional :: user_data

    real(dp), allocatable :: step(:), xp(:), xm(:), lo(:), hi(:)
    real(dp) :: fp, fm, hp, hm
    integer :: n, i, status
    character(len=:), allocatable :: message

    n = size(x)
    if (size(gradient) /= n) error stop 'gradient has wrong size'
    allocate(step(n), xp(n), xm(n), lo(n), hi(n))
    step = 1.0e-3_dp
    if (present(ndeps)) then
      if (size(ndeps) /= n) error stop 'ndeps has wrong size'
      step = ndeps
    end if
    call make_bounds(n, lower, upper, lo, hi, status, message)
    if (status /= roptim_success) error stop 'invalid derivative bounds'

    do i = 1, n
      hp = min(step(i), hi(i)-x(i))
      hm = min(step(i), x(i)-lo(i))
      if (.not. ieee_is_finite(hp)) hp = step(i)
      if (.not. ieee_is_finite(hm)) hm = step(i)
      if (hp > 0.0_dp .and. hm > 0.0_dp) then
        xp = x
        xm = x
        xp(i) = x(i) + hp
        xm(i) = x(i) - hm
        fp = call_objective(objective, xp, user_data)
        fm = call_objective(objective, xm, user_data)
        gradient(i) = (fp-fm)/(hp+hm)
      else if (hp > 0.0_dp) then
        xp = x
        xp(i) = x(i) + hp
        fp = call_objective(objective, xp, user_data)
        fm = call_objective(objective, x, user_data)
        gradient(i) = (fp-fm)/hp
      else if (hm > 0.0_dp) then
        xm = x
        xm(i) = x(i) - hm
        fp = call_objective(objective, x, user_data)
        fm = call_objective(objective, xm, user_data)
        gradient(i) = (fp-fm)/hm
      else
        gradient(i) = 0.0_dp
      end if
    end do
  end subroutine roptim_approximate_gradient

  subroutine roptim_approximate_hessian(x, objective, hessian, gradient, ndeps, &
                                        user_data)
    real(dp), intent(in) :: x(:)
    procedure(objective_interface) :: objective
    real(dp), intent(out) :: hessian(:, :)
    procedure(gradient_interface), optional :: gradient
    real(dp), intent(in), optional :: ndeps(:)
    class(*), intent(inout), optional :: user_data

    real(dp), allocatable :: step(:), xp(:), xm(:), gp(:), gm(:)
    integer :: n, i, j
    real(dp) :: h

    n = size(x)
    if (size(hessian,1) /= n .or. size(hessian,2) /= n) &
         error stop 'hessian has wrong shape'
    allocate(step(n), xp(n), xm(n), gp(n), gm(n))
    step = 1.0e-3_dp
    if (present(ndeps)) then
      if (size(ndeps) /= n) error stop 'ndeps has wrong size'
      step = ndeps
    end if

    do i = 1, n
      h = step(i)
      xp = x
      xm = x
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      if (present(gradient)) then
        call call_gradient(gradient, xp, gp, user_data)
        call call_gradient(gradient, xm, gm, user_data)
      else
        if (present(user_data)) then
          call roptim_approximate_gradient(xp, objective, gp, ndeps=step, &
                                           user_data=user_data)
          call roptim_approximate_gradient(xm, objective, gm, ndeps=step, &
                                           user_data=user_data)
        else
          call roptim_approximate_gradient(xp, objective, gp, ndeps=step)
          call roptim_approximate_gradient(xm, objective, gm, ndeps=step)
        end if
      end if
      hessian(:, i) = (gp-gm)/(2.0_dp*h)
    end do

    do i = 1, n
      do j = i+1, n
        h = 0.5_dp*(hessian(i,j)+hessian(j,i))
        hessian(i,j) = h
        hessian(j,i) = h
      end do
    end do
  end subroutine roptim_approximate_hessian

  function call_objective(objective, x, user_data) result(f)
    procedure(objective_interface) :: objective
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: f

    if (present(user_data)) then
      f = objective(x, user_data)
    else
      f = objective(x)
    end if
  end function call_objective

  subroutine call_gradient(gradient, x, g, user_data)
    procedure(gradient_interface) :: gradient
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data

    if (present(user_data)) then
      call gradient(x, g, user_data)
    else
      call gradient(x, g)
    end if
  end subroutine call_gradient

  subroutine gaussian_proposal(current, candidate, scale)
    real(dp), intent(in) :: current(:)
    real(dp), intent(out) :: candidate(:)
    real(dp), intent(in) :: scale

    integer :: i
    real(dp) :: u1, u2, radius, angle

    i = 1
    do while (i <= size(current))
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      radius = sqrt(-2.0_dp*log(u1))
      angle = 2.0_dp*acos(-1.0_dp)*u2
      candidate(i) = current(i) + scale*radius*cos(angle)
      if (i+1 <= size(current)) then
        candidate(i+1) = current(i+1) + scale*radius*sin(angle)
      end if
      i = i + 2
    end do
  end subroutine gaussian_proposal

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)
    integer(int64) :: state

    call random_seed(size=n)
    allocate(values(n))
    state = int(abs(seed), int64) + 104729_int64
    do i = 1, n
      state = modulo(1664525_int64*state + 1013904223_int64, 2147483647_int64)
      values(i) = int(max(1_int64, state))
    end do
    call random_seed(put=values)
  end subroutine set_random_seed

  subroutine make_bounds(n, lower, upper, lo, hi, status, message)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: lower(:), upper(:)
    real(dp), intent(out) :: lo(n), hi(n)
    integer, intent(out) :: status
    character(len=:), allocatable, intent(out) :: message

    lo = ieee_value(0.0_dp, ieee_negative_inf)
    hi = ieee_value(0.0_dp, ieee_positive_inf)
    status = roptim_success
    message = 'valid bounds'

    if (present(lower)) then
      if (size(lower) == 1) then
        lo = lower(1)
      else if (size(lower) == n) then
        lo = lower
      else
        status = roptim_invalid_input
        message = 'lower must have length one or size(x)'
        return
      end if
    end if
    if (present(upper)) then
      if (size(upper) == 1) then
        hi = upper(1)
      else if (size(upper) == n) then
        hi = upper
      else
        status = roptim_invalid_input
        message = 'upper must have length one or size(x)'
        return
      end if
    end if
    if (any(lo > hi)) then
      status = roptim_invalid_input
      message = 'each lower bound must be <= upper bound'
    end if
  end subroutine make_bounds

  logical function valid_control(ctrl, method, n, result)
    type(roptim_control_t), intent(in) :: ctrl
    character(len=*), intent(in) :: method
    integer, intent(in) :: n
    type(roptim_result_t), intent(inout) :: result

    valid_control = .false.
    if (.not. ieee_is_finite(ctrl%fnscale) .or. &
        abs(ctrl%fnscale) <= tiny(1.0_dp)) then
      call set_result_failure(result, roptim_invalid_input, &
           'fnscale must be finite and nonzero')
      return
    end if
    if (ctrl%reltol < 0.0_dp .or. .not. ieee_is_finite(ctrl%reltol)) then
      call set_result_failure(result, roptim_invalid_input, &
           'reltol must be finite and nonnegative')
      return
    end if
    if (ctrl%max_iterations < 0) then
      call set_result_failure(result, roptim_invalid_input, &
           'max_iterations must be nonnegative')
      return
    end if
    if (ctrl%report < 1 .and. ctrl%trace > 0) then
      call set_result_failure(result, roptim_invalid_input, &
           'report must be positive when trace is enabled')
      return
    end if
    if (method == method_nelder_mead) then
      if (ctrl%alpha <= 0.0_dp .or. ctrl%beta <= 0.0_dp .or. &
          ctrl%beta >= 1.0_dp .or. ctrl%gamma <= 1.0_dp) then
        call set_result_failure(result, roptim_invalid_input, &
             'invalid Nelder-Mead alpha, beta, or gamma')
        return
      end if
    end if
    if (method == method_cg .and. &
        (ctrl%cg_type < 1 .or. ctrl%cg_type > 3)) then
      call set_result_failure(result, roptim_invalid_input, &
           'cg_type must be 1, 2, or 3')
      return
    end if
    if (method == method_lbfgsb) then
      if (ctrl%lbfgsb_memory < 1 .or. ctrl%factr < 0.0_dp .or. &
          ctrl%pgtol < 0.0_dp) then
        call set_result_failure(result, roptim_invalid_input, &
             'invalid L-BFGS-B memory, factr, or pgtol')
        return
      end if
    end if
    if (method == method_sann) then
      if (ctrl%temperature <= 0.0_dp .or. ctrl%tmax < 1) then
        call set_result_failure(result, roptim_invalid_input, &
             'temperature and tmax must be positive')
        return
      end if
    end if
    if (n < 1) return
    valid_control = .true.
  end function valid_control

  integer function default_max_iterations(ctrl, method)
    type(roptim_control_t), intent(in) :: ctrl
    character(len=*), intent(in) :: method

    if (ctrl%max_iterations > 0) then
      default_max_iterations = ctrl%max_iterations
    else
      select case (method)
      case (method_nelder_mead)
        default_max_iterations = 500
      case (method_sann)
        default_max_iterations = 10000
      case default
        default_max_iterations = 100
      end select
    end if
  end function default_max_iterations

  real(dp) function convergence_threshold(f, abstol, reltol)
    real(dp), intent(in) :: f, abstol, reltol
    real(dp) :: at

    at = 0.0_dp
    if (ieee_is_finite(abstol) .and. abstol > 0.0_dp) at = abstol
    convergence_threshold = max(at, reltol*(abs(f)+reltol))
  end function convergence_threshold

  logical function finite_state(f, g)
    real(dp), intent(in) :: f, g(:)
    finite_state = ieee_is_finite(f) .and. all(ieee_is_finite(g))
  end function finite_state

  subroutine identity_matrix(a)
    real(dp), intent(out) :: a(:, :)
    integer :: i

    a = 0.0_dp
    do i = 1, min(size(a,1), size(a,2))
      a(i,i) = 1.0_dp
    end do
  end subroutine identity_matrix

  pure function outer_product(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))
    c = spread(a, 2, size(b))*spread(b, 1, size(a))
  end function outer_product

  subroutine sort_indices(values, order)
    real(dp), intent(in) :: values(:)
    integer, intent(out) :: order(:)
    integer :: i, j, key

    do i = 1, size(values)
      order(i) = i
    end do
    do i = 2, size(values)
      key = order(i)
      j = i-1
      do while (j >= 1)
        if (values(order(j)) <= values(key)) exit
        order(j+1) = order(j)
        j = j-1
      end do
      order(j+1) = key
    end do
  end subroutine sort_indices

  function canonical_method(method) result(canonical)
    character(len=*), intent(in) :: method
    character(len=:), allocatable :: canonical
    character(len=:), allocatable :: upper

    upper = uppercase(trim(adjustl(method)))
    select case (upper)
    case ('NELDER-MEAD', 'NELDER_MEAD', 'NM')
      canonical = method_nelder_mead
    case ('BFGS')
      canonical = method_bfgs
    case ('CG', 'CONJUGATE-GRADIENT', 'CONJUGATE_GRADIENT')
      canonical = method_cg
    case ('L-BFGS-B', 'LBFGSB', 'L_BFGS_B')
      canonical = method_lbfgsb
    case ('SANN', 'SIMULATED-ANNEALING', 'SIMULATED_ANNEALING')
      canonical = method_sann
    case default
      canonical = ''
    end select
  end function canonical_method

  pure function uppercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        out(i:i) = achar(code-32)
      else
        out(i:i) = text(i:i)
      end if
    end do
  end function uppercase

  subroutine initialize_result(result, n)
    type(roptim_result_t), intent(out) :: result
    integer, intent(in) :: n

    allocate(result%par(max(0,n)))
    result%par = 0.0_dp
    result%value = huge(1.0_dp)
    result%function_evaluations = 0
    result%gradient_evaluations = 0
    result%iterations = 0
    result%convergence = roptim_invalid_input
    result%success = .false.
    result%message = 'not started'
  end subroutine initialize_result

  subroutine set_result_failure(result, code, message)
    type(roptim_result_t), intent(inout) :: result
    integer, intent(in) :: code
    character(len=*), intent(in) :: message

    result%convergence = code
    result%success = .false.
    result%message = message
  end subroutine set_result_failure

  function roptim_message(code) result(message)
    integer, intent(in) :: code
    character(len=:), allocatable :: message

    select case (code)
    case (roptim_success)
      message = 'converged'
    case (roptim_iteration_limit)
      message = 'maximum iterations reached'
    case (roptim_invalid_input)
      message = 'invalid input'
    case (roptim_nonfinite)
      message = 'non-finite objective or gradient'
    case (roptim_line_search_failure)
      message = 'line search failure'
    case (roptim_user_stop)
      message = 'stopped by user callback'
    case default
      message = 'unknown status'
    end select
  end function roptim_message

  function lbfgsb_task_message(task) result(message)
    integer, intent(in) :: task
    character(len=:), allocatable :: message

    select case (task)
    case (3)
      message = 'L-BFGS-B stopped'
    case (5)
      message = 'L-BFGS-B abnormal termination'
    case (6)
      message = 'converged: relative reduction of objective'
    case (7)
      message = 'converged: projected gradient'
    case (8)
      message = 'converged'
    case (9)
      message = 'converged'
    case (10:19)
      message = 'L-BFGS-B warning'
    case (22:26)
      message = 'converged'
    case default
      block
        character(len=64) :: buffer
        write(buffer, '(a,i0)') 'L-BFGS-B task ', task
        message = trim(buffer)
      end block
    end select
  end function lbfgsb_task_message

end module roptim_mod
