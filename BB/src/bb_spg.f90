! SPDX-License-Identifier: GPL-3.0-only
module bb_spg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bb_kinds, only: dp
  use bb_interfaces, only: bb_scalar_fn, bb_gradient_fn, bb_projection_fn
  use bb_projection, only: project_box, project_linear
  use bb_types, only: spg_control, spg_result, bb_success, bb_max_iterations, &
    bb_max_evaluations, bb_function_error, bb_gradient_error, &
    bb_projection_error, bb_invalid_input
  implicit none
  private

  public :: spg, spg_box, spg_projected, spg_linear

contains

  function spg(par, fn, control, gr) result(res)
    real(dp), intent(in) :: par(:)
    procedure(bb_scalar_fn) :: fn
    type(spg_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: res
    type(spg_control) :: ctrl

    ctrl = spg_control()
    if (present(control)) ctrl = control
    if (present(gr)) then
      res = spg_core(par, fn, identity_projection, .false., ctrl, gr)
    else
      res = spg_core(par, fn, identity_projection, .false., ctrl)
    end if
  end function spg

  function spg_box(par, fn, lower, upper, control, gr) result(res)
    real(dp), intent(in) :: par(:), lower(:), upper(:)
    procedure(bb_scalar_fn) :: fn
    type(spg_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: res
    type(spg_control) :: ctrl

    ctrl = spg_control()
    if (present(control)) ctrl = control
    if (size(lower) /= size(par) .or. size(upper) /= size(par)) then
      call set_spg_error(res, par, bb_invalid_input, 'Bounds must match parameter size.')
      return
    end if
    if (any(lower > upper)) then
      call set_spg_error(res, par, bb_invalid_input, 'Every lower bound must be <= upper bound.')
      return
    end if

    if (present(gr)) then
      res = spg_core(par, fn, box_projection, .true., ctrl, gr)
    else
      res = spg_core(par, fn, box_projection, .true., ctrl)
    end if
  contains
    subroutine box_projection(x, projected, ok)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: projected(:)
      logical, intent(out) :: ok
      call project_box(x, lower, upper, projected, ok)
    end subroutine box_projection
  end function spg_box

  function spg_projected(par, fn, project, control, gr) result(res)
    real(dp), intent(in) :: par(:)
    procedure(bb_scalar_fn) :: fn
    procedure(bb_projection_fn) :: project
    type(spg_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: res
    type(spg_control) :: ctrl

    ctrl = spg_control()
    if (present(control)) ctrl = control
    if (present(gr)) then
      res = spg_core(par, fn, project, .true., ctrl, gr)
    else
      res = spg_core(par, fn, project, .true., ctrl)
    end if
  end function spg_projected

  function spg_linear(par, fn, a, b, meq, control, gr, lower, upper) result(res)
    real(dp), intent(in) :: par(:)
    procedure(bb_scalar_fn) :: fn
    real(dp), intent(in) :: a(:, :), b(:)
    integer, intent(in) :: meq
    type(spg_control), intent(in), optional :: control
    procedure(bb_gradient_fn), optional :: gr
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(spg_result) :: res
    type(spg_control) :: ctrl
    real(dp), allocatable :: aa(:, :), bb(:)
    integer :: n, m, nlower, nupper, i, k

    ctrl = spg_control()
    if (present(control)) ctrl = control
    n = size(par)
    m = size(a, 1)
    if (size(a, 2) /= n .or. size(b) /= m .or. meq < 0 .or. meq > m) then
      call set_spg_error(res, par, bb_invalid_input, 'Invalid linear-constraint dimensions.')
      return
    end if
    if (present(lower) .neqv. present(upper)) then
      call set_spg_error(res, par, bb_invalid_input, 'Supply both lower and upper bounds or neither.')
      return
    end if

    if (present(lower)) then
      if (size(lower) /= n .or. size(upper) /= n .or. any(lower > upper)) then
        call set_spg_error(res, par, bb_invalid_input, 'Invalid bound dimensions or values.')
        return
      end if
      nlower = count(ieee_is_finite(lower))
      nupper = count(ieee_is_finite(upper))
      allocate(aa(m + nlower + nupper, n), bb(m + nlower + nupper))
      if (m > 0) then
        aa(1:m, :) = a
        bb(1:m) = b
      end if
      k = m
      do i = 1, n
        if (ieee_is_finite(lower(i))) then
          k = k + 1
          aa(k, :) = 0.0_dp
          aa(k, i) = 1.0_dp
          bb(k) = lower(i)
        end if
      end do
      do i = 1, n
        if (ieee_is_finite(upper(i))) then
          k = k + 1
          aa(k, :) = 0.0_dp
          aa(k, i) = -1.0_dp
          bb(k) = -upper(i)
        end if
      end do
    else
      allocate(aa(m, n), bb(m))
      aa = a
      bb = b
    end if

    if (present(gr)) then
      res = spg_core(par, fn, linear_projection, .true., ctrl, gr)
    else
      res = spg_core(par, fn, linear_projection, .true., ctrl)
    end if
  contains
    subroutine linear_projection(x, projected, ok)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: projected(:)
      logical, intent(out) :: ok
      call project_linear(x, aa, bb, meq, projected, ok)
    end subroutine linear_projection
  end function spg_linear

  function spg_core(par0, fn, project, do_project, ctrl, gr) result(res)
    real(dp), intent(in) :: par0(:)
    procedure(bb_scalar_fn) :: fn
    procedure(bb_projection_fn) :: project
    logical, intent(in) :: do_project
    type(spg_control), intent(in) :: ctrl
    procedure(bb_gradient_fn), optional :: gr
    type(spg_result) :: res

    real(dp), allocatable :: par(:), pnew(:), g(:), gnew(:), pg(:), d(:)
    real(dp), allocatable :: s(:), y(:), lastfv(:), work(:)
    real(dp) :: f, f0, fbest, fchg, pginfn, gbest, lambda
    real(dp) :: gtd, sts, yty, sty, lmin, lmax, fnew
    integer :: n, iter, feval, lsflag, idx
    logical :: ok

    n = size(par0)
    if (n < 1 .or. ctrl%m < 1 .or. ctrl%maxit < 0 .or. ctrl%maxfeval < 1 .or. &
        ctrl%method < 1 .or. ctrl%method > 3 .or. ctrl%eps <= 0.0_dp) then
      call set_spg_error(res, par0, bb_invalid_input, 'Invalid SPG control or empty parameter vector.')
      return
    end if

    allocate(par(n), pnew(n), g(n), gnew(n), pg(n), d(n), s(n), y(n), work(n))
    allocate(lastfv(ctrl%m))
    par = par0
    lmin = 1.0e-30_dp
    lmax = 1.0e30_dp
    iter = 0
    feval = 0
    fchg = huge(1.0_dp)
    lastfv = -1.0e99_dp

    ! The R implementation evaluates the objective before projecting the
    ! initial point.  We preserve that evaluation count and baseline.
    f = signed_objective(fn, par, ctrl%maximize)
    feval = 1
    if (.not. ieee_is_finite(f)) then
      call set_spg_error(res, par0, bb_function_error, 'Failure in initial function evaluation.')
      res%feval = feval
      return
    end if
    f0 = f

    if (do_project) then
      call project(par, work, ok)
      if (.not. ok .or. .not. all(ieee_is_finite(work))) then
        call set_spg_error(res, par0, bb_projection_error, 'Failure in projecting initial guess.')
        res%feval = feval
        return
      end if
      par = work
      ! R keeps the pre-projection f. This re-evaluation is intentionally not
      ! done to maintain BB's iteration semantics.
    end if

    if (present(gr)) then
      call gr(par, g)
      if (ctrl%maximize) g = -g
    else
      call numerical_gradient(fn, par, f, ctrl%maximize, ctrl%eps, g, feval, ok)
      if (.not. ok) then
        call set_spg_error(res, par, bb_gradient_error, 'Failure in initial numerical gradient evaluation.')
        res%feval = feval
        return
      end if
    end if
    if (.not. all(ieee_is_finite(g))) then
      call set_spg_error(res, par, bb_gradient_error, 'Failure in initial gradient evaluation.')
      res%feval = feval
      return
    end if

    fbest = f
    lastfv(1) = f
    pg = par - g
    if (do_project) then
      call project(pg, work, ok)
      if (.not. ok .or. .not. all(ieee_is_finite(work))) then
        call set_spg_error(res, par, bb_projection_error, 'Failure in initial projection.')
        res%feval = feval
        return
      end if
      pg = work
    end if
    pg = pg - par
    pginfn = maxval(abs(pg))
    gbest = pginfn
    lambda = 1.0_dp
    if (pginfn /= 0.0_dp) lambda = min(lmax, max(lmin, 1.0_dp / pginfn))

    if (ctrl%trace) write(*,'(a,i0,a,es16.8,a,es12.4)') &
      'iter: ', 0, ' f-value: ', display_value(f, ctrl%maximize), ' pgrad: ', pginfn

    lsflag = 0
    do while (pginfn > ctrl%gtol .and. iter <= ctrl%maxit .and. fchg > ctrl%ftol)
      iter = iter + 1
      d = par - lambda * g
      if (do_project) then
        call project(d, work, ok)
        if (.not. ok .or. .not. all(ieee_is_finite(work))) then
          lsflag = bb_projection_error
          exit
        end if
        d = work
      end if
      d = d - par
      gtd = dot_product(g, d)
      if (.not. ieee_is_finite(gtd)) then
        lsflag = bb_projection_error
        exit
      end if

      call nonmonotone_line_search(par, f, d, gtd, lastfv, fn, ctrl%maximize, &
        ctrl%maxfeval, feval, pnew, fnew, lsflag)
      if (lsflag /= 0) exit

      fchg = abs(f - fnew)
      f = fnew
      idx = mod(iter, ctrl%m) + 1
      lastfv(idx) = f

      if (present(gr)) then
        call gr(pnew, gnew)
        if (ctrl%maximize) gnew = -gnew
      else
        call numerical_gradient(fn, pnew, f, ctrl%maximize, ctrl%eps, gnew, feval, ok)
        if (.not. ok) then
          lsflag = bb_gradient_error
          exit
        end if
      end if
      if (.not. all(ieee_is_finite(gnew))) then
        lsflag = bb_gradient_error
        exit
      end if

      s = pnew - par
      y = gnew - g
      sts = dot_product(s, s)
      yty = dot_product(y, y)
      sty = dot_product(s, y)
      select case (ctrl%method)
      case (1)
        if (sts == 0.0_dp .or. sty < 0.0_dp) then
          lambda = lmax
        else
          lambda = min(lmax, max(lmin, sts / sty))
        end if
      case (2)
        if (sty < 0.0_dp .or. yty == 0.0_dp) then
          lambda = lmax
        else
          lambda = min(lmax, max(lmin, sty / yty))
        end if
      case (3)
        if (sts == 0.0_dp .or. yty == 0.0_dp) then
          lambda = lmax
        else
          lambda = min(lmax, max(lmin, sqrt(sts / yty)))
        end if
      end select

      par = pnew
      g = gnew
      pg = par - g
      if (do_project) then
        call project(pg, work, ok)
        if (.not. ok .or. .not. all(ieee_is_finite(work))) then
          lsflag = bb_projection_error
          exit
        end if
        pg = work
      end if
      pg = pg - par
      pginfn = maxval(abs(pg))

      if (ctrl%trace .and. mod(iter, max(1, ctrl%triter)) == 0) then
        write(*,'(a,i0,a,es16.8,a,es12.4)') 'iter: ', iter, ' f-value: ', &
          display_value(f, ctrl%maximize), ' pgrad: ', pginfn
      end if

      if (f < fbest) then
        fbest = f
        res%par = pnew
        gbest = pginfn
      end if
    end do

    if (.not. allocated(res%par)) res%par = par
    res%iter = iter
    res%feval = feval

    if (lsflag == 0) then
      if (pginfn <= ctrl%gtol .or. fchg <= ctrl%ftol) then
        res%convergence = bb_success
        res%message = 'Successful convergence'
      else if (iter > ctrl%maxit) then
        res%convergence = bb_max_iterations
        res%message = 'Maximum number of iterations exceeded'
      else
        res%convergence = bb_success
        res%message = 'Successful convergence'
      end if
      res%value = display_value(fbest, ctrl%maximize)
      res%gradient = gbest
    else
      select case (lsflag)
      case (bb_function_error)
        res%convergence = bb_function_error
        res%message = 'Failure: error in function evaluation'
      case (bb_max_evaluations)
        res%convergence = bb_max_evaluations
        res%message = 'Maximum function evaluations exceeded'
      case (bb_gradient_error)
        res%convergence = bb_gradient_error
        res%message = 'Failure: error in gradient evaluation'
      case default
        res%convergence = bb_projection_error
        res%message = 'Failure: error in projection or search direction'
      end select
      res%value = display_value(fbest, ctrl%maximize)
      res%gradient = gbest
    end if
    res%fn_reduction = display_value(f0 - fbest, ctrl%maximize)
  end function spg_core

  subroutine nonmonotone_line_search(p, f, d, gtd, lastfv, fn, maximize, &
      maxfeval, feval, pnew, fnew, flag)
    real(dp), intent(in) :: p(:), f, d(:), gtd, lastfv(:)
    procedure(bb_scalar_fn) :: fn
    logical, intent(in) :: maximize
    integer, intent(in) :: maxfeval
    integer, intent(inout) :: feval
    real(dp), intent(out) :: pnew(:), fnew
    integer, intent(out) :: flag
    real(dp) :: gamma, fmax, alpha, atemp, denom

    gamma = 1.0e-4_dp
    fmax = maxval(lastfv)
    alpha = 1.0_dp
    flag = 0

    do
      pnew = p + alpha * d
      fnew = signed_objective(fn, pnew, maximize)
      feval = feval + 1
      if (.not. ieee_is_finite(fnew)) then
        flag = bb_function_error
        return
      end if
      if (fnew <= fmax + gamma * alpha * gtd) return
      if (feval > maxfeval) then
        flag = bb_max_evaluations
        return
      end if

      if (alpha <= 0.1_dp) then
        alpha = alpha / 2.0_dp
      else
        denom = 2.0_dp * (fnew - f - alpha * gtd)
        if (denom == 0.0_dp) then
          atemp = alpha / 2.0_dp
        else
          atemp = -(gtd * alpha * alpha) / denom
        end if
        if (.not. ieee_is_finite(atemp) .or. atemp < 0.1_dp .or. &
            atemp > 0.9_dp * alpha) atemp = alpha / 2.0_dp
        alpha = atemp
      end if
    end do
  end subroutine nonmonotone_line_search

  subroutine numerical_gradient(fn, x, fbase, maximize, eps, g, feval, ok)
    procedure(bb_scalar_fn) :: fn
    real(dp), intent(in) :: x(:), fbase, eps
    logical, intent(in) :: maximize
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: feval
    logical, intent(out) :: ok
    real(dp) :: dx(size(x)), fx
    integer :: i

    ok = .false.
    do i = 1, size(x)
      dx = x
      dx(i) = dx(i) + eps
      fx = signed_objective(fn, dx, maximize)
      feval = feval + 1
      if (.not. ieee_is_finite(fx)) return
      g(i) = (fx - fbase) / eps
    end do
    ok = all(ieee_is_finite(g))
  end subroutine numerical_gradient

  real(dp) function signed_objective(fn, x, maximize) result(value)
    procedure(bb_scalar_fn) :: fn
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: maximize
    value = fn(x)
    if (maximize) value = -value
  end function signed_objective

  real(dp) function display_value(value, maximize) result(out)
    real(dp), intent(in) :: value
    logical, intent(in) :: maximize
    if (maximize) then
      out = -value
    else
      out = value
    end if
  end function display_value

  subroutine identity_projection(x, projected, ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: projected(:)
    logical, intent(out) :: ok
    projected = x
    ok = all(ieee_is_finite(x))
  end subroutine identity_projection

  subroutine set_spg_error(res, par, status, message)
    type(spg_result), intent(out) :: res
    real(dp), intent(in) :: par(:)
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    res%par = par
    res%convergence = status
    res%message = message
  end subroutine set_spg_error

end module bb_spg
