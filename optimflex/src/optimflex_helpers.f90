! SPDX-License-Identifier: MIT
module optimflex_helpers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use optimflex_types, only : dp, optim_control, objective_fn, gradient_fn, hessian_fn, &
       residual_fn, jacobian_fn
  use optimflex_diff, only : fast_grad, fast_hess, fast_jac
  use optimflex_linalg, only : is_pd_fast, maxabs, symmetrize, solve_spd, solve_linear
  implicit none
  private
  public :: eval_objective, eval_gradient, eval_hessian, eval_jacobian
  public :: converged_basic, project_bounds, fill_bounds, projected_gradient
  public :: strong_wolfe, armijo_search, solve_with_ridge, dogleg_boundary_tau
  public :: start_timer, stop_timer

contains

  subroutine eval_objective(fn, x, f)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    f = fn(x)
  end subroutine eval_objective

  subroutine eval_gradient(fn, x, method, g, gradient)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: method
    real(dp), intent(out) :: g(:)
    procedure(gradient_fn), optional :: gradient
    if (present(gradient)) then
      call gradient(x, g)
    else
      call fast_grad(fn, x, g, method)
    end if
  end subroutine eval_gradient

  subroutine eval_hessian(fn, x, method, h, hessian)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: method
    real(dp), intent(out) :: h(:,:)
    procedure(hessian_fn), optional :: hessian
    if (present(hessian)) then
      call hessian(x, h)
    else
      call fast_hess(fn, x, h, method)
    end if
    call symmetrize(h)
  end subroutine eval_hessian

  subroutine eval_jacobian(residual, x, method, jac, jacobian)
    procedure(residual_fn) :: residual
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: method
    real(dp), intent(out) :: jac(:,:)
    procedure(jacobian_fn), optional :: jacobian
    real(dp), allocatable :: tmp(:,:)
    if (present(jacobian)) then
      tmp = jacobian(x)
      if (size(tmp,1) == size(jac,1) .and. size(tmp,2) == size(jac,2)) then
        jac = tmp
      else
        jac = 0.0_dp
      end if
    else
      call fast_jac(residual, x, jac, method)
    end if
  end subroutine eval_jacobian

  logical function converged_basic(c, ginf, f, fold, x, xold, iter, have_old) result(ok)
    type(optim_control), intent(in) :: c
    real(dp), intent(in) :: ginf, f, fold, x(:), xold(:)
    integer, intent(in) :: iter
    logical, intent(in) :: have_old
    real(dp) :: denom
    ok = .true.
    if (c%use_grad) ok = ok .and. (ginf <= c%tol_grad)
    if (c%use_abs_f .and. have_old) ok = ok .and. (abs(f-fold) <= c%tol_abs_f)
    if (c%use_rel_f .and. have_old) then
      denom = max(1.0_dp, abs(fold))
      ok = ok .and. (abs((f-fold)/denom) <= c%tol_rel_f)
    end if
    if (c%use_abs_x .and. iter > 1) ok = ok .and. (maxabs(x-xold) <= c%tol_abs_x)
    if (c%use_rel_x .and. iter > 1) then
      denom = max(1.0_dp, maxabs(xold))
      ok = ok .and. (maxabs(x-xold)/denom <= c%tol_rel_x)
    end if
  end function converged_basic

  subroutine fill_bounds(n, lower_in, upper_in, lower, upper)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: lower_in(:), upper_in(:)
    real(dp), intent(out) :: lower(n), upper(n)
    lower = -huge(1.0_dp)/16.0_dp
    upper = huge(1.0_dp)/16.0_dp
    if (present(lower_in)) then
      if (size(lower_in) == 1) then
        lower = lower_in(1)
      else if (size(lower_in) == n) then
        lower = lower_in
      end if
    end if
    if (present(upper_in)) then
      if (size(upper_in) == 1) then
        upper = upper_in(1)
      else if (size(upper_in) == n) then
        upper = upper_in
      end if
    end if
  end subroutine fill_bounds

  pure subroutine project_bounds(x, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    x = max(lower, min(x, upper))
  end subroutine project_bounds

  pure subroutine projected_gradient(x, g, lower, upper, epsb, pg)
    real(dp), intent(in) :: x(:), g(:), lower(:), upper(:), epsb
    real(dp), intent(out) :: pg(:)
    integer :: i
    pg = g
    do i = 1, size(x)
      if (x(i) <= lower(i)+epsb .and. g(i) > 0.0_dp) pg(i) = 0.0_dp
      if (x(i) >= upper(i)-epsb .and. g(i) < 0.0_dp) pg(i) = 0.0_dp
    end do
  end subroutine projected_gradient

  subroutine armijo_search(fn, x, f, g, p, c, alpha, xnew, fnew, ok, lower, upper)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:), f, g(:), p(:)
    type(optim_control), intent(in) :: c
    real(dp), intent(out) :: alpha, xnew(:), fnew
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: lower(:), upper(:)
    integer :: k
    real(dp) :: dphi
    dphi = dot_product(g,p)
    alpha = c%ls_alpha0
    ok = .false.
    do k = 1, c%ls_max_steps
      xnew = x + alpha*p
      if (present(lower) .and. present(upper)) call project_bounds(xnew,lower,upper)
      call eval_objective(fn, xnew, fnew)
      if (ieee_is_finite(fnew)) then
        if (fnew <= f + c%wolfe_c1*alpha*dphi) then
          ok = .true.
          return
        end if
      end if
      alpha = alpha*c%ls_shrink
      if (alpha < c%ls_min_step) exit
    end do
  end subroutine armijo_search

  subroutine strong_wolfe(fn, x, f, g, p, c, method, alpha, xnew, fnew, gnew, ok, gradient)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:), f, g(:), p(:)
    type(optim_control), intent(in) :: c
    integer, intent(in) :: method
    real(dp), intent(out) :: alpha, xnew(:), fnew, gnew(:)
    logical, intent(out) :: ok
    procedure(gradient_fn), optional :: gradient
    real(dp) :: aprev, fprev, a, dphi0, dphi, alo, ahi, flo, amid
    real(dp), allocatable :: xt(:), gt(:), xlo(:)
    integer :: k, z

    allocate(xt(size(x)), gt(size(x)), xlo(size(x)))
    dphi0 = dot_product(g,p)
    aprev = 0.0_dp
    fprev = f
    a = c%ls_alpha0
    ok = .false.
    do k = 1, c%ls_max_steps
      xt = x + a*p
      call eval_objective(fn,xt,fnew)
      if ((fnew > f + c%wolfe_c1*a*dphi0) .or. (k > 1 .and. fnew >= fprev)) then
        alo = aprev
        ahi = a
        flo = fprev
        xlo = x + alo*p
        do z = 1, c%zoom_max_steps
          amid = 0.5_dp*(alo+ahi)
          xt = x + amid*p
          call eval_objective(fn,xt,fnew)
          if ((fnew > f + c%wolfe_c1*amid*dphi0) .or. fnew >= flo) then
            ahi = amid
          else
            if (present(gradient)) then
              call eval_gradient(fn,xt,method,gt,gradient)
            else
              call eval_gradient(fn,xt,method,gt)
            end if
            dphi = dot_product(gt,p)
            if (abs(dphi) <= -c%wolfe_c2*dphi0) then
              alpha = amid
              xnew = xt
              gnew = gt
              ok = .true.
              return
            end if
            if (dphi*(ahi-alo) >= 0.0_dp) ahi = alo
            alo = amid
            flo = fnew
            xlo = xt
          end if
          if (abs(ahi-alo) < 1.0e-15_dp) exit
        end do
        return
      end if
      if (present(gradient)) then
        call eval_gradient(fn,xt,method,gt,gradient)
      else
        call eval_gradient(fn,xt,method,gt)
      end if
      dphi = dot_product(gt,p)
      if (abs(dphi) <= -c%wolfe_c2*dphi0) then
        alpha = a
        xnew = xt
        gnew = gt
        ok = .true.
        return
      end if
      if (dphi >= 0.0_dp) then
        alo = a
        ahi = aprev
        flo = fnew
        do z = 1, c%zoom_max_steps
          amid = 0.5_dp*(alo+ahi)
          xt = x + amid*p
          call eval_objective(fn,xt,fnew)
          if ((fnew > f + c%wolfe_c1*amid*dphi0) .or. fnew >= flo) then
            ahi = amid
          else
            if (present(gradient)) then
              call eval_gradient(fn,xt,method,gt,gradient)
            else
              call eval_gradient(fn,xt,method,gt)
            end if
            dphi = dot_product(gt,p)
            if (abs(dphi) <= -c%wolfe_c2*dphi0) then
              alpha = amid
              xnew = xt
              gnew = gt
              ok = .true.
              return
            end if
            if (dphi*(ahi-alo) >= 0.0_dp) ahi = alo
            alo = amid
            flo = fnew
          end if
        end do
        return
      end if
      aprev = a
      fprev = fnew
      a = 2.0_dp*a
    end do
  end subroutine strong_wolfe

  subroutine solve_with_ridge(h, rhs, c, step, hused, ok)
    real(dp), intent(in) :: h(:,:), rhs(:)
    type(optim_control), intent(in) :: c
    real(dp), intent(out) :: step(:), hused(:,:)
    logical, intent(out) :: ok
    real(dp) :: tau
    integer :: i, k

    hused = 0.5_dp*(h+transpose(h))
    call solve_spd(hused,rhs,step,ok)
    if (ok) return
    tau = max(c%ridge_offset,1.0e-12_dp)
    do k = 1, c%ridge_max_tries
      hused = 0.5_dp*(h+transpose(h))
      do i = 1, size(rhs)
        hused(i,i) = hused(i,i)+tau
      end do
      call solve_spd(hused,rhs,step,ok)
      if (ok) return
      tau = tau*c%ridge_mult
    end do
    call solve_linear(hused,rhs,step,ok)
  end subroutine solve_with_ridge

  pure real(dp) function dogleg_boundary_tau(p0, d, delta) result(tau)
    real(dp), intent(in) :: p0(:), d(:), delta
    real(dp) :: aa, bb, cc, disc
    aa = dot_product(d,d)
    bb = 2.0_dp*dot_product(p0,d)
    cc = dot_product(p0,p0)-delta*delta
    disc = max(0.0_dp,bb*bb-4.0_dp*aa*cc)
    tau = (-bb+sqrt(disc))/(2.0_dp*(aa+epsilon(1.0_dp)))
  end function dogleg_boundary_tau

  subroutine start_timer(cpu0, tick0, rate)
    real(dp), intent(out) :: cpu0
    integer, intent(out) :: tick0, rate
    call cpu_time(cpu0)
    call system_clock(tick0, rate)
  end subroutine start_timer

  subroutine stop_timer(cpu0, tick0, rate, cpu, elapsed)
    real(dp), intent(in) :: cpu0
    integer, intent(in) :: tick0, rate
    real(dp), intent(out) :: cpu, elapsed
    real(dp) :: c1
    integer :: t1
    call cpu_time(c1)
    call system_clock(t1)
    cpu = c1-cpu0
    if (rate > 0) then
      elapsed = real(t1-tick0,dp)/real(rate,dp)
    else
      elapsed = cpu
    end if
  end subroutine stop_timer

end module optimflex_helpers
