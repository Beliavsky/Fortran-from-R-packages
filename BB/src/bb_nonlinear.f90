! SPDX-License-Identifier: GPL-3.0-only
module bb_nonlinear
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bb_kinds, only: dp
  use bb_interfaces, only: bb_vector_fn
  use bb_types, only: sane_control, sane_result
  use bb_aux_optim, only: nelder_mead_residual, bfgs_residual
  implicit none
  private

  public :: sane, dfsane

contains

  function dfsane(par0, fn, control) result(res)
    real(dp), intent(in) :: par0(:)
    procedure(bb_vector_fn) :: fn
    type(sane_control), intent(in), optional :: control
    type(sane_result) :: res

    type(sane_control) :: ctrl
    real(dp), allocatable :: par(:), pnew(:), pbest(:), f(:), fnew(:), lastfv(:)
    real(dp) :: f0, normf, normf_best, fun, fune, alfa, alfa1, alfa2
    real(dp) :: eta, eps, pf, pp, ff, temp, normf_new
    integer :: n, fcnt, iter, bl, flag, knoimp, nm_eval, bfgs_eval
    logical :: ok

    ctrl = sane_control()
    if (present(control)) ctrl = control
    n = size(par0)
    if (.not. valid_sane_control(n, ctrl)) then
      call set_sane_error(res, par0, 2, 'Invalid DF-SANE control or empty parameter vector.')
      return
    end if

    allocate(par(n), pnew(n), pbest(n), f(n), fnew(n), lastfv(ctrl%m))
    par = par0
    fcnt = 0
    iter = 0
    bl = 0
    alfa = 1.0_dp
    alfa1 = 1.0_dp
    alfa2 = 1.0_dp
    eta = 1.0_dp
    eps = 1.0e-10_dp
    lastfv = 0.0_dp

    if (ctrl%nm) then
      call nelder_mead_residual(par, fn, 100, nm_eval, ok)
      fcnt = fcnt + nm_eval
      if (.not. ok) then
        call set_sane_error(res, par0, 3, 'Failure in Nelder-Mead start.')
        res%feval = fcnt
        return
      end if
    end if

    call fn(par, f)
    fcnt = fcnt + 1
    if (.not. all(ieee_is_finite(f))) then
      call set_sane_error(res, par, 3, 'Failure in initial functional evaluation.')
      res%feval = fcnt
      return
    end if

    f0 = sqrt(dot_product(f, f))
    normf = f0
    if (ctrl%trace) write(*,'(a,i0,a,es16.8)') &
      'Iteration: ', 0, ' ||F(x0)||: ', f0 / sqrt(real(n, dp))
    pbest = par
    normf_best = normf
    lastfv(1) = normf * normf
    flag = 0
    knoimp = 0

    do while (normf / sqrt(real(n, dp)) > ctrl%tol .and. iter <= ctrl%maxit)
      if (abs(alfa) <= eps .or. abs(alfa) >= 1.0_dp / eps) then
        if (normf > 1.0_dp) then
          alfa = 1.0_dp
        else if (normf >= 1.0e-5_dp) then
          alfa = 1.0_dp / normf
        else
          alfa = 1.0e5_dp
        end if
      end if

      if (iter == 0) then
        alfa = min(1.0_dp / max(normf, tiny(1.0_dp)), 1.0_dp)
        alfa1 = alfa
        alfa2 = alfa
      end if
      temp = alfa2
      alfa2 = alfa
      if (normf <= 0.01_dp) alfa = alfa1
      alfa1 = temp

      call dfsane_line_search(par, f, normf * normf, alfa, lastfv, eta, fn, &
        fcnt, bl, pnew, fnew, fune, flag)
      if (flag > 0) exit

      pf = dot_product(pnew - par, fnew - f)
      pp = dot_product(pnew - par, pnew - par)
      ff = dot_product(fnew - f, fnew - f)
      select case (ctrl%method)
      case (1)
        if (pf == 0.0_dp) then
          alfa = eps
        else
          alfa = pp / pf
        end if
      case (2)
        if (ff == 0.0_dp) then
          alfa = eps
        else
          alfa = pf / ff
        end if
      case (3)
        if (ff == 0.0_dp) then
          alfa = eps
        else
          alfa = sign(1.0_dp, pf) * sqrt(pp / ff)
        end if
      end select
      if (.not. ieee_is_finite(alfa)) alfa = eps

      par = pnew
      f = fnew
      fun = fune
      normf = sqrt(max(fun, 0.0_dp))
      if (normf < normf_best) then
        pbest = par
        normf_best = normf
        knoimp = 0
      else
        knoimp = knoimp + 1
      end if

      iter = iter + 1
      lastfv(1 + mod(iter, ctrl%m)) = fun
      eta = f0 / real(iter + 1, dp)**2
      if (ctrl%trace .and. mod(iter, max(1, ctrl%triter)) == 0) &
        write(*,'(a,i0,a,es16.8)') 'iteration: ', iter, ' ||F(xn)|| = ', normf

      if (knoimp == ctrl%noimp) then
        flag = 3
        exit
      end if
    end do

    if (flag == 0) then
      if (normf_best / sqrt(real(n, dp)) <= ctrl%tol) then
        res%convergence = 0
        res%message = 'Successful convergence'
      else if (iter > ctrl%maxit) then
        res%convergence = 1
        res%message = 'Maximum limit for iterations exceeded'
      else
        res%convergence = 2
        res%message = 'Method stagnated'
      end if
    else if (flag == 1) then
      res%convergence = 3
      res%message = 'Failure: error in function evaluation'
    else if (flag == 2) then
      res%convergence = 4
      res%message = 'Failure: maximum limit on steplength reductions exceeded'
    else
      res%convergence = 5
      res%message = 'Lack of improvement in objective function'
    end if

    if (ctrl%bfgs .and. (res%convergence == 2 .or. res%convergence == 5)) then
      par = pbest
      call bfgs_residual(par, fn, 200, bfgs_eval, ok)
      fcnt = fcnt + bfgs_eval
      if (ok) then
        call fn(par, fnew)
        fcnt = fcnt + 1
        if (all(ieee_is_finite(fnew))) then
          normf_new = sqrt(dot_product(fnew, fnew))
          if (normf_new < normf_best) then
            normf_best = normf_new
            pbest = par
          end if
        end if
      end if
      if (normf_best / sqrt(real(n, dp)) <= ctrl%tol) then
        res%convergence = 0
        res%message = 'Successful convergence'
      end if
    end if

    res%par = pbest
    res%residual = normf_best / sqrt(real(n, dp))
    res%fn_reduction = f0 - normf_best
    res%feval = fcnt
    res%iter = iter
  end function dfsane

  function sane(par0, fn, control) result(res)
    real(dp), intent(in) :: par0(:)
    procedure(bb_vector_fn) :: fn
    type(sane_control), intent(in), optional :: control
    type(sane_result) :: res

    type(sane_control) :: ctrl
    real(dp), allocatable :: par(:), pnew(:), pbest(:), f(:), fnew(:), fa(:), lastfv(:)
    real(dp) :: f0, normf, normf_best, fun, fune, alfa, alfa1, alfa2
    real(dp) :: eps, h, dg, sgn, lambda, temp, normf_new, denom
    integer :: n, fcnt, iter, bl, flag, knoimp, nm_eval, bfgs_eval
    logical :: ok

    ctrl = sane_control()
    if (present(control)) ctrl = control
    n = size(par0)
    if (.not. valid_sane_control(n, ctrl)) then
      call set_sane_error(res, par0, 2, 'Invalid SANE control or empty parameter vector.')
      return
    end if

    allocate(par(n), pnew(n), pbest(n), f(n), fnew(n), fa(n), lastfv(ctrl%m))
    par = par0
    fcnt = 0
    iter = 0
    bl = 0
    alfa = 1.0_dp
    alfa1 = 1.0_dp
    alfa2 = 1.0_dp
    eps = 1.0e-10_dp
    h = 1.0e-7_dp
    lastfv = 0.0_dp

    if (ctrl%nm) then
      call nelder_mead_residual(par, fn, 100, nm_eval, ok)
      fcnt = fcnt + nm_eval
      if (.not. ok) then
        call set_sane_error(res, par0, 2, 'Failure in Nelder-Mead start.')
        res%feval = fcnt
        return
      end if
    end if

    call fn(par, f)
    fcnt = fcnt + 1
    if (.not. all(ieee_is_finite(f))) then
      call set_sane_error(res, par, 2, 'Failure in initial functional evaluation.')
      res%feval = fcnt
      return
    end if

    f0 = sqrt(dot_product(f, f))
    normf = f0
    if (ctrl%trace) write(*,'(a,i0,a,es16.8)') &
      'Iteration: ', 0, ' ||F(x0)||: ', f0 / sqrt(real(n, dp))
    pbest = par
    normf_best = normf
    lastfv(1) = normf * normf
    flag = 0
    knoimp = 0

    do while (normf / sqrt(real(n, dp)) > ctrl%tol .and. iter <= ctrl%maxit)
      call fn(par + h * f, fa)
      fcnt = fcnt + 1
      if (.not. all(ieee_is_finite(fa))) then
        flag = 1
        exit
      end if
      dg = (dot_product(f, fa) - normf * normf) / h
      denom = normf * normf
      if (denom == 0.0_dp .or. .not. ieee_is_finite(dg)) then
        flag = 3
        exit
      end if
      if (abs(dg / denom) < eps) then
        flag = 3
        exit
      end if

      if (alfa <= eps .or. alfa >= 1.0_dp / eps) then
        if (normf > 1.0_dp) then
          alfa = 1.0_dp
        else if (normf >= 1.0e-5_dp) then
          alfa = normf
        else
          alfa = 1.0e-5_dp
        end if
      end if
      if (dg > 0.0_dp) then
        sgn = -1.0_dp
      else
        sgn = 1.0_dp
      end if

      if (iter == 0) then
        alfa = max(normf, 1.0_dp)
        alfa1 = alfa
        alfa2 = alfa
      end if
      temp = alfa2
      alfa2 = alfa
      if (normf <= 0.01_dp) alfa = alfa1
      alfa1 = temp

      lambda = 1.0_dp / alfa
      call sane_line_search(par, f, normf * normf, dg, lastfv, sgn, lambda, fn, &
        fcnt, bl, pnew, fnew, fune, flag)
      if (flag > 0) exit

      select case (ctrl%method)
      case (1)
        denom = lambda * dot_product(f, f)
        if (denom == 0.0_dp) then
          alfa = eps
        else
          alfa = dot_product(f, f - fnew) / denom
        end if
      case (2)
        denom = lambda * dot_product(f, f - fnew)
        if (denom == 0.0_dp) then
          alfa = eps
        else
          alfa = dot_product(f - fnew, f - fnew) / denom
        end if
      case (3)
        denom = lambda * lambda * dot_product(f, f)
        if (denom <= 0.0_dp) then
          alfa = eps
        else
          alfa = sqrt(dot_product(f - fnew, f - fnew) / denom)
        end if
      end select
      if (.not. ieee_is_finite(alfa)) alfa = eps

      par = pnew
      f = fnew
      fun = fune
      normf = sqrt(max(fun, 0.0_dp))
      if (normf < normf_best) then
        pbest = par
        normf_best = normf
        knoimp = 0
      else
        knoimp = knoimp + 1
      end if

      if (knoimp == ctrl%noimp) then
        flag = 4
        exit
      end if
      iter = iter + 1
      lastfv(1 + mod(iter, ctrl%m)) = fun
      if (ctrl%trace .and. mod(iter, max(1, ctrl%triter)) == 0) &
        write(*,'(a,i0,a,es16.8)') 'iteration: ', iter, ' ||F(xn)|| = ', normf
    end do

    if (flag == 0) then
      if (normf_best / sqrt(real(n, dp)) <= ctrl%tol) then
        res%convergence = 0
        res%message = 'Successful convergence'
      else if (iter > ctrl%maxit) then
        res%convergence = 1
        res%message = 'Maximum number of iterations exceeded'
      else
        res%convergence = 4
        res%message = 'Anomalous iteration'
      end if
    else if (flag == 1) then
      res%convergence = 2
      res%message = 'Error in function evaluation'
    else if (flag == 2) then
      res%convergence = 3
      res%message = 'Maximum limit on steplength reductions exceeded'
    else if (flag == 3) then
      res%convergence = 4
      res%message = 'Anomalous iteration'
    else
      res%convergence = 5
      res%message = 'Lack of improvement in objective function'
    end if

    if (ctrl%bfgs .and. (res%convergence == 4 .or. res%convergence == 5)) then
      par = pbest
      call bfgs_residual(par, fn, 200, bfgs_eval, ok)
      fcnt = fcnt + bfgs_eval
      if (ok) then
        call fn(par, fnew)
        fcnt = fcnt + 1
        if (all(ieee_is_finite(fnew))) then
          normf_new = sqrt(dot_product(fnew, fnew))
          if (normf_new < normf_best) then
            normf_best = normf_new
            pbest = par
          end if
        end if
      end if
      if (normf_best / sqrt(real(n, dp)) <= ctrl%tol) then
        res%convergence = 0
        res%message = 'Successful convergence'
      end if
    end if

    res%par = pbest
    res%residual = normf_best / sqrt(real(n, dp))
    res%fn_reduction = f0 - normf_best
    res%feval = fcnt
    res%iter = iter
  end function sane

  subroutine dfsane_line_search(x, f, fval, alfa, lastfv, eta, fn, fcnt, bl, &
      xnew, fnew, fune, flag)
    real(dp), intent(in) :: x(:), f(:), fval, alfa, lastfv(:), eta
    procedure(bb_vector_fn) :: fn
    integer, intent(inout) :: fcnt, bl
    real(dp), intent(out) :: xnew(:), fnew(:), fune
    integer, intent(out) :: flag

    integer :: cbl
    real(dp) :: gamma, sigma1, sigma2, lam1, lam2, fmax, fune1, fune2
    real(dp) :: lamc, c1, c2, d(size(x)), denom

    gamma = 1.0e-4_dp
    sigma1 = 0.1_dp
    sigma2 = 0.5_dp
    lam1 = 1.0_dp
    lam2 = 1.0_dp
    fmax = maxval(lastfv)
    d = -alfa * f
    flag = 2
    fune = huge(1.0_dp)

    do cbl = 0, 99
      xnew = x + lam1 * d
      call fn(xnew, fnew)
      fcnt = fcnt + 1
      if (.not. all(ieee_is_finite(fnew))) then
        flag = 1
        return
      end if
      fune1 = dot_product(fnew, fnew)
      if (fune1 <= fmax + eta - lam1 * lam1 * gamma * fval) then
        if (cbl >= 1) bl = bl + 1
        fune = fune1
        flag = 0
        return
      end if

      xnew = x - lam2 * d
      call fn(xnew, fnew)
      fcnt = fcnt + 1
      if (.not. all(ieee_is_finite(fnew))) then
        flag = 1
        return
      end if
      fune2 = dot_product(fnew, fnew)
      if (fune2 <= fmax + eta - lam2 * lam2 * gamma * fval) then
        if (cbl >= 1) bl = bl + 1
        fune = fune2
        flag = 0
        return
      end if

      denom = 2.0_dp * (fune1 + (2.0_dp * lam1 - 1.0_dp) * fval)
      if (denom == 0.0_dp) then
        lamc = sigma2 * lam1
      else
        lamc = 2.0_dp * fval * lam1 * lam1 / denom
      end if
      c1 = sigma1 * lam1
      c2 = sigma2 * lam1
      lam1 = min(c2, max(c1, lamc))

      denom = 2.0_dp * (fune2 + (2.0_dp * lam2 - 1.0_dp) * fval)
      if (denom == 0.0_dp) then
        lamc = sigma2 * lam2
      else
        lamc = 2.0_dp * fval * lam2 * lam2 / denom
      end if
      c1 = sigma1 * lam2
      c2 = sigma2 * lam2
      lam2 = min(c2, max(c1, lamc))
    end do
  end subroutine dfsane_line_search

  subroutine sane_line_search(x, f, fval, dg, lastfv, sgn, lambda, fn, fcnt, bl, &
      xnew, fnew, fune, flag)
    real(dp), intent(in) :: x(:), f(:), fval, dg, lastfv(:), sgn
    real(dp), intent(inout) :: lambda
    procedure(bb_vector_fn) :: fn
    integer, intent(inout) :: fcnt, bl
    real(dp), intent(out) :: xnew(:), fnew(:), fune
    integer, intent(out) :: flag

    integer :: cbl
    real(dp) :: gamma, sigma1, sigma2, fmax, gpd, lamc, c1, c2, denom

    gamma = 1.0e-4_dp
    sigma1 = 0.1_dp
    sigma2 = 0.5_dp
    fmax = maxval(lastfv)
    gpd = -2.0_dp * abs(dg)
    flag = 2
    fune = huge(1.0_dp)

    do cbl = 0, 99
      xnew = x + lambda * sgn * f
      call fn(xnew, fnew)
      fcnt = fcnt + 1
      if (.not. all(ieee_is_finite(fnew))) then
        flag = 1
        return
      end if
      fune = dot_product(fnew, fnew)
      if (fune <= fmax + lambda * gpd * gamma) then
        if (cbl >= 1) bl = bl + 1
        flag = 0
        return
      end if

      denom = 2.0_dp * (fune - fval - lambda * gpd)
      if (denom == 0.0_dp) then
        lamc = sigma2 * lambda
      else
        lamc = -(gpd * lambda * lambda) / denom
      end if
      c1 = sigma1 * lambda
      c2 = sigma2 * lambda
      lambda = min(c2, max(c1, lamc))
    end do
  end subroutine sane_line_search

  logical function valid_sane_control(n, ctrl) result(ok)
    integer, intent(in) :: n
    type(sane_control), intent(in) :: ctrl
    ok = n > 0 .and. ctrl%m > 0 .and. ctrl%maxit >= 0 .and. ctrl%method >= 1 .and. &
      ctrl%method <= 3 .and. ctrl%tol >= 0.0_dp .and. ctrl%noimp > 0
  end function valid_sane_control

  subroutine set_sane_error(res, par, status, message)
    type(sane_result), intent(out) :: res
    real(dp), intent(in) :: par(:)
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    res%par = par
    res%convergence = status
    res%message = message
  end subroutine set_sane_error

end module bb_nonlinear
