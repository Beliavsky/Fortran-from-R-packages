! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_roots
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : root_func, scalar_func, root_result, steady_options
  use rootsolve_derivatives, only : perturb_value
  use rootsolve_linalg, only : solve_linear
  implicit none
  private
  public :: multiroot, uniroot_all, brent_root
contains

  function multiroot(f, start, options) result(res)
    procedure(root_func) :: f
    real(dp), intent(in) :: start(:)
    type(steady_options), intent(in), optional :: options
    type(root_result) :: res
    type(steady_options) :: opt
    real(dp), allocatable :: x(:), fx(:), fp(:), jac(:,:), delta(:), ewt(:), step(:)
    integer :: n, i, j, info
    real(dp) :: d, p

    if (present(options)) opt = options
    n = size(start)
    allocate(x(n), fx(n), fp(n), jac(n,n), delta(n), ewt(n), step(n))
    allocate(res%root(n), res%f_root(n), res%precision(opt%maxiter))
    res%precision = 0.0_dp
    x = start
    call f(x, fx)

    do i = 1, opt%maxiter
      res%iterations = i
      p = sum(abs(fx)) / real(max(1,n), dp)
      res%precision(i) = p
      call make_weights(x, opt, ewt)
      if (maxval(abs(fx) / ewt) < 1.0_dp) then
        res%converged = .true.
        exit
      end if

      do j = 1, n
        d = perturb_value(x(j), opt%pert)
        delta = x
        delta(j) = delta(j) + d
        call f(delta, fp)
        jac(:,j) = (fp - fx) / d
      end do
      call solve_linear(jac, -fx, step, info)
      if (info /= 0) then
        res%status = -2
        exit
      end if
      if (maxval(abs(step)) < opt%ctol) then
        x = x + step
        call enforce_positive(x, opt)
        call f(x, fx)
        res%converged = .true.
        exit
      end if
      x = x + step
      call enforce_positive(x, opt)
      call f(x, fx)
    end do

    res%root = x
    res%f_root = fx
    if (res%iterations > 0) res%estimated_precision = res%precision(res%iterations)
    if (res%converged) res%status = 1
    if (.not. res%converged .and. res%status == 0) res%status = -1
    if (res%iterations < size(res%precision)) then
      res%precision = res%precision(:res%iterations)
    end if
  end function multiroot

  function uniroot_all(f, lower, upper, ngrid, tol, maxiter) result(roots)
    procedure(scalar_func) :: f
    real(dp), intent(in) :: lower, upper
    integer, intent(in), optional :: ngrid, maxiter
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: roots(:)
    real(dp), allocatable :: x(:), fx(:), tmp(:)
    integer :: n, i, nr, mit
    real(dp) :: tt, r

    if (lower >= upper) error stop 'uniroot_all: lower must be less than upper'
    n = 100
    if (present(ngrid)) n = max(1, ngrid)
    tt = epsilon(1.0_dp)**0.2_dp
    if (present(tol)) tt = tol
    mit = 1000
    if (present(maxiter)) mit = maxiter
    allocate(x(n+1), fx(n+1), tmp(2*n+2))
    nr = 0
    do i = 1, n+1
      x(i) = lower + real(i-1,dp) * (upper-lower) / real(n,dp)
      fx(i) = f(x(i))
      if (abs(fx(i)) <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)))) then
        call append_unique(tmp, nr, x(i), 10.0_dp*tt)
      end if
    end do
    do i = 1, n
      if (fx(i) * fx(i+1) < 0.0_dp) then
        r = brent_root(f, x(i), x(i+1), tt, mit)
        call append_unique(tmp, nr, r, 10.0_dp*tt)
      end if
    end do
    allocate(roots(nr))
    if (nr > 0) roots = tmp(:nr)
  end function uniroot_all

  function brent_root(f, ax, bx, tol, maxiter) result(root)
    procedure(scalar_func) :: f
    real(dp), intent(in) :: ax, bx
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: root
    real(dp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, tt
    real(dp) :: min1, min2
    integer :: iter, mit

    a = ax
    b = bx
    fa = f(a)
    fb = f(b)
    if (fa * fb > 0.0_dp) error stop 'brent_root: root not bracketed'
    tt = epsilon(1.0_dp)**0.2_dp
    if (present(tol)) tt = tol
    mit = 1000
    if (present(maxiter)) mit = maxiter
    c = b
    fc = fb
    d = b-a
    e = d

    do iter = 1, mit
      if (fb*fc > 0.0_dp) then
        c = a
        fc = fa
        d = b-a
        e = d
      end if
      if (abs(fc) < abs(fb)) then
        a = b
        fa = fb
        b = c
        fb = fc
        c = a
        fc = fa
      end if
      tol1 = 2.0_dp*epsilon(1.0_dp)*abs(b) + 0.5_dp*tt
      xm = 0.5_dp*(c-b)
      if (abs(xm) <= tol1 .or. abs(fb) <= tiny(1.0_dp)) then
        root = b
        return
      end if
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (abs(a-c) <= epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(c))) then
          p = 2.0_dp*xm*s
          q = 1.0_dp-s
        else
          q = fa/fc
          r = fb/fc
          p = s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
          q = (q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
        end if
        if (p > 0.0_dp) q = -q
        p = abs(p)
        min1 = 3.0_dp*xm*q-abs(tol1*q)
        min2 = abs(e*q)
        if (2.0_dp*p < min(min1,min2)) then
          e = d
          d = p/q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      a = b
      fa = fb
      if (abs(d) > tol1) then
        b = b+d
      else
        b = b+sign(tol1,xm)
      end if
      fb = f(b)
    end do
    root = b
  end function brent_root

  subroutine append_unique(a, n, x, tol)
    real(dp), intent(inout) :: a(:)
    integer, intent(inout) :: n
    real(dp), intent(in) :: x, tol
    integer :: i
    do i = 1, n
      if (abs(a(i)-x) <= tol) return
    end do
    n = n+1
    a(n) = x
  end subroutine append_unique

  subroutine make_weights(x, opt, ewt)
    real(dp), intent(in) :: x(:)
    type(steady_options), intent(in) :: opt
    real(dp), intent(out) :: ewt(:)
    if (allocated(opt%rtol_vec)) then
      if (size(opt%rtol_vec) /= size(x)) error stop 'rtol_vec size mismatch'
      ewt = opt%rtol_vec * abs(x)
    else
      ewt = opt%rtol * abs(x)
    end if
    if (allocated(opt%atol_vec)) then
      if (size(opt%atol_vec) /= size(x)) error stop 'atol_vec size mismatch'
      ewt = ewt + opt%atol_vec
    else
      ewt = ewt + opt%atol
    end if
    ewt = max(ewt, tiny(1.0_dp))
  end subroutine make_weights

  subroutine enforce_positive(x, opt)
    real(dp), intent(inout) :: x(:)
    type(steady_options), intent(in) :: opt
    integer :: i
    if (opt%positive) x = max(x, 0.0_dp)
    if (allocated(opt%positive_index)) then
      do i = 1, size(opt%positive_index)
        if (opt%positive_index(i) >= 1 .and. opt%positive_index(i) <= size(x)) then
          x(opt%positive_index(i)) = max(0.0_dp, x(opt%positive_index(i)))
        end if
      end do
    end if
  end subroutine enforce_positive
end module rootsolve_roots
