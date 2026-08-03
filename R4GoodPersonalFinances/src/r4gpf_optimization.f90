! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_optimization
  use r4gpf_kinds, only: dp
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_not_converged
  use r4gpf_linalg, only: vector_norm2
  implicit none
  private
  public :: scalar_function, vector_function
  public :: golden_section_minimize, bisection_root, nelder_mead_minimize

  abstract interface
    function scalar_function(x) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: value
    end function scalar_function

    function vector_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function vector_function
  end interface
contains

  subroutine golden_section_minimize(f, lower, upper, x_min, f_min, status, tolerance, max_iterations)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: lower, upper
    real(dp), intent(out) :: x_min, f_min
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp) :: a, b, c, d, fc, fd, tol, gr
    integer :: iter, maxit

    if (upper <= lower) then
      x_min = lower
      f_min = huge(1.0_dp)
      status = r4gpf_invalid_argument
      return
    end if
    tol = 1.0e-9_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 500
    if (present(max_iterations)) maxit = max(1, max_iterations)
    gr = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
    a = lower
    b = upper
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = f(c)
    fd = f(d)
    do iter = 1, maxit
      if (abs(b - a) <= tol * (1.0_dp + abs(a) + abs(b))) exit
      if (fc <= fd) then
        b = d
        d = c
        fd = fc
        c = b - gr * (b - a)
        fc = f(c)
      else
        a = c
        c = d
        fc = fd
        d = a + gr * (b - a)
        fd = f(d)
      end if
    end do
    x_min = 0.5_dp * (a + b)
    f_min = f(x_min)
    if (iter <= maxit) then
      status = r4gpf_success
    else
      status = r4gpf_not_converged
    end if
  end subroutine golden_section_minimize

  subroutine bisection_root(f, lower, upper, root, status, tolerance, max_iterations)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: lower, upper
    real(dp), intent(out) :: root
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp) :: a, b, fa, fb, fm, tol
    integer :: iter, maxit, expand

    tol = 1.0e-9_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 500
    if (present(max_iterations)) maxit = max(1, max_iterations)
    a = lower
    b = upper
    fa = f(a)
    fb = f(b)
    do expand = 1, 20
      if (fa * fb <= 0.0_dp) exit
      b = b + 1.5_dp * (b - a)
      fb = f(b)
    end do
    if (fa * fb > 0.0_dp) then
      root = 0.5_dp * (a + b)
      status = r4gpf_invalid_argument
      return
    end if
    do iter = 1, maxit
      root = 0.5_dp * (a + b)
      fm = f(root)
      if (abs(fm) <= tol .or. abs(b - a) <= tol * (1.0_dp + abs(root))) then
        status = r4gpf_success
        return
      end if
      if (fa * fm <= 0.0_dp) then
        b = root
        fb = fm
      else
        a = root
        fa = fm
      end if
    end do
    status = r4gpf_not_converged
  end subroutine bisection_root

  subroutine nelder_mead_minimize(f, x0, x_best, f_best, status, step, lower, upper, tolerance, max_iterations)
    procedure(vector_function) :: f
    real(dp), intent(in) :: x0(:)
    real(dp), allocatable, intent(out) :: x_best(:)
    real(dp), intent(out) :: f_best
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step(:), lower(:), upper(:), tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:), tmp(:)
    real(dp) :: alpha, gamma, rho, sigma, tol, spread, fr, fe, fc, s
    integer :: n, j, iter, maxit, ilo, ihi, inhi

    n = size(x0)
    allocate(x_best(n))
    if (n < 1) then
      f_best = huge(1.0_dp)
      status = r4gpf_invalid_argument
      return
    end if
    alpha = 1.0_dp
    gamma = 2.0_dp
    rho = 0.5_dp
    sigma = 0.5_dp
    tol = 1.0e-9_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 3000
    if (present(max_iterations)) maxit = max(1, max_iterations)
    allocate(simplex(n, n + 1), values(n + 1), centroid(n), xr(n), xe(n), xc(n), tmp(n))
    simplex(:, 1) = x0
    call clamp_point(simplex(:, 1))
    do j = 2, n + 1
      simplex(:, j) = x0
      if (present(step)) then
        s = step(j - 1)
      else
        s = 0.05_dp * max(1.0_dp, abs(x0(j - 1)))
      end if
      simplex(j - 1, j) = simplex(j - 1, j) + s
      call clamp_point(simplex(:, j))
    end do
    do j = 1, n + 1
      values(j) = f(simplex(:, j))
    end do

    do iter = 1, maxit
      call find_indices(values, ilo, ihi, inhi)
      spread = 0.0_dp
      do j = 1, n + 1
        spread = max(spread, vector_norm2(simplex(:, j) - simplex(:, ilo)))
      end do
      if (spread <= tol * (1.0_dp + vector_norm2(simplex(:, ilo)))) exit
      centroid = 0.0_dp
      do j = 1, n + 1
        if (j /= ihi) centroid = centroid + simplex(:, j)
      end do
      centroid = centroid / real(n, dp)
      xr = centroid + alpha * (centroid - simplex(:, ihi))
      call clamp_point(xr)
      fr = f(xr)
      if (fr < values(ilo)) then
        xe = centroid + gamma * (xr - centroid)
        call clamp_point(xe)
        fe = f(xe)
        if (fe < fr) then
          simplex(:, ihi) = xe
          values(ihi) = fe
        else
          simplex(:, ihi) = xr
          values(ihi) = fr
        end if
      else if (fr < values(inhi)) then
        simplex(:, ihi) = xr
        values(ihi) = fr
      else
        if (fr < values(ihi)) then
          xc = centroid + rho * (xr - centroid)
        else
          xc = centroid - rho * (centroid - simplex(:, ihi))
        end if
        call clamp_point(xc)
        fc = f(xc)
        if (fc < min(fr, values(ihi))) then
          simplex(:, ihi) = xc
          values(ihi) = fc
        else
          do j = 1, n + 1
            if (j == ilo) cycle
            simplex(:, j) = simplex(:, ilo) + sigma * (simplex(:, j) - simplex(:, ilo))
            call clamp_point(simplex(:, j))
            values(j) = f(simplex(:, j))
          end do
        end if
      end if
    end do
    call find_indices(values, ilo, ihi, inhi)
    x_best = simplex(:, ilo)
    f_best = values(ilo)
    if (iter <= maxit) then
      status = r4gpf_success
    else
      status = r4gpf_not_converged
    end if

  contains
    subroutine clamp_point(x)
      real(dp), intent(inout) :: x(:)
      if (present(lower)) x = max(x, lower)
      if (present(upper)) x = min(x, upper)
    end subroutine clamp_point

    subroutine find_indices(v, lo, hi, nhi)
      real(dp), intent(in) :: v(:)
      integer, intent(out) :: lo, hi, nhi
      integer :: k
      lo = 1
      hi = 1
      do k = 2, size(v)
        if (v(k) < v(lo)) lo = k
        if (v(k) > v(hi)) hi = k
      end do
      nhi = lo
      do k = 1, size(v)
        if (k == hi) cycle
        if (nhi == hi .or. v(k) > v(nhi)) nhi = k
      end do
    end subroutine find_indices
  end subroutine nelder_mead_minimize

end module r4gpf_optimization
