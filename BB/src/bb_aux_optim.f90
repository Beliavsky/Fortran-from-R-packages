! SPDX-License-Identifier: GPL-3.0-only
module bb_aux_optim
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bb_kinds, only: dp
  use bb_interfaces, only: bb_vector_fn
  implicit none
  private
  public :: nelder_mead_residual, bfgs_residual

contains

  subroutine nelder_mead_residual(x, fn, maxit, feval, ok)
    real(dp), intent(inout) :: x(:)
    procedure(bb_vector_fn) :: fn
    integer, intent(in) :: maxit
    integer, intent(out) :: feval
    logical, intent(out) :: ok

    integer :: n, i, iter, ilo, ihi, inhi
    real(dp), allocatable :: simplex(:, :), vals(:), centroid(:), trial(:), fvec(:)
    real(dp) :: alpha, gamma, rho, sigma, fr, fe, fc, scale, spread

    n = size(x)
    feval = 0
    ok = .false.
    if (n < 1) return
    allocate(simplex(n, n + 1), vals(n + 1), centroid(n), trial(n), fvec(n))
    alpha = 1.0_dp
    gamma = 2.0_dp
    rho = 0.5_dp
    sigma = 0.5_dp

    simplex(:, 1) = x
    do i = 1, n
      simplex(:, i + 1) = x
      scale = 0.05_dp * max(abs(x(i)), 1.0_dp)
      simplex(i, i + 1) = simplex(i, i + 1) + scale
    end do
    do i = 1, n + 1
      vals(i) = residual_objective(simplex(:, i), fn, fvec)
      feval = feval + 1
      if (.not. ieee_is_finite(vals(i))) return
    end do

    do iter = 1, maxit
      call extrema(vals, ilo, ihi, inhi)
      spread = maxval(abs(vals - vals(ilo)))
      if (spread <= 1.0e-12_dp * (1.0_dp + abs(vals(ilo)))) exit

      centroid = 0.0_dp
      do i = 1, n + 1
        if (i /= ihi) centroid = centroid + simplex(:, i)
      end do
      centroid = centroid / real(n, dp)

      trial = centroid + alpha * (centroid - simplex(:, ihi))
      fr = residual_objective(trial, fn, fvec)
      feval = feval + 1
      if (.not. ieee_is_finite(fr)) return

      if (fr < vals(ilo)) then
        x = centroid + gamma * (trial - centroid)
        fe = residual_objective(x, fn, fvec)
        feval = feval + 1
        if (.not. ieee_is_finite(fe)) return
        if (fe < fr) then
          simplex(:, ihi) = x
          vals(ihi) = fe
        else
          simplex(:, ihi) = trial
          vals(ihi) = fr
        end if
      else if (fr < vals(inhi)) then
        simplex(:, ihi) = trial
        vals(ihi) = fr
      else
        if (fr < vals(ihi)) then
          x = centroid + rho * (trial - centroid)
        else
          x = centroid + rho * (simplex(:, ihi) - centroid)
        end if
        fc = residual_objective(x, fn, fvec)
        feval = feval + 1
        if (.not. ieee_is_finite(fc)) return
        if (fc < min(fr, vals(ihi))) then
          simplex(:, ihi) = x
          vals(ihi) = fc
        else
          do i = 1, n + 1
            if (i == ilo) cycle
            simplex(:, i) = simplex(:, ilo) + sigma * (simplex(:, i) - simplex(:, ilo))
            vals(i) = residual_objective(simplex(:, i), fn, fvec)
            feval = feval + 1
            if (.not. ieee_is_finite(vals(i))) return
          end do
        end if
      end if
    end do

    ilo = minloc(vals, dim=1)
    x = simplex(:, ilo)
    ok = .true.
  end subroutine nelder_mead_residual

  subroutine bfgs_residual(x, fn, maxit, feval, ok)
    real(dp), intent(inout) :: x(:)
    procedure(bb_vector_fn) :: fn
    integer, intent(in) :: maxit
    integer, intent(out) :: feval
    logical, intent(out) :: ok

    integer :: n, i, iter
    real(dp), allocatable :: h(:, :), g(:), gnew(:), p(:), xnew(:), s(:), y(:), fvec(:)
    real(dp), allocatable :: eye(:, :), v(:, :)
    real(dp) :: f, fnew, alpha, ys, eps, c1
    logical :: gok

    n = size(x)
    feval = 0
    ok = .false.
    if (n < 1) return
    allocate(h(n,n), g(n), gnew(n), p(n), xnew(n), s(n), y(n), fvec(n), eye(n,n), v(n,n))
    h = 0.0_dp
    eye = 0.0_dp
    do i = 1, n
      h(i,i) = 1.0_dp
      eye(i,i) = 1.0_dp
    end do
    eps = sqrt(epsilon(1.0_dp))
    c1 = 1.0e-4_dp

    f = residual_objective(x, fn, fvec)
    feval = feval + 1
    if (.not. ieee_is_finite(f)) return
    call residual_gradient(x, fn, f, eps, g, feval, gok)
    if (.not. gok) return

    do iter = 1, maxit
      if (maxval(abs(g)) <= 1.0e-8_dp) exit
      p = -matmul(h, g)
      if (dot_product(g, p) >= 0.0_dp) p = -g
      alpha = 1.0_dp
      do
        xnew = x + alpha * p
        fnew = residual_objective(xnew, fn, fvec)
        feval = feval + 1
        if (ieee_is_finite(fnew)) then
          if (fnew <= f + c1 * alpha * dot_product(g, p)) exit
        end if
        alpha = alpha * 0.5_dp
        if (alpha < 1.0e-12_dp) then
          xnew = x
          fnew = f
          exit
        end if
      end do
      if (all(xnew == x)) exit

      call residual_gradient(xnew, fn, fnew, eps, gnew, feval, gok)
      if (.not. gok) return
      s = xnew - x
      y = gnew - g
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * sqrt(dot_product(y,y) * dot_product(s,s))) then
        v = eye - outer_product(s, y) / ys
        h = matmul(matmul(v, h), transpose(v)) + outer_product(s, s) / ys
      else
        h = eye
      end if
      x = xnew
      f = fnew
      g = gnew
    end do
    ok = all(ieee_is_finite(x))
  end subroutine bfgs_residual

  real(dp) function residual_objective(x, fn, fvec) result(value)
    real(dp), intent(in) :: x(:)
    procedure(bb_vector_fn) :: fn
    real(dp), intent(out) :: fvec(:)
    call fn(x, fvec)
    if (.not. all(ieee_is_finite(fvec))) then
      value = huge(1.0_dp)
    else
      value = dot_product(fvec, fvec)
    end if
  end function residual_objective

  subroutine residual_gradient(x, fn, fbase, eps, g, feval, ok)
    real(dp), intent(in) :: x(:), fbase, eps
    procedure(bb_vector_fn) :: fn
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: feval
    logical, intent(out) :: ok
    real(dp) :: xp(size(x)), fv(size(x)), fp
    integer :: i

    ok = .false.
    do i = 1, size(x)
      xp = x
      xp(i) = xp(i) + eps * max(1.0_dp, abs(x(i)))
      fp = residual_objective(xp, fn, fv)
      feval = feval + 1
      if (.not. ieee_is_finite(fp)) return
      g(i) = (fp - fbase) / (xp(i) - x(i))
    end do
    ok = all(ieee_is_finite(g))
  end subroutine residual_gradient

  pure function outer_product(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))
    integer :: i
    do i = 1, size(a)
      c(i, :) = a(i) * b
    end do
  end function outer_product

  subroutine extrema(vals, ilo, ihi, inhi)
    real(dp), intent(in) :: vals(:)
    integer, intent(out) :: ilo, ihi, inhi
    integer :: i
    ilo = minloc(vals, dim=1)
    ihi = maxloc(vals, dim=1)
    inhi = merge(2, 1, ihi == 1)
    do i = 1, size(vals)
      if (i /= ihi .and. vals(i) > vals(inhi)) inhi = i
    end do
  end subroutine extrema

end module bb_aux_optim
