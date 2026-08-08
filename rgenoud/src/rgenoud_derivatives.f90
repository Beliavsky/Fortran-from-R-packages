! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_derivatives
  use rgenoud_kinds, only : dp
  use rgenoud_types, only : objective_fn, gradient_fn
  implicit none
  private
  public :: numerical_gradient, numerical_hessian, bfgs_optimize
contains
  subroutine numerical_gradient(fn, x, g, lower, upper, bounded)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp), intent(in), optional :: lower(:), upper(:)
    logical, intent(in), optional :: bounded
    real(dp) :: xp(size(x)), xm(size(x)), h, fp, fm
    logical :: use_bounds
    integer :: i

    use_bounds = .false.
    if (present(bounded)) use_bounds = bounded
    do i = 1, size(x)
      h = sqrt(epsilon(1.0_dp)) * max(abs(x(i)), 1.0_dp)
      xp = x
      xm = x
      xp(i) = x(i) + h
      xm(i) = x(i) - h
      if (use_bounds .and. present(lower) .and. present(upper)) then
        xp(i) = min(xp(i), upper(i))
        xm(i) = max(xm(i), lower(i))
      end if
      if (abs(xp(i) - xm(i)) <= tiny(1.0_dp)) then
        g(i) = 0.0_dp
      else
        fp = fn(xp)
        fm = fn(xm)
        g(i) = (fp - fm) / (xp(i) - xm(i))
      end if
    end do
  end subroutine numerical_gradient

  subroutine numerical_hessian(fn, x, hess, lower, upper, bounded)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(:, :)
    real(dp), intent(in), optional :: lower(:), upper(:)
    logical, intent(in), optional :: bounded
    real(dp) :: xp(size(x)), xm(size(x)), gp(size(x)), gm(size(x)), h
    logical :: use_bounds
    integer :: j

    use_bounds = .false.
    if (present(bounded)) use_bounds = bounded
    do j = 1, size(x)
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(abs(x(j)), 1.0_dp)
      xp = x
      xm = x
      xp(j) = x(j) + h
      xm(j) = x(j) - h
      if (use_bounds .and. present(lower) .and. present(upper)) then
        xp(j) = min(xp(j), upper(j))
        xm(j) = max(xm(j), lower(j))
      end if
      call numerical_gradient(fn, xp, gp, lower, upper, use_bounds)
      call numerical_gradient(fn, xm, gm, lower, upper, use_bounds)
      if (abs(xp(j) - xm(j)) <= tiny(1.0_dp)) then
        hess(:, j) = 0.0_dp
      else
        hess(:, j) = (gp - gm) / (xp(j) - xm(j))
      end if
    end do
    hess = 0.5_dp * (hess + transpose(hess))
  end subroutine numerical_hessian

  subroutine bfgs_optimize(fn, gradfn, x, maximize, lower, upper, bounded, &
      max_iter, gtol, f, iterations)
    procedure(objective_fn) :: fn
    procedure(gradient_fn) :: gradfn
    real(dp), intent(inout) :: x(:)
    logical, intent(in) :: maximize, bounded
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: max_iter
    real(dp), intent(in) :: gtol
    real(dp), intent(out) :: f
    integer, intent(out) :: iterations

    real(dp) :: hinv(size(x), size(x)), eye(size(x), size(x))
    real(dp) :: g(size(x)), gn(size(x)), p(size(x)), xn(size(x))
    real(dp) :: s(size(x)), y(size(x)), hy(size(x))
    real(dp) :: phi, phin, alpha, ys, yhy, signv
    integer :: i, iter, ls

    eye = 0.0_dp
    do i = 1, size(x)
      eye(i, i) = 1.0_dp
    end do
    hinv = eye
    signv = merge(-1.0_dp, 1.0_dp, maximize)
    f = fn(x)
    phi = signv * f
    call gradfn(x, g)
    g = signv * g
    iterations = 0

    do iter = 1, max_iter
      iterations = iter
      if (maxval(abs(g)) <= gtol) exit
      p = -matmul(hinv, g)
      if (dot_product(p, g) >= -epsilon(1.0_dp)) p = -g
      alpha = 1.0_dp
      do ls = 1, 40
        xn = x + alpha * p
        if (bounded) xn = min(max(xn, lower), upper)
        phin = signv * fn(xn)
        if (phin <= phi + 1.0e-4_dp * dot_product(g, xn - x)) exit
        alpha = 0.5_dp * alpha
      end do
      if (alpha <= 2.0_dp**(-40)) exit
      call gradfn(xn, gn)
      gn = signv * gn
      s = xn - x
      y = gn - g
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, norm2(s) * norm2(y))) then
        hy = matmul(hinv, y)
        yhy = dot_product(y, hy)
        hinv = hinv + ((ys + yhy) / (ys * ys)) * outer(s, s) &
          - (outer(hy, s) + outer(s, hy)) / ys
      else
        hinv = eye
      end if
      x = xn
      g = gn
      phi = phin
      f = signv * phi
    end do
    f = fn(x)
  contains
    pure function outer(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      integer :: ii
      do ii = 1, size(a)
        c(ii, :) = a(ii) * b
      end do
    end function outer
  end subroutine bfgs_optimize
end module rgenoud_derivatives
