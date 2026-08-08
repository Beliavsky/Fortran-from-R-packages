! SPDX-License-Identifier: MIT
module optimflex_diff
  use optimflex_types, only : dp, diff_forward, diff_central, diff_richardson, &
       objective_fn, residual_fn
  implicit none
  private
  public :: get_eps, fast_grad, fast_hess, fast_jac

contains

  pure real(dp) function get_eps(xi, kind_code) result(h)
    real(dp), intent(in) :: xi
    integer, intent(in) :: kind_code
    real(dp) :: scale
    scale = max(abs(xi), 1.0_dp)
    select case(kind_code)
    case(1)
      h = sqrt(epsilon(1.0_dp)) * scale
    case(2)
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp) * scale
    case default
      h = epsilon(1.0_dp)**0.25_dp * scale
    end select
  end function get_eps

  subroutine fast_grad(fn, x, g, method)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(in), optional :: method
    integer :: i, m
    real(dp) :: f0, fp, fm, h, g1, g2
    real(dp), allocatable :: z(:)

    m = diff_forward
    if (present(method)) m = method
    allocate(z(size(x)))
    if (m == diff_forward) then
      f0 = fn(x)
      do i = 1, size(x)
        h = get_eps(x(i), 1)
        z = x
        z(i) = x(i) + h
        g(i) = (fn(z) - f0) / h
      end do
    else if (m == diff_central) then
      do i = 1, size(x)
        h = get_eps(x(i), 2)
        z = x
        z(i) = x(i) + h
        fp = fn(z)
        z(i) = x(i) - h
        fm = fn(z)
        g(i) = (fp - fm) / (2.0_dp*h)
      end do
    else
      do i = 1, size(x)
        h = get_eps(x(i), 2)
        z = x
        z(i) = x(i) + h
        fp = fn(z)
        z(i) = x(i) - h
        fm = fn(z)
        g1 = (fp - fm) / (2.0_dp*h)
        h = 0.5_dp*h
        z = x
        z(i) = x(i) + h
        fp = fn(z)
        z(i) = x(i) - h
        fm = fn(z)
        g2 = (fp - fm) / (2.0_dp*h)
        g(i) = (4.0_dp*g2 - g1) / 3.0_dp
      end do
    end if
  end subroutine fast_grad

  subroutine hess_central(fn, x, hscale, hess)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:), hscale
    real(dp), intent(out) :: hess(:,:)
    integer :: i, j, n
    real(dp) :: f0, fpp, fpn, fnp, fnn, fp, fm, hi, hj
    real(dp), allocatable :: z(:)

    n = size(x)
    allocate(z(n))
    f0 = fn(x)
    hess = 0.0_dp
    do i = 1, n
      hi = hscale * get_eps(x(i), 3)
      z = x
      z(i) = x(i) + hi
      fp = fn(z)
      z(i) = x(i) - hi
      fm = fn(z)
      hess(i,i) = (fp - 2.0_dp*f0 + fm) / (hi*hi)
      do j = i + 1, n
        hj = hscale * get_eps(x(j), 3)
        z = x
        z(i) = x(i) + hi
        z(j) = x(j) + hj
        fpp = fn(z)
        z(j) = x(j) - hj
        fpn = fn(z)
        z(i) = x(i) - hi
        fnn = fn(z)
        z(j) = x(j) + hj
        fnp = fn(z)
        hess(i,j) = (fpp - fpn - fnp + fnn) / (4.0_dp*hi*hj)
        hess(j,i) = hess(i,j)
      end do
    end do
  end subroutine hess_central

  subroutine fast_hess(fn, x, hess, method)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(:,:)
    integer, intent(in), optional :: method
    integer :: i, j, n, m
    real(dp) :: f0, fi, fj, fij, hi, hj
    real(dp), allocatable :: z(:), fgr(:), h1(:,:), h2(:,:)

    m = diff_forward
    if (present(method)) m = method
    n = size(x)
    if (m == diff_forward) then
      allocate(z(n), fgr(n))
      f0 = fn(x)
      do i = 1, n
        hi = get_eps(x(i), 3)
        z = x
        z(i) = x(i) + hi
        fgr(i) = fn(z)
      end do
      hess = 0.0_dp
      do i = 1, n
        hi = get_eps(x(i), 3)
        z = x
        z(i) = x(i) + 2.0_dp*hi
        hess(i,i) = (fn(z) - 2.0_dp*fgr(i) + f0)/(hi*hi)
        do j = i + 1, n
          hj = get_eps(x(j), 3)
          z = x
          z(i) = x(i) + hi
          z(j) = x(j) + hj
          fij = fn(z)
          fi = fgr(i)
          fj = fgr(j)
          hess(i,j) = (fij-fi-fj+f0)/(hi*hj)
          hess(j,i) = hess(i,j)
        end do
      end do
    else if (m == diff_central) then
      call hess_central(fn, x, 1.0_dp, hess)
    else
      allocate(h1(n,n), h2(n,n))
      call hess_central(fn, x, 1.0_dp, h1)
      call hess_central(fn, x, 0.5_dp, h2)
      hess = (4.0_dp*h2 - h1)/3.0_dp
    end if
  end subroutine fast_hess

  subroutine jac_central(resfn, x, hscale, jac)
    procedure(residual_fn) :: resfn
    real(dp), intent(in) :: x(:), hscale
    real(dp), intent(out) :: jac(:,:)
    integer :: j
    real(dp) :: h
    real(dp), allocatable :: z(:), rp(:), rm(:)
    allocate(z(size(x)))
    do j = 1, size(x)
      h = hscale * get_eps(x(j), 2)
      z = x
      z(j) = x(j) + h
      rp = resfn(z)
      z(j) = x(j) - h
      rm = resfn(z)
      jac(:,j) = (rp-rm)/(2.0_dp*h)
    end do
  end subroutine jac_central

  subroutine fast_jac(resfn, x, jac, method)
    procedure(residual_fn) :: resfn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: jac(:,:)
    integer, intent(in), optional :: method
    integer :: j, m
    real(dp) :: h
    real(dp), allocatable :: z(:), r0(:), rp(:), j1(:,:), j2(:,:)

    m = diff_forward
    if (present(method)) m = method
    allocate(z(size(x)))
    if (m == diff_forward) then
      r0 = resfn(x)
      do j = 1, size(x)
        h = get_eps(x(j), 1)
        z = x
        z(j) = x(j) + h
        rp = resfn(z)
        jac(:,j) = (rp-r0)/h
      end do
    else if (m == diff_central) then
      call jac_central(resfn, x, 1.0_dp, jac)
    else
      allocate(j1(size(jac,1),size(jac,2)), j2(size(jac,1),size(jac,2)))
      call jac_central(resfn, x, 1.0_dp, j1)
      call jac_central(resfn, x, 0.5_dp, j2)
      jac = (4.0_dp*j2-j1)/3.0_dp
    end if
  end subroutine fast_jac

end module optimflex_diff
