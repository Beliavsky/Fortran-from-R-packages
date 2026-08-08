! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_dense
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result, rq_multi_result
  use quantreg_linalg, only : spd_solve
  implicit none
  private
  public :: rq_fit_fnb, rq_fit_fnc, rq_fit_qfnb, rq_wfit_fnb
  public :: rq_fit_lasso, rq_fit_scad, check_loss, check_loss_sum
contains

  pure elemental real(dp) function check_loss(r, tau) result(v)
    real(dp), intent(in) :: r, tau
    if (r >= 0.0_dp) then
      v = tau * r
    else
      v = (tau - 1.0_dp) * r
    end if
  end function check_loss

  pure real(dp) function check_loss_sum(r, tau) result(v)
    real(dp), intent(in) :: r(:), tau
    integer :: i
    v = 0.0_dp
    do i = 1, size(r)
      v = v + check_loss(r(i), tau)
    end do
  end function check_loss_sum

  subroutine rq_fit_fnb(x, y, tau, result, beta, eps)
    real(dp), intent(in) :: x(:,:), y(:), tau
    type(rq_result), intent(out) :: result
    real(dp), intent(in), optional :: beta, eps
    real(dp) :: bta, tol
    real(dp), allocatable :: rhs(:), xvar(:), upper(:), coef_work(:), dual_obs(:)
    integer :: nit(3), info
    integer :: n, p, j

    n = size(x,1)
    p = size(x,2)
    result%tau = tau
    result%info = 0
    if (size(y) /= n .or. tau <= 0.0_dp .or. tau >= 1.0_dp) then
      result%info = -1
      return
    end if
    bta = 0.99995_dp
    if (present(beta)) bta = beta
    tol = 1.0e-6_dp
    if (present(eps)) tol = eps
    if (tau < tol .or. tau > 1.0_dp - tol) then
      result%info = -2
      return
    end if

    allocate(rhs(p), xvar(n), upper(n), coef_work(p), dual_obs(n))
    do j = 1, p
      rhs(j) = (1.0_dp - tau) * sum(x(:,j))
    end do
    xvar = 1.0_dp - tau
    upper = 1.0_dp
    call fnb_core(x, -y, rhs, xvar, upper, bta, tol, coef_work, dual_obs, nit, info)

    allocate(result%coefficients(p), result%residuals(n), result%dual(n))
    result%coefficients = -coef_work
    result%residuals = y - matmul(x, result%coefficients)
    result%dual = dual_obs
    result%iterations = nit(1)
    result%corrector_steps = nit(2)
    result%info = info
  end subroutine rq_fit_fnb

  subroutine rq_wfit_fnb(x, y, weights, tau, result, beta, eps)
    real(dp), intent(in) :: x(:,:), y(:), weights(:), tau
    type(rq_result), intent(out) :: result
    real(dp), intent(in), optional :: beta, eps
    real(dp), allocatable :: wx(:,:), wy(:)
    integer :: j

    if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) then
      result%info = -1
      return
    end if
    allocate(wx(size(x,1),size(x,2)), wy(size(y)))
    wy = y * weights
    do j = 1, size(x,2)
      wx(:,j) = x(:,j) * weights
    end do
    if (present(beta) .and. present(eps)) then
      call rq_fit_fnb(wx, wy, tau, result, beta, eps)
    else if (present(beta)) then
      call rq_fit_fnb(wx, wy, tau, result, beta=beta)
    else if (present(eps)) then
      call rq_fit_fnb(wx, wy, tau, result, eps=eps)
    else
      call rq_fit_fnb(wx, wy, tau, result)
    end if
    if (result%info == 0) then
      result%residuals = y - matmul(x, result%coefficients)
    end if
  end subroutine rq_wfit_fnb

  subroutine rq_fit_qfnb(x, y, taus, result, beta, eps)
    real(dp), intent(in) :: x(:,:), y(:), taus(:)
    type(rq_multi_result), intent(out) :: result
    real(dp), intent(in), optional :: beta, eps
    type(rq_result) :: fit
    integer :: k, p

    p = size(x,2)
    allocate(result%coefficients(p,size(taus)), result%tau(size(taus)))
    result%tau = taus
    result%info = 0
    do k = 1, size(taus)
      if (present(beta) .and. present(eps)) then
        call rq_fit_fnb(x, y, taus(k), fit, beta, eps)
      else if (present(beta)) then
        call rq_fit_fnb(x, y, taus(k), fit, beta=beta)
      else if (present(eps)) then
        call rq_fit_fnb(x, y, taus(k), fit, eps=eps)
      else
        call rq_fit_fnb(x, y, taus(k), fit)
      end if
      if (fit%info /= 0) then
        result%info = fit%info
        return
      end if
      result%coefficients(:,k) = fit%coefficients
      result%iterations = result%iterations + fit%iterations
      result%corrector_steps = result%corrector_steps + fit%corrector_steps
    end do
  end subroutine rq_fit_qfnb

  subroutine rq_fit_fnc(x, y, rmat, rvec, tau, result, beta, eps)
    real(dp), intent(in) :: x(:,:), y(:), rmat(:,:), rvec(:), tau
    type(rq_result), intent(out) :: result
    real(dp), intent(in), optional :: beta, eps
    real(dp) :: bta, tol
    real(dp), allocatable :: rhs(:), x1(:), x2(:), coef_work(:), dual_obs(:)
    integer :: nit(3), info, j, p, n1, n2

    n1 = size(x,1)
    p = size(x,2)
    n2 = size(rmat,1)
    result%tau = tau
    if (size(y) /= n1 .or. size(rvec) /= n2 .or. size(rmat,2) /= p) then
      result%info = -1
      return
    end if
    bta = 0.9995_dp
    if (present(beta)) bta = beta
    tol = 1.0e-6_dp
    if (present(eps)) tol = eps
    if (tau < tol .or. tau > 1.0_dp - tol) then
      result%info = -2
      return
    end if

    allocate(rhs(p), x1(n1), x2(n2), coef_work(p), dual_obs(n1))
    do j = 1, p
      rhs(j) = (1.0_dp - tau) * sum(x(:,j))
    end do
    x1 = 1.0_dp - tau
    x2 = 1.0_dp
    call fnc_core(x, -y, rmat, -rvec, rhs, x1, x2, bta, tol, coef_work, dual_obs, nit, info)

    allocate(result%coefficients(p), result%residuals(n1), result%dual(n1))
    result%coefficients = -coef_work
    result%residuals = y - matmul(x, result%coefficients)
    result%dual = dual_obs
    result%iterations = nit(1)
    result%corrector_steps = nit(2)
    result%info = info
  end subroutine rq_fit_fnc

  subroutine rq_fit_lasso(x, y, tau, lambda, result, beta, eps)
    real(dp), intent(in) :: x(:,:), y(:), tau
    real(dp), intent(in) :: lambda(:)
    type(rq_result), intent(out) :: result
    real(dp), intent(in), optional :: beta, eps
    real(dp), allocatable :: xa(:,:), ya(:), rhs(:), xv(:), upper(:), workcoef(:), dual(:)
    integer, allocatable :: keep(:)
    integer :: p, n, m, j, k, nit(3), info
    real(dp) :: bta, tol

    n = size(x,1)
    p = size(x,2)
    if (size(lambda) /= p .or. any(lambda < 0.0_dp)) then
      result%info = -1
      return
    end if
    m = count(lambda > 0.0_dp)
    allocate(keep(m))
    k = 0
    do j = 1, p
      if (lambda(j) > 0.0_dp) then
        k = k + 1
        keep(k) = j
      end if
    end do
    allocate(xa(n+m,p), ya(n+m), rhs(p), xv(n+m), upper(n+m), workcoef(p), dual(n+m))
    xa(1:n,:) = x
    ya(1:n) = y
    xa(n+1:n+m,:) = 0.0_dp
    ya(n+1:n+m) = 0.0_dp
    do k = 1, m
      xa(n+k,keep(k)) = lambda(keep(k))
    end do
    do j = 1, p
      rhs(j) = (1.0_dp - tau) * sum(x(:,j)) + 0.5_dp * sum(xa(n+1:n+m,j))
    end do
    xv = 0.5_dp
    upper = 1.0_dp
    bta = 0.99995_dp
    if (present(beta)) bta = beta
    tol = 1.0e-6_dp
    if (present(eps)) tol = eps
    call fnb_core(xa, -ya, rhs, xv, upper, bta, tol, workcoef, dual, nit, info)
    allocate(result%coefficients(p), result%residuals(n), result%dual(n))
    result%coefficients = -workcoef
    result%residuals = y - matmul(x, result%coefficients)
    result%dual = dual(1:n)
    result%tau = tau
    result%iterations = nit(1)
    result%corrector_steps = nit(2)
    result%info = info
  end subroutine rq_fit_lasso

  subroutine rq_fit_scad(x, y, tau, lambda, result, alpha, beta, eps, max_outer)
    real(dp), intent(in) :: x(:,:), y(:), tau
    real(dp), intent(in) :: lambda(:)
    type(rq_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, beta, eps
    integer, intent(in), optional :: max_outer
    real(dp) :: a, bta, tol
    integer :: maxit, p, n, m, j, k, outer, nit(3), info
    integer, allocatable :: keep(:)
    real(dp), allocatable :: xa(:,:), ya(:), rhs(:), vrhs(:), xv(:), upper(:)
    real(dp), allocatable :: coef(:), oldcoef(:), dpen(:), dual(:)
    type(rq_result) :: init

    n = size(x,1)
    p = size(x,2)
    if (size(lambda) /= p .or. any(lambda < 0.0_dp)) then
      result%info = -1
      return
    end if
    a = 3.2_dp
    if (present(alpha)) a = alpha
    bta = 0.9995_dp
    if (present(beta)) bta = beta
    tol = 1.0e-6_dp
    if (present(eps)) tol = eps
    maxit = 100
    if (present(max_outer)) maxit = max_outer

    m = count(lambda > 0.0_dp)
    allocate(keep(m))
    k = 0
    do j = 1, p
      if (lambda(j) > 0.0_dp) then
        k = k + 1
        keep(k) = j
      end if
    end do
    allocate(xa(n+m,p), ya(n+m), rhs(p), vrhs(p), xv(n+m), upper(n+m))
    allocate(coef(p), oldcoef(p), dpen(p), dual(n+m))
    xa(1:n,:) = x
    ya(1:n) = y
    xa(n+1:n+m,:) = 0.0_dp
    ya(n+1:n+m) = 0.0_dp
    do k = 1, m
      xa(n+k,keep(k)) = lambda(keep(k))
    end do
    do j = 1, p
      rhs(j) = (1.0_dp - tau) * sum(x(:,j)) + sum(xa(n+1:n+m,j))
    end do

    call rq_fit_fnb(x, y, tau, init, eps=tol)
    if (init%info /= 0) then
      result%info = init%info
      return
    end if
    coef = init%coefficients
    oldcoef = huge(1.0_dp)
    do outer = 1, maxit
      if (sum(abs(coef - oldcoef)) <= tol) exit
      oldcoef = coef
      dpen = scad_derivative(coef, lambda, a)
      dpen(1) = 0.0_dp
      vrhs = rhs - dpen * sign(1.0_dp, coef)
      xv(1:n) = 1.0_dp - tau
      if (m > 0) xv(n+1:n+m) = 0.5_dp
      upper = 1.0_dp
      call fnb_core(xa, -ya, vrhs, xv, upper, bta, tol, coef, dual, nit, info)
      coef = -coef
      if (info /= 0) exit
    end do
    allocate(result%coefficients(p), result%residuals(n), result%dual(n))
    result%coefficients = coef
    result%residuals = y - matmul(x, coef)
    result%dual = dual(1:n)
    result%tau = tau
    result%iterations = outer - 1
    result%corrector_steps = nit(2)
    result%info = info
  contains
    pure function scad_derivative(b, lam, aa) result(d)
      real(dp), intent(in) :: b(:), lam(:), aa
      real(dp) :: d(size(b))
      integer :: ii
      do ii = 1, size(b)
        if (abs(b(ii)) <= lam(ii)) then
          d(ii) = lam(ii)
        else if (abs(b(ii)) <= aa * lam(ii)) then
          d(ii) = (aa * lam(ii) - abs(b(ii))) / (aa - 1.0_dp)
        else
          d(ii) = 0.0_dp
        end if
      end do
    end function scad_derivative
  end subroutine rq_fit_scad

  subroutine fnb_core(a, c, b, xvar, upper, beta, eps, coef_work, dual_obs, nit, info)
    real(dp), intent(in) :: a(:,:), c(:), b(:), upper(:), beta, eps
    real(dp), intent(inout) :: xvar(:)
    real(dp), intent(out) :: coef_work(:), dual_obs(:)
    integer, intent(out) :: nit(3), info
    integer, parameter :: maxit = 500
    integer :: n, p, i
    real(dp) :: gap, deltap, deltad, mu, g, bigv, dxdz, dsdw
    real(dp), allocatable :: s(:), z(:), w(:), d(:), dx(:), ds(:), dz(:), dw(:), dr(:)
    real(dp), allocatable :: dual(:), dy(:), rhs(:), uwork(:)

    n = size(a,1)
    p = size(a,2)
    allocate(s(n), z(n), w(n), d(n), dx(n), ds(n), dz(n), dw(n), dr(n))
    allocate(dual(p), dy(p), rhs(p), uwork(n))
    nit = [0, 0, n]
    info = 0
    bigv = 1.0e20_dp

    dual = matmul(transpose(a), c)
    d = 1.0_dp
    call stepy(a, d, dual, info)
    if (info /= 0) return
    s = c - matmul(a, dual)
    do i = 1, n
      if (abs(s(i)) < eps) then
        z(i) = max(s(i), 0.0_dp) + eps
        w(i) = max(-s(i), 0.0_dp) + eps
      else
        z(i) = max(s(i), 0.0_dp)
        w(i) = max(-s(i), 0.0_dp)
      end if
      s(i) = upper(i) - xvar(i)
    end do
    gap = dot_product(z, xvar) + dot_product(w, s)

    do while (gap > eps .and. nit(1) < maxit)
      nit(1) = nit(1) + 1
      do i = 1, n
        d(i) = 1.0_dp / (z(i) / xvar(i) + w(i) / s(i))
        ds(i) = z(i) - w(i)
        dz(i) = d(i) * ds(i)
      end do
      dy = b - matmul(transpose(a), xvar) + matmul(transpose(a), dz)
      rhs = dy
      call stepy(a, d, dy, info)
      if (info /= 0) return
      ds = matmul(a, dy) - ds
      deltap = bigv
      deltad = bigv
      do i = 1, n
        dx(i) = d(i) * ds(i)
        ds(i) = -dx(i)
        dz(i) = -z(i) * (dx(i) / xvar(i) + 1.0_dp)
        dw(i) = -w(i) * (ds(i) / s(i) + 1.0_dp)
        if (dx(i) < 0.0_dp) deltap = min(deltap, -xvar(i) / dx(i))
        if (ds(i) < 0.0_dp) deltap = min(deltap, -s(i) / ds(i))
        if (dz(i) < 0.0_dp) deltad = min(deltad, -z(i) / dz(i))
        if (dw(i) < 0.0_dp) deltad = min(deltad, -w(i) / dw(i))
      end do
      deltap = min(beta * deltap, 1.0_dp)
      deltad = min(beta * deltad, 1.0_dp)
      if (min(deltap, deltad) < 1.0_dp) then
        nit(2) = nit(2) + 1
        mu = dot_product(xvar, z) + dot_product(s, w)
        g = mu + deltap * dot_product(dx, z) + deltad * dot_product(dz, xvar)
        g = g + deltap * deltad * dot_product(dz, dx)
        g = g + deltap * dot_product(ds, w) + deltad * dot_product(dw, s)
        g = g + deltap * deltad * dot_product(ds, dw)
        mu = mu * (g / mu)**3 / real(2*n, dp)
        do i = 1, n
          dr(i) = d(i) * (mu * (1.0_dp / s(i) - 1.0_dp / xvar(i)) &
            + dx(i) * dz(i) / xvar(i) - ds(i) * dw(i) / s(i))
        end do
        dy = rhs + matmul(transpose(a), dr)
        call stepy_solved(a, d, dy, info)
        if (info /= 0) return
        uwork = matmul(a, dy)
        deltap = bigv
        deltad = bigv
        do i = 1, n
          dxdz = dx(i) * dz(i)
          dsdw = ds(i) * dw(i)
          dx(i) = d(i) * (uwork(i) - z(i) + w(i)) - dr(i)
          ds(i) = -dx(i)
          dz(i) = -z(i) + (mu - z(i) * dx(i) - dxdz) / xvar(i)
          dw(i) = -w(i) + (mu - w(i) * ds(i) - dsdw) / s(i)
          if (dx(i) < 0.0_dp) deltap = min(deltap, -xvar(i) / dx(i))
          if (ds(i) < 0.0_dp) deltap = min(deltap, -s(i) / ds(i))
          if (dz(i) < 0.0_dp) deltad = min(deltad, -z(i) / dz(i))
          if (dw(i) < 0.0_dp) deltad = min(deltad, -w(i) / dw(i))
        end do
        deltap = min(beta * deltap, 1.0_dp)
        deltad = min(beta * deltad, 1.0_dp)
      end if
      xvar = xvar + deltap * dx
      s = s + deltap * ds
      dual = dual + deltad * dy
      z = z + deltad * dz
      w = w + deltad * dw
      gap = dot_product(z, xvar) + dot_product(w, s)
    end do
    dual_obs = z - w
    coef_work = dual
  end subroutine fnb_core

  subroutine stepy(a, d, b, info)
    real(dp), intent(in) :: a(:,:), d(:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ada(:,:)
    integer :: i, p

    p = size(a,2)
    allocate(ada(p,p))
    ada = 0.0_dp
    do i = 1, size(a,1)
      ada = ada + d(i) * outer(a(i,:), a(i,:))
    end do
    call spd_solve(ada, b, info)
  end subroutine stepy

  subroutine stepy_solved(a, d, b, info)
    real(dp), intent(in) :: a(:,:), d(:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info
    call stepy(a, d, b, info)
  end subroutine stepy_solved

  subroutine fnc_core(a1, c1, a2, c2, b, x1, x2, beta, eps, coef_work, dual_obs, nit, info)
    real(dp), intent(in) :: a1(:,:), c1(:), a2(:,:), c2(:), b(:), beta, eps
    real(dp), intent(inout) :: x1(:), x2(:)
    real(dp), intent(out) :: coef_work(:), dual_obs(:)
    integer, intent(out) :: nit(3), info
    integer, parameter :: maxit = 500
    integer :: n1, n2, p, i
    real(dp) :: gap, deltap, deltad, mu, g, bigv, dxdz1, dxdz2, dsdw
    real(dp), allocatable :: s(:), z1(:), z2(:), w(:), d1(:), d2(:)
    real(dp), allocatable :: dx1(:), dx2(:), ds(:), dz1(:), dz2(:), dw(:)
    real(dp), allocatable :: dr1(:), dr2(:), r2(:), dual(:), dy(:), rhs(:), u1(:), u2(:)

    n1 = size(a1,1)
    n2 = size(a2,1)
    p = size(a1,2)
    allocate(s(n1), z1(n1), z2(n2), w(n1), d1(n1), d2(n2))
    allocate(dx1(n1), dx2(n2), ds(n1), dz1(n1), dz2(n2), dw(n1))
    allocate(dr1(n1), dr2(n2), r2(n2), dual(p), dy(p), rhs(p), u1(n1), u2(n2))
    nit = [0, 0, n1]
    info = 0
    bigv = 1.0e20_dp

    dual = matmul(transpose(a1), c1)
    d1 = 1.0_dp
    d2 = 0.0_dp
    z2 = 1.0_dp
    call stepy2(a1, d1, a2, d2, dual, info)
    if (info /= 0) return
    s = c1 - matmul(a1, dual)
    do i = 1, n1
      if (abs(s(i)) < eps) then
        z1(i) = max(s(i), 0.0_dp) + eps
        w(i) = max(-s(i), 0.0_dp) + eps
      else
        z1(i) = max(s(i), 0.0_dp)
        w(i) = max(-s(i), 0.0_dp)
      end if
      s(i) = 1.0_dp - x1(i)
    end do
    gap = dot_product(z1,x1) + dot_product(z2,x2) + dot_product(w,s)

    do while (gap > eps .and. nit(1) < maxit)
      nit(1) = nit(1) + 1
      r2 = c2 - matmul(a2, dual)
      dy = b - matmul(transpose(a1),x1) - matmul(transpose(a2),x2)
      do i = 1, n1
        d1(i) = 1.0_dp / (z1(i)/x1(i) + w(i)/s(i))
        ds(i) = z1(i) - w(i)
        dz1(i) = d1(i) * ds(i)
      end do
      do i = 1, n2
        d2(i) = x2(i) / z2(i)
        dz2(i) = d2(i) * r2(i)
      end do
      dy = dy + matmul(transpose(a1),dz1) + matmul(transpose(a2),dz2)
      rhs = dy
      call stepy2(a1,d1,a2,d2,dy,info)
      if (info /= 0) return
      ds = matmul(a1,dy) - ds
      deltap = bigv
      deltad = bigv
      do i = 1, n1
        dx1(i) = d1(i) * ds(i)
        ds(i) = -dx1(i)
        dz1(i) = -z1(i) * (dx1(i)/x1(i) + 1.0_dp)
        dw(i) = -w(i) * (ds(i)/s(i) + 1.0_dp)
        if (dx1(i) < 0.0_dp) deltap = min(deltap,-x1(i)/dx1(i))
        if (ds(i) < 0.0_dp) deltap = min(deltap,-s(i)/ds(i))
        if (dz1(i) < 0.0_dp) deltad = min(deltad,-z1(i)/dz1(i))
        if (dw(i) < 0.0_dp) deltad = min(deltad,-w(i)/dw(i))
      end do
      dx2 = r2 - matmul(a2,dy)
      do i = 1, n2
        dx2(i) = d2(i) * dx2(i)
        dz2(i) = -z2(i) * (dx2(i)/x2(i) + 1.0_dp)
        if (dx2(i) < 0.0_dp) deltap = min(deltap,-x2(i)/dx2(i))
        if (dz2(i) < 0.0_dp) deltad = min(deltad,-z2(i)/dz2(i))
      end do
      deltap = min(beta*deltap,1.0_dp)
      deltad = min(beta*deltad,1.0_dp)
      if (min(deltap,deltad) < 1.0_dp) then
        nit(2) = nit(2) + 1
        mu = dot_product(x1,z1) + dot_product(x2,z2) + dot_product(s,w)
        g = mu + deltap*dot_product(dx1,z1) + deltad*dot_product(dz1,x1)
        g = g + deltap*deltad*dot_product(dz1,dx1)
        g = g + deltap*dot_product(dx2,z2) + deltad*dot_product(dz2,x2)
        g = g + deltap*deltad*dot_product(dz2,dx2)
        g = g + deltap*dot_product(ds,w) + deltad*dot_product(dw,s)
        g = g + deltap*deltad*dot_product(ds,dw)
        mu = mu * (g/mu)**3 / real(2*n1+n2,dp)
        do i = 1, n1
          dsdw = ds(i)*dw(i)
          dr1(i) = d1(i) * (mu*(1.0_dp/s(i)-1.0_dp/x1(i)) &
            + dx1(i)*dz1(i)/x1(i)-dsdw/s(i))
        end do
        do i = 1, n2
          dr2(i) = d2(i) * (dx2(i)*dz2(i)/x2(i)-mu/x2(i))
        end do
        dy = rhs + matmul(transpose(a1),dr1) + matmul(transpose(a2),dr2)
        call stepy2(a1,d1,a2,d2,dy,info)
        if (info /= 0) return
        u1 = matmul(a1,dy)
        u2 = matmul(a2,dy)
        deltap = bigv
        deltad = bigv
        do i = 1, n1
          dsdw = ds(i)*dw(i)
          dxdz1 = dx1(i)*dz1(i)
          dx1(i) = d1(i)*(u1(i)-z1(i)+w(i))-dr1(i)
          ds(i) = -dx1(i)
          dz1(i) = -z1(i)+(mu-z1(i)*dx1(i)-dxdz1)/x1(i)
          dw(i) = -w(i)+(mu-w(i)*ds(i)-dsdw)/s(i)
          if (dx1(i) < 0.0_dp) deltap = min(deltap,-x1(i)/dx1(i))
          if (ds(i) < 0.0_dp) deltap = min(deltap,-s(i)/ds(i))
          if (dz1(i) < 0.0_dp) deltad = min(deltad,-z1(i)/dz1(i))
          if (dw(i) < 0.0_dp) deltad = min(deltad,-w(i)/dw(i))
        end do
        do i = 1, n2
          dxdz2 = dx2(i)*dz2(i)
          dx2(i) = d2(i)*(u2(i)-r2(i))-dr2(i)
          dz2(i) = -z2(i)+(mu-z2(i)*dx2(i)-dxdz2)/x2(i)
          if (dx2(i) < 0.0_dp) deltap = min(deltap,-x2(i)/dx2(i))
          if (dz2(i) < 0.0_dp) deltad = min(deltad,-z2(i)/dz2(i))
        end do
        deltap = min(beta*deltap,1.0_dp)
        deltad = min(beta*deltad,1.0_dp)
      end if
      x1 = x1 + deltap*dx1
      x2 = x2 + deltap*dx2
      s = s + deltap*ds
      dual = dual + deltad*dy
      z1 = z1 + deltad*dz1
      z2 = z2 + deltad*dz2
      w = w + deltad*dw
      gap = dot_product(z1,x1)+dot_product(z2,x2)+dot_product(w,s)
    end do
    dual_obs = z1 - w
    coef_work = dual
  end subroutine fnc_core

  subroutine stepy2(a1,d1,a2,d2,b,info)
    real(dp), intent(in) :: a1(:,:), d1(:), a2(:,:), d2(:)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ada(:,:)
    integer :: i, p

    p = size(a1,2)
    allocate(ada(p,p))
    ada = 0.0_dp
    do i = 1, size(a1,1)
      ada = ada + d1(i)*outer(a1(i,:),a1(i,:))
    end do
    do i = 1, size(a2,1)
      ada = ada + d2(i)*outer(a2(i,:),a2(i,:))
    end do
    call spd_solve(ada,b,info)
  end subroutine stepy2

  pure function outer(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i = 1, size(x)
      a(i,:) = x(i)*y
    end do
  end function outer

end module quantreg_dense
