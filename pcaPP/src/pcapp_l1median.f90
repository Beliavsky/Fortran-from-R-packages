module pcapp_l1median
  use pcapp_kinds, only: dp
  use pcapp_types, only: median_result
  use pcapp_stats, only: median_vec, spatial_objective, spatial_gradient, spatial_hessian
  implicit none
  private

  public :: l1median, l1median_bfgs, l1median_cg, l1median_hocr
  public :: l1median_nlm, l1median_nm, l1median_vazh

contains

  pure subroutine default_start(x, start)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix.
    real(dp), intent(out) :: start(size(x, 2)) !! Coordinatewise-median starting point.
    integer :: j

    do j = 1, size(x, 2)
      start(j) = median_vec(x(:, j))
    end do
  end subroutine default_start

  pure subroutine choose_start(x, m_init, start)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix.
    real(dp), intent(in), optional :: m_init(:) !! Optional initial spatial-median estimate of length ncol(x).
    real(dp), intent(out) :: start(size(x, 2)) !! Validated starting point used by an optimizer.

    if (present(m_init)) then
      if (size(m_init) == size(x, 2)) then
        start = m_init
      else
        call default_start(x, start)
      end if
    else
      call default_start(x, start)
    end if
  end subroutine choose_start

  pure subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:) !! Square coefficient matrix for a small dense linear system.
    real(dp), intent(in) :: b(:) !! Right-hand-side vector matching the matrix order.
    real(dp), intent(out) :: x(size(b)) !! Solution vector when ok is true.
    logical, intent(out) :: ok !! True when Gaussian elimination found a nonsingular pivot sequence.
    real(dp), allocatable :: m(:,:), rhs(:), row_tmp(:)
    real(dp) :: factor, pivot
    integer :: i, j, k, best, n

    n = size(b)
    allocate(m(n, n), rhs(n), row_tmp(n))
    m = a
    rhs = b
    ok = .true.
    do k = 1, n - 1
      best = k
      do i = k + 1, n
        if (abs(m(i, k)) > abs(m(best, k))) best = i
      end do
      if (abs(m(best, k)) <= tiny(1.0_dp)) then
        ok = .false.
        x = 0.0_dp
        return
      end if
      if (best /= k) then
        row_tmp = m(k, :)
        m(k, :) = m(best, :)
        m(best, :) = row_tmp
        pivot = rhs(k)
        rhs(k) = rhs(best)
        rhs(best) = pivot
      end if
      do i = k + 1, n
        factor = m(i, k) / m(k, k)
        m(i, k:n) = m(i, k:n) - factor * m(k, k:n)
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do
    if (n > 0 .and. abs(m(n, n)) <= tiny(1.0_dp)) then
      ok = .false.
      x = 0.0_dp
      return
    end if
    x = 0.0_dp
    do i = n, 1, -1
      pivot = rhs(i)
      do j = i + 1, n
        pivot = pivot - m(i, j) * x(j)
      end do
      x(i) = pivot / m(i, i)
    end do
  end subroutine solve_linear

  pure function l1median(x, maxstep, ittol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix whose L1/spatial median is requested.
    integer, intent(in), optional :: maxstep !! Maximum optimization iterations; pcaPP defaults to 200.
    real(dp), intent(in), optional :: ittol !! Convergence tolerance; pcaPP defaults to 1e-8.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    integer :: mi
    real(dp) :: tt

    mi = 200
    if (present(maxstep)) mi = maxstep
    tt = 1.0e-8_dp
    if (present(ittol)) tt = ittol
    res = l1median_nlm(x, mi, tt, m_init=m_init)
  end function l1median

  pure function l1median_vazh(x, maxit, tol, zero_tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix for the Vardi-Zhang spatial median.
    integer, intent(in), optional :: maxit !! Maximum Weiszfeld/Vardi-Zhang iterations; default 200.
    real(dp), intent(in), optional :: tol !! Relative L1 change tolerance; default 1e-8.
    real(dp), intent(in), optional :: zero_tol !! Coincidence threshold for observations; default 1e-15.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: m(:), old(:), tpoint(:), rvec(:), dist(:)
    real(dp) :: tolerance, ztol, sumw, rnorm, eta, diff, scale
    integer :: i, iter, max_iter, nz

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    ztol = 1.0e-15_dp
    if (present(zero_tol)) ztol = zero_tol
    allocate(m(size(x, 2)), old(size(x, 2)), tpoint(size(x, 2)), rvec(size(x, 2)))
    allocate(dist(size(x, 1)))
    call choose_start(x, m_init, m)
    res%code = 1

    do iter = 1, max_iter
      old = m
      nz = 0
      do i = 1, size(x, 1)
        dist(i) = sqrt(sum((x(i, :) - m)**2))
        if (dist(i) <= ztol) nz = nz + 1
      end do
      if (2 * nz > size(x, 1)) then
        res%code = 3
        exit
      end if

      tpoint = 0.0_dp
      rvec = 0.0_dp
      sumw = 0.0_dp
      do i = 1, size(x, 1)
        if (dist(i) <= ztol) cycle
        tpoint = tpoint + x(i, :) / dist(i)
        rvec = rvec + (x(i, :) - m) / dist(i)
        sumw = sumw + 1.0_dp / dist(i)
      end do
      if (sumw <= tiny(1.0_dp)) then
        res%code = 3
        exit
      end if
      tpoint = tpoint / sumw
      if (nz > 0) then
        rnorm = sqrt(sum(rvec**2))
        if (rnorm <= real(nz, dp)) then
          m = old
        else
          eta = real(nz, dp) / rnorm
          m = (1.0_dp - eta) * tpoint + eta * old
        end if
      else
        m = tpoint
      end if

      diff = sum(abs(m - old))
      scale = max(sum(abs(m)), 1.0_dp)
      if (diff < tolerance * scale) then
        res%code = 0
        exit
      end if
    end do

    res%par = m
    res%value = spatial_objective(x, m)
    res%iterations = min(iter, max_iter)
  end function l1median_vazh

  pure function l1median_hocr(x, maxit, tol, zero_tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable data matrix for the HoCr spatial-median algorithm.
    integer, intent(in), optional :: maxit !! Maximum iterations; default 200.
    real(dp), intent(in), optional :: tol !! Euclidean step tolerance; default 1e-8.
    real(dp), intent(in), optional :: zero_tol !! Coincidence threshold; default 1e-15.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: m(:), old(:), delta(:), dist(:), weights(:)
    real(dp) :: tolerance, ztol, obj_old, obj_new, sumw, nd
    integer :: i, iter, max_iter, nzero, half, max_half

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    ztol = 1.0e-15_dp
    if (present(zero_tol)) ztol = zero_tol
    allocate(m(size(x, 2)), old(size(x, 2)), delta(size(x, 2)))
    allocate(dist(size(x, 1)), weights(size(x, 1)))
    call choose_start(x, m_init, m)
    obj_old = spatial_objective(x, m)
    res%code = 1

    do iter = 1, max_iter
      nzero = 0
      weights = 0.0_dp
      do i = 1, size(x, 1)
        dist(i) = sqrt(sum((x(i, :) - m)**2))
        if (dist(i) <= ztol) then
          nzero = nzero + 1
        else
          weights(i) = 1.0_dp / dist(i)
        end if
      end do
      if (2 * nzero > size(x, 1)) then
        res%code = 3
        exit
      end if
      sumw = sum(weights)
      if (sumw <= tiny(1.0_dp)) then
        res%code = 3
        exit
      end if
      delta = 0.0_dp
      do i = 1, size(x, 1)
        delta = delta + weights(i) * (x(i, :) - m)
      end do
      delta = delta / sumw
      nd = sqrt(sum(delta**2))
      old = m
      m = old + delta
      obj_new = spatial_objective(x, m)
      if (nd < tolerance) then
        res%code = 0
        obj_old = obj_new
        exit
      end if
      max_half = max(1, ceiling(log(max(nd / tolerance, 1.0_dp)) / log(2.0_dp)) + 2)
      do half = 1, max_half
        if (obj_new < obj_old) exit
        delta = 0.5_dp * delta
        m = old + delta
        obj_new = spatial_objective(x, m)
      end do
      if (obj_new >= obj_old) then
        m = old
        res%code = 2
        exit
      end if
      obj_old = obj_new
    end do

    res%par = m
    res%value = spatial_objective(x, m)
    res%iterations = min(iter, max_iter)
  end function l1median_hocr

  pure function l1median_nlm(x, maxit, tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for damped-Newton spatial-median fitting.
    integer, intent(in), optional :: maxit !! Maximum Newton iterations; default 200.
    real(dp), intent(in), optional :: tol !! Gradient/step convergence tolerance; default 1e-8.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: m(:), grad(:), hess(:,:), step(:), trial(:)
    real(dp) :: tolerance, obj, obj_trial, alpha, gnorm, ridge
    integer :: i, iter, max_iter
    logical :: ok

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    allocate(m(size(x, 2)), grad(size(x, 2)), hess(size(x, 2), size(x, 2)))
    allocate(step(size(x, 2)), trial(size(x, 2)))
    call choose_start(x, m_init, m)
    obj = spatial_objective(x, m)
    res%code = 1

    do iter = 1, max_iter
      call spatial_gradient(x, m, sqrt(epsilon(1.0_dp)), grad)
      gnorm = sqrt(sum(grad**2))
      if (gnorm <= tolerance) then
        res%code = 0
        exit
      end if
      call spatial_hessian(x, m, sqrt(epsilon(1.0_dp)), hess)
      ridge = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(hess)))
      do i = 1, size(m)
        hess(i, i) = hess(i, i) + ridge
      end do
      call solve_linear(hess, -grad, step, ok)
      if (.not. ok .or. dot_product(step, grad) >= 0.0_dp) step = -grad
      alpha = 1.0_dp
      do i = 1, 40
        trial = m + alpha * step
        obj_trial = spatial_objective(x, trial)
        if (obj_trial <= obj + 1.0e-4_dp * alpha * dot_product(grad, step)) exit
        alpha = 0.5_dp * alpha
      end do
      if (alpha * sqrt(sum(step**2)) <= tolerance * max(1.0_dp, sqrt(sum(m**2)))) then
        m = trial
        obj = obj_trial
        res%code = 0
        exit
      end if
      m = trial
      obj = obj_trial
    end do

    res%par = m
    res%value = obj
    res%iterations = min(iter, max_iter)
  end function l1median_nlm

  pure function l1median_bfgs(x, maxit, tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for BFGS spatial-median fitting.
    integer, intent(in), optional :: maxit !! Maximum BFGS iterations; default 200.
    real(dp), intent(in), optional :: tol !! Gradient/step convergence tolerance; default 1e-8.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: m(:), grad(:), grad_new(:), hinv(:,:), direction(:), trial(:), s(:), y(:)
    real(dp), allocatable :: identity(:,:), v(:,:)
    real(dp) :: tolerance, obj, obj_trial, alpha, ys, rho
    integer :: i, iter, max_iter, p

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    p = size(x, 2)
    allocate(m(p), grad(p), grad_new(p), hinv(p, p), direction(p), trial(p), s(p), y(p))
    allocate(identity(p, p), v(p, p))
    call choose_start(x, m_init, m)
    identity = 0.0_dp
    do i = 1, p
      identity(i, i) = 1.0_dp
    end do
    hinv = identity
    obj = spatial_objective(x, m)
    call spatial_gradient(x, m, sqrt(epsilon(1.0_dp)), grad)
    res%code = 1

    do iter = 1, max_iter
      if (sqrt(sum(grad**2)) <= tolerance) then
        res%code = 0
        exit
      end if
      direction = -matmul(hinv, grad)
      if (dot_product(direction, grad) >= 0.0_dp) direction = -grad
      alpha = 1.0_dp
      do i = 1, 40
        trial = m + alpha * direction
        obj_trial = spatial_objective(x, trial)
        if (obj_trial <= obj + 1.0e-4_dp * alpha * dot_product(grad, direction)) exit
        alpha = 0.5_dp * alpha
      end do
      s = trial - m
      call spatial_gradient(x, trial, sqrt(epsilon(1.0_dp)), grad_new)
      y = grad_new - grad
      ys = dot_product(y, s)
      if (ys > sqrt(epsilon(1.0_dp)) * sqrt(sum(y**2) * sum(s**2))) then
        rho = 1.0_dp / ys
        v = identity - rho * spread(s, 2, p) * spread(y, 1, p)
        hinv = matmul(matmul(v, hinv), transpose(v)) + rho * spread(s, 2, p) * spread(s, 1, p)
      else
        hinv = identity
      end if
      m = trial
      grad = grad_new
      obj = obj_trial
      if (sqrt(sum(s**2)) <= tolerance * max(1.0_dp, sqrt(sum(m**2)))) then
        res%code = 0
        exit
      end if
    end do

    res%par = m
    res%value = obj
    res%iterations = min(iter, max_iter)
  end function l1median_bfgs

  pure function l1median_cg(x, maxit, tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for conjugate-gradient spatial-median fitting.
    integer, intent(in), optional :: maxit !! Maximum nonlinear-CG iterations; default 200.
    real(dp), intent(in), optional :: tol !! Gradient/step convergence tolerance; default 1e-8.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: m(:), grad(:), grad_new(:), direction(:), trial(:), step(:)
    real(dp) :: tolerance, obj, obj_trial, alpha, beta, denom
    integer :: i, iter, max_iter, p

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    p = size(x, 2)
    allocate(m(p), grad(p), grad_new(p), direction(p), trial(p), step(p))
    call choose_start(x, m_init, m)
    obj = spatial_objective(x, m)
    call spatial_gradient(x, m, sqrt(epsilon(1.0_dp)), grad)
    direction = -grad
    res%code = 1

    do iter = 1, max_iter
      if (sqrt(sum(grad**2)) <= tolerance) then
        res%code = 0
        exit
      end if
      if (dot_product(direction, grad) >= 0.0_dp) direction = -grad
      alpha = 1.0_dp
      do i = 1, 40
        trial = m + alpha * direction
        obj_trial = spatial_objective(x, trial)
        if (obj_trial <= obj + 1.0e-4_dp * alpha * dot_product(grad, direction)) exit
        alpha = 0.5_dp * alpha
      end do
      step = trial - m
      call spatial_gradient(x, trial, sqrt(epsilon(1.0_dp)), grad_new)
      denom = max(dot_product(grad, grad), tiny(1.0_dp))
      beta = max(0.0_dp, dot_product(grad_new, grad_new - grad) / denom)
      direction = -grad_new + beta * direction
      m = trial
      grad = grad_new
      obj = obj_trial
      if (sqrt(sum(step**2)) <= tolerance * max(1.0_dp, sqrt(sum(m**2)))) then
        res%code = 0
        exit
      end if
    end do

    res%par = m
    res%value = obj
    res%iterations = min(iter, max_iter)
  end function l1median_cg

  pure subroutine sort_simplex(simplex, f)
    real(dp), intent(inout) :: simplex(:,:) !! Nelder-Mead simplex vertices stored by column.
    real(dp), intent(inout) :: f(:) !! Objective values associated with simplex columns.
    real(dp), allocatable :: temp(:)
    real(dp) :: tf
    integer :: i, j, best

    allocate(temp(size(simplex, 1)))
    do i = 1, size(f) - 1
      best = i
      do j = i + 1, size(f)
        if (f(j) < f(best)) best = j
      end do
      if (best /= i) then
        tf = f(i)
        f(i) = f(best)
        f(best) = tf
        temp = simplex(:, i)
        simplex(:, i) = simplex(:, best)
        simplex(:, best) = temp
      end if
    end do
  end subroutine sort_simplex

  pure function l1median_nm(x, maxit, tol, m_init) result(res)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for Nelder-Mead spatial-median fitting.
    integer, intent(in), optional :: maxit !! Maximum simplex iterations; default 200.
    real(dp), intent(in), optional :: tol !! Simplex convergence tolerance; default 1e-8.
    real(dp), intent(in), optional :: m_init(:) !! Optional starting point; default is coordinatewise medians.
    type(median_result) :: res
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: tolerance, fr, fe, fc, delta
    integer :: iter, j, max_iter, p

    max_iter = 200
    if (present(maxit)) max_iter = maxit
    tolerance = 1.0e-8_dp
    if (present(tol)) tolerance = tol
    p = size(x, 2)
    allocate(simplex(p, p + 1), f(p + 1), centroid(p), xr(p), xe(p), xc(p))
    call choose_start(x, m_init, simplex(:, 1))
    do j = 1, p
      simplex(:, j + 1) = simplex(:, 1)
      delta = 0.05_dp * max(1.0_dp, abs(simplex(j, 1)))
      simplex(j, j + 1) = simplex(j, j + 1) + delta
    end do
    do j = 1, p + 1
      f(j) = spatial_objective(x, simplex(:, j))
    end do
    res%code = 1

    do iter = 1, max_iter
      call sort_simplex(simplex, f)
      if (maxval(abs(f - f(1))) <= tolerance * max(1.0_dp, abs(f(1)))) then
        res%code = 0
        exit
      end if
      centroid = sum(simplex(:, 1:p), dim=2) / real(p, dp)
      xr = centroid + (centroid - simplex(:, p + 1))
      fr = spatial_objective(x, xr)
      if (fr < f(1)) then
        xe = centroid + 2.0_dp * (xr - centroid)
        fe = spatial_objective(x, xe)
        if (fe < fr) then
          simplex(:, p + 1) = xe
          f(p + 1) = fe
        else
          simplex(:, p + 1) = xr
          f(p + 1) = fr
        end if
      else if (fr < f(p)) then
        simplex(:, p + 1) = xr
        f(p + 1) = fr
      else
        if (fr < f(p + 1)) then
          xc = centroid + 0.5_dp * (xr - centroid)
        else
          xc = centroid + 0.5_dp * (simplex(:, p + 1) - centroid)
        end if
        fc = spatial_objective(x, xc)
        if (fc < min(fr, f(p + 1))) then
          simplex(:, p + 1) = xc
          f(p + 1) = fc
        else
          do j = 2, p + 1
            simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
            f(j) = spatial_objective(x, simplex(:, j))
          end do
        end if
      end if
    end do

    call sort_simplex(simplex, f)
    res%par = simplex(:, 1)
    res%value = f(1)
    res%iterations = min(iter, max_iter)
  end function l1median_nm

end module pcapp_l1median
